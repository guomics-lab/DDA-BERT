"""Base Transformer models for working with mass spectra and peptides"""
import re
import einops
import pandas as pd
import numpy as np
from typing import Tuple, List

import torch
from torch import nn


class NLinearMemoryEfficient(nn.Module):
    """Memory efficient implementation of parallel linear layers"""
    def __init__(self, n_features: int, d_in: int, d_out: int):
        super().__init__()
        self.weight = nn.Parameter(torch.randn(n_features, d_in, d_out))
        self.bias = nn.Parameter(torch.zeros(n_features, d_out))
        
    def forward(self, x):
        # x: (B, n_features, d_in)
        return torch.einsum('bfi,fio->bfo', x, self.weight) + self.bias


class NumEmbeddings(nn.Module):
    def __init__(
            self,
            n_features: int,
            d_embedding: int,
            embedding_arch: list,
            d_feature: int,
    ) -> None:
        super().__init__()
        assert embedding_arch
        assert set(embedding_arch) <= {
            'linear',
            'shared_linear',
            'relu',
            'layernorm',
            'batchnorm',
        }

        layers: list[nn.Module] = []

        if embedding_arch[0] == 'linear':
            assert d_embedding is not None
            layers.append(
                NLinearMemoryEfficient(n_features, d_feature, d_embedding)
            )
        elif embedding_arch[0] == 'shared_linear':
            layers.append(
                nn.Linear(d_feature, d_embedding)
            )
        d_current = d_embedding

        for x in embedding_arch[1:]:
            layers.append(
                nn.ReLU()
                if x == 'relu'
                else NLinearMemoryEfficient(n_features, d_current, d_embedding)
                if x == 'linear'
                else nn.Linear(d_current, d_embedding)
                if x == 'shared_linear'
                else nn.LayerNorm([n_features, d_current])
                if x == 'layernorm'
                else nn.BatchNorm1d(n_features)
                if x == 'batchnorm'
                else nn.Identity()
            )
            if x in ['linear']:
                d_current = d_embedding
            assert not isinstance(layers[-1], nn.Identity)
        self.d_embedding = d_current
        self.layers = nn.Sequential(*layers)

    def forward(self, x):
        return self.layers(x)


class MultiScalePeakEmbedding(nn.Module):
    """Multi-scale sinusoidal embedding based on Voronov et. al."""

    def __init__(self, h_size: int, dropout: float = 0) -> None:
        super().__init__()
        self.h_size = h_size

        self.mlp = nn.Sequential(
            nn.Linear(h_size, h_size),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(h_size, h_size),
            nn.Dropout(dropout),
        )

        self.head = nn.Sequential(
            nn.Linear(h_size + 1, h_size),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(h_size, h_size),
            nn.Dropout(dropout),
        )

        freqs = 2 * np.pi / torch.logspace(-2, -3, int(h_size / 2), dtype=torch.float64)
        self.register_buffer("freqs", freqs)

    def forward(self, mz_values: Tensor, intensities: Tensor) -> Tensor:
        """Encode peaks."""
        x = self.encode_mass(mz_values)
        if x.dtype != mz_values.dtype:
            x = x.to(mz_values.dtype)
        
        x = self.mlp(x)
        x = torch.cat([x, intensities], axis=2)
        return self.head(x)

    def encode_mass(self, x: Tensor) -> Tensor:
        """Encode mz."""
        x = self.freqs[None, None, :] * x
        x = torch.cat([torch.sin(x), torch.cos(x)], axis=2)
        return x.float()
    

class MassEncoder(torch.nn.Module):
    """Encode mass values using sine and cosine waves."""

    def __init__(self, dim_model, min_wavelength=0.001, max_wavelength=10000):
        """Initialize the MassEncoder"""
        super().__init__()

        n_sin = int(dim_model / 2)
        n_cos = dim_model - n_sin

        if min_wavelength:
            base = min_wavelength / (2 * np.pi)
            scale = max_wavelength / min_wavelength
        else:
            base = 1
            scale = max_wavelength / (2 * np.pi)

        sin_term = base * scale ** (
            torch.arange(0, n_sin).float() / (n_sin - 1)
        )
        cos_term = base * scale ** (
            torch.arange(0, n_cos).float() / (n_cos - 1)
        )

        self.register_buffer("sin_term", sin_term)
        self.register_buffer("cos_term", cos_term)

    def forward(self, X):
        """Encode m/z values."""
        sin_mz = torch.sin(X / self.sin_term)
        cos_mz = torch.cos(X / self.cos_term)
        return torch.cat([sin_mz, cos_mz], axis=-1)


class PositionalEncoder(torch.nn.Module):
    """The positional encoder for sequences."""

    def __init__(self, dim_model, max_wavelength=10000):
        """Initialize the positional encoder."""
        super().__init__()

        n_sin = int(dim_model / 2)
        n_cos = dim_model - n_sin
        scale = max_wavelength / (2 * np.pi)

        sin_term = scale ** (torch.arange(0, n_sin).float() / (n_sin - 1))
        cos_term = scale ** (torch.arange(0, n_cos).float() / (n_cos - 1))
        self.register_buffer("sin_term", sin_term)
        self.register_buffer("cos_term", cos_term)

    def forward(self, X):
        """Encode positions in a sequence."""
        pos = torch.arange(X.shape[1]).type_as(self.sin_term)
        pos = einops.repeat(pos, "n -> b n", b=X.shape[0])
        sin_in = einops.repeat(pos, "b n -> b n f", f=len(self.sin_term))
        cos_in = einops.repeat(pos, "b n -> b n f", f=len(self.cos_term))

        sin_pos = torch.sin(sin_in / self.sin_term)
        cos_pos = torch.cos(cos_in / self.cos_term)
        encoded = torch.cat([sin_pos, cos_pos], axis=2)
        return encoded + X


class PeptideDecoder(nn.Module):
    """
    Updated decoder (removes the "replace masked tokens with mass" logic):

    - Input tokens are already a sequence with contiguous masking applied
      (i.e., masked positions are <mask> token_ids).
    - The decoder performs standard Transformer decoding over the full sequence,
      producing hidden states for every token position. The upstream
      MaskedLanguageModel can directly predict the original tokens.
    - No additional mass vectors or special mass encodings are injected for masked tokens.
    """

    def __init__(
        self,
        dim_model: int = 768,
        n_head: int = 8,
        dim_feedforward: int = 1024,
        n_layers: int = 1,
        dropout: float = 0.1,
        vocab: list | None = None,
        max_charge: int = 5,
        hidden_size: int = 50,
        id2aa=None,
    ):
        super().__init__()

        if vocab is None:
            vocab = []
        self.vocab = vocab
        self.vocab_size = len(vocab)

        self.dim_model = dim_model
        self.hidden_size = hidden_size

        self.pos_encoder = PositionalEncoder(self.dim_model)
        self.charge_encoder = nn.Embedding(max_charge, self.dim_model)

        # Pure semantic token embeddings (PAD=0, MASK uses mask_token_id in vocab)
        self.aa_encoder = nn.Embedding(
            self.vocab_size,
            dim_model,
            padding_idx=0,
        )

        # Still use MassEncoder + numerical feature embedding for precursors;
        # it does not participate in token-level replacement
        self.mass_encoder = MassEncoder(self.dim_model)

        layer = nn.TransformerDecoderLayer(
            d_model=self.dim_model,
            nhead=n_head,
            dim_feedforward=dim_feedforward,
            batch_first=True,
            dropout=dropout,
        )
        self.transformer_decoder = nn.TransformerDecoder(layer, num_layers=n_layers)

        embedding_arch = ["shared_linear", "batchnorm", "relu"]
        self.num_embeddings = NumEmbeddings(
            n_features=dim_model,
            d_embedding=dim_model,
            embedding_arch=embedding_arch,
            d_feature=2,
        )

    def forward(
        self,
        memory: torch.Tensor,                  # (B, 1+N, D)
        memory_key_padding_mask: torch.Tensor, # (B, 1+N) bool
        precursors: torch.Tensor,              # (B, 4)
        tokens: torch.Tensor,                  # (B, L) already contains contiguous <mask> tokens
        token_mask: torch.Tensor,              # (B, L) bool, True=masked (only used for selecting positions in the upstream loss)
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Returns:
            decoder_output: (B, 1+L, D)
                Position 0 is the precursor token; the remaining positions are the hidden states
                for each text token.
            mask_positions_tensor: (B, L)
                Identity mapping [0..L-1]. Indices are provided for all positions; the upstream
                module decides which positions participate in MLM based on token_mask.
            peptide_embedding: (B, D)
                decoder_output[:, 0, :], used for the DDA task.
        """
        target_dtype = memory.dtype
        device = memory.device

        # ===== 1. Precursor embedding (mass + charge + RT and other numeric features) =====
        precursors = precursors.to(target_dtype)
        masses = self.mass_encoder(precursors[:, None, [0]])          # (B,1,D)
        charges = self.charge_encoder(precursors[:, 1].int() - 1)     # (B,D)
        rt = self.num_embeddings(precursors[:, 2:])                   # (B,dim_model)
        precursors_emb = masses + charges[:, None, :] + rt[:, None, :]  # (B,1,D)

        # ===== 2. Token embedding (no mass replacement) =====
        tokens = tokens.long()
        aa_embed = self.aa_encoder(tokens)                            # (B,L,D)

        # Do not compress masked tokens; keep the full sequence so the upstream head
        # can take logits at any token_mask positions
        tgt_tokens_emb = aa_embed                                     # (B,L,D)

        # ===== 3. Concatenate precursor token + token sequence =====
        tgt = torch.cat([precursors_emb, tgt_tokens_emb], dim=1)      # (B,1+L,D)

        B, total_len, _ = tgt.size()

        # key_padding_mask: position 0 (precursor) is always valid; the rest depends on PAD
        is_padding = (tokens == 0)                                    # (B,L)
        precursor_pad = torch.zeros(B, 1, dtype=torch.bool, device=device)
        tgt_key_padding_mask = torch.cat([precursor_pad, is_padding], dim=1)  # (B,1+L)

        # ===== 4. Positional encoding =====
        tgt = self.pos_encoder(tgt)                                   # (B,1+L,D)

        # ===== 5. Transformer decoder (no autoregressive mask) =====
        tgt_mask = generate_no_mask(total_len).to(device).bool()

        decoder_output = self.transformer_decoder(
            tgt=tgt,
            memory=memory,
            tgt_mask=tgt_mask,
            tgt_key_padding_mask=tgt_key_padding_mask,
            memory_key_padding_mask=memory_key_padding_mask.bool(),
        )                                                              # (B,1+L,D)

        peptide_embedding = decoder_output[:, 0, :]                   # (B,D)

        # mask_positions_tensor: the index of each original token position in decoder_output
        # decoder_output[:, 1 + i, :] corresponds to tokens[:, i]
        L = tokens.size(1)
        pos = torch.arange(L, device=device, dtype=torch.long)        # [0..L-1]
        mask_positions_tensor = einops.repeat(pos, "l -> b l", b=B)   # (B,L)

        return decoder_output, mask_positions_tensor, peptide_embedding


def generate_no_mask(sz):
    """Generate a no mask for the sequence."""
    mask = torch.zeros(sz, sz).float()
    return mask


class MaskedLanguageModel(nn.Module):
    """Predicting origin token from masked input sequence"""

    def __init__(self, hidden, vocab_size):
        super().__init__()
        self.linear = nn.Linear(hidden, vocab_size)

    def forward(self, x):
        return self.linear(x)
