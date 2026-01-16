from __future__ import annotations

import torch
from torch import nn
from torch import Tensor
import torch.nn.functional as F
from typing import List

import yaml
import numpy as np
import glob
import os

from transformer.decoder import PeptideDecoder, MultiScalePeakEmbedding, MaskedLanguageModel


class DDABert(nn.Module):
    """The DDABert model."""

    def __init__(
            self,
            dim_model: int = 768,
            n_head: int = 16,
            dim_feedforward: int = 1024,
            n_layers: int = 9,
            dropout: float = 0,
            max_length: int = 50,
            vocab: list = [],
            max_charge: int = 10,
    ) -> None:
        super().__init__()
        self.dim_model = dim_model
        self.max_length = max_length
        self.vocab_size = len(vocab)

        # The latent representations for the spectrum
        self.latent_spectrum = nn.Parameter(torch.randn(1, 1, dim_model))

        # Encoder
        self.peak_encoder = MultiScalePeakEmbedding(dim_model, dropout=dropout)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim_model,
            nhead=n_head,
            dim_feedforward=dim_feedforward,
            batch_first=True,
            dropout=dropout,
        )
        self.encoder = nn.TransformerEncoder(
            encoder_layer,
            num_layers=n_layers,
        )

        # Decoder
        self.spectrum_sequence_encoder = PeptideDecoder(
            dim_model=dim_model,
            n_head=n_head,
            dim_feedforward=dim_feedforward,
            n_layers=n_layers,
            dropout=dropout,
            vocab=vocab,
            max_charge=max_charge,
            hidden_size=max_length
        )

        # Network heads for the DDA task
        self.psm_1 = nn.Linear(self.dim_model, 64)
        self.psm_2 = nn.Linear(64, 1)

        # Network head for the peptide masking (MLM) task
        self.mask_lm = MaskedLanguageModel(self.dim_model, self.vocab_size)

        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(p=dropout)

    def forward(
            self,
            spectra: Tensor,
            spectra_mask: Tensor,
            precursors: Tensor,
            tokens: Tensor,
            token_mask: Tensor,
    ) -> tuple:
        """Model forward pass with dtype handling."""

        # Ensure input dtypes are correct:
        # - spectra and precursors should already be in target_dtype (e.g., bfloat16)
        # - tokens must be long
        tokens = tokens.long()
        token_mask = token_mask.bool()

        # Convert spectra_mask to bool (avoid dtype-related issues)
        if spectra_mask.dtype != torch.bool:
            spectra_mask = spectra_mask.bool()

        # ============ 1. Encode spectrum ============
        spectra_encoded, spectra_mask_encoded = self._encoder(spectra, spectra_mask)
        spectrum_emb = spectra_encoded[:, 0, :]

        # ============ 2. Decode with peptide ============
        decoder_output, mask_positions, peptide_emb = self._decoder(
            spectra_encoded, spectra_mask_encoded, precursors, tokens, token_mask
        )

        # ============ 3. DDA prediction ============
        dda_pred = self.dropout(self.relu(self.psm_1(peptide_emb)))
        dda_pred = self.psm_2(dda_pred).squeeze(-1)

        # ============ 4. MLM prediction ============
        tokens_output = decoder_output[:, 1:, :]
        mask_pred = self.mask_lm(tokens_output)
        mask_pred = mask_pred.transpose(1, 2)

        return dda_pred, mask_pred, mask_positions

    def _encoder(self, spectra: Tensor, spectra_mask: Tensor) -> tuple[Tensor, Tensor]:
        """Spectrum encoder

        Args:
            spectra: (B, N, 2)
            spectra_mask: (B, N) - bool, True=padding

        Returns:
            spectra_encoded: (B, 1+N, D) - [latent_token, peak_1, .. ., peak_N]
            spectra_mask_encoded: (B, 1+N) - [False, mask_1, ..., mask_N]
        """
        # Ensure spectra_mask is bool
        spectra_mask = spectra_mask.bool()

        # Peak encoding
        spectra = self.peak_encoder(spectra[:, :, [0]], spectra[:, :, [1]])

        # Add latent spectrum token
        latent_spectra = self.latent_spectrum.expand(spectra.shape[0], -1, -1)
        spectra = torch.cat([latent_spectra, spectra], dim=1)

        # Add latent mask
        latent_mask = torch.zeros((spectra_mask.shape[0], 1), dtype=torch.bool, device=spectra_mask.device)
        spectra_mask = torch.cat([latent_mask, spectra_mask], dim=1)

        # Self-attention
        spectra = self.encoder(spectra, src_key_padding_mask=spectra_mask)

        return spectra, spectra_mask

    @torch.no_grad()
    def pred(
            self,
            spectra: Tensor,
            spectra_mask: Tensor,
            precursors: Tensor,
            tokens: Tensor,
    ) -> Tensor:
        """Inference mode: only return DDA scores (for inference/prediction).

        Args:
            spectra: (B, N, 2) - [mz, intensity]
            spectra_mask: (B, N) - bool, True=padding
            precursors: (B, 4) - [mass, charge, delta_rt, pred_rt]
            tokens: (B, L) - token ids (original tokens; no masking needed)

        Returns:
            dda_scores: (B,) - PSM scores in [0, 1] after sigmoid
        """
        # In inference, no masking is required
        token_mask = torch.zeros_like(tokens, dtype=torch.bool)  # all False, no masked positions

        # ============ 1. Encode spectrum ============
        spectra_encoded, spectra_mask_encoded = self._encoder(spectra, spectra_mask)

        # ============ 2. Decode with peptide ============
        decoder_output, _, peptide_emb = self._decoder(
            spectra_encoded, spectra_mask_encoded, precursors, tokens, token_mask
        )

        # ============ 3. DDA prediction ============
        # Do not use dropout (automatically disabled in eval mode)
        dda_logits = self.relu(self.psm_1(peptide_emb))
        dda_logits = self.psm_2(dda_logits).squeeze(-1)  # (B,)

        # ============ 4. Apply sigmoid to convert logits to probabilities ============
        dda_scores = torch.sigmoid(dda_logits)  # (B,)

        return dda_scores

    def _decoder(
            self,
            memory: Tensor,                      # (B, 1+N, D)
            memory_key_padding_mask: Tensor,     # (B, 1+N)
            precursors: Tensor,                  # (B, 4)
            tokens: Tensor,                      # (B, L)
            token_mask: Tensor,                  # (B, L) - bool, True=masked positions
    ) -> tuple[Tensor, List[List[int]], Tensor]:
        """Peptide decoder

        Args:
            memory: (B, 1+N, D) - encoded spectrum
            memory_key_padding_mask: (B, 1+N)
            precursors: (B, 4)
            tokens: (B, L)
            token_mask: (B, L) - bool, True=masked positions

        Returns:
            decoder_output: (B, 1+L_compressed, D)
            mask_positions: List[List[int]] - indices of masked positions in tokens_emb
            peptide_emb: (B, D) - decoder_output[:, 0, :]
        """
        decoder_output, mask_positions, peptide_emb = self.spectrum_sequence_encoder(
            memory, memory_key_padding_mask, precursors, tokens, token_mask
        )

        return decoder_output, mask_positions, peptide_emb

    # Load a .pt checkpoint
    @classmethod
    def load_pt(
        cls,
        path: str,
        config: dict,
        dtype: torch.dtype | None = None,
        map_location: str = "cpu",
    ) -> nn.Module:
        """
        More robust weight loading:
        - Supports DeepSpeed mp_rank_00_model_states.pt / pytorch_model.bin, etc.
        - Supports legacy formats like {'module': state_dict} or a raw state_dict
        - Automatically handles possible 'module.' prefixes
        - Prints detailed missing/unexpected keys to help debugging
        """
        import glob
        import os

        def _find_actual_file(p: str) -> str:
            # If p is a directory, try to locate common checkpoint filenames inside
            if os.path.isdir(p):
                # 1) Prefer DeepSpeed mp_rank_00
                candidates = sorted(
                    glob.glob(os.path.join(p, "mp_rank_00*_model_states.pt"))
                )
                if candidates:
                    return candidates[0]

                # 2) Common consolidated model filenames
                for name in ["pytorch_model.bin", "model.pt", "model.bin", "checkpoint.pt"]:
                    f = os.path.join(p, name)
                    if os.path.isfile(f):
                        return f

                raise FileNotFoundError(
                    f"Cannot find any model file in directory: {p}. "
                    f"Tried mp_rank_00*_model_states.pt, pytorch_model.bin, model.pt, model.bin, checkpoint.pt."
                )

            # If p is a file, return it directly
            if os.path.isfile(p):
                return p

            raise FileNotFoundError(f"Checkpoint path does not exist: {p}")

        ckpt_file = _find_actual_file(path)
        print(f"[load_pt] loading checkpoint file: {ckpt_file}")
        ckpt = torch.load(ckpt_file, map_location=map_location)

        # Extract state_dict
        if isinstance(ckpt, dict):
            if "module" in ckpt and isinstance(ckpt["module"], dict):
                raw_state_dict = ckpt["module"]
                print("[load_pt] Use ckpt['module'] as state_dict.")
            elif "state_dict" in ckpt and isinstance(ckpt["state_dict"], dict):
                raw_state_dict = ckpt["state_dict"]
                print("[load_pt] Use ckpt['state_dict'] as state_dict.")
            else:
                # Default: treat ckpt itself as a state_dict
                raw_state_dict = ckpt
                print("[load_pt] Use ckpt itself as state_dict.")
        else:
            raise ValueError(f"Unexpected checkpoint format in file: {ckpt_file}")

        # Check if there is a common prefix (e.g., 'module.') that should be stripped
        example_key = next(iter(raw_state_dict.keys()))
        print(f"[load_pt] example key in ckpt: {repr(example_key)}")

        def _maybe_strip_prefix(state_dict: dict, prefix: str) -> dict:
            if not example_key.startswith(prefix):
                return state_dict
            print(f"[load_pt] Detected prefix '{prefix}' in keys, stripping it.")
            return {k[len(prefix):]: v for k, v in state_dict.items()}

        # Enable prefix stripping as needed; if keys do not start with 'module.', it is skipped automatically
        state_dict = _maybe_strip_prefix(raw_state_dict, "module.")
        # You can extend here for other prefixes if needed:
        # state_dict = _maybe_strip_prefix(state_dict, "model.")

        # Build the model architecture (must match training)
        model = cls(
            dim_model=config["dim_model"],
            n_head=config["n_head"],
            dim_feedforward=config["dim_feedforward"],
            n_layers=config["n_layers"],
            dropout=config["dropout"],
            max_length=config["max_length"],
            vocab=config["vocab"],
            max_charge=config["max_charge"],
        )

        model_sd = model.state_dict()

        model_keys = set(model_sd.keys())
        ckpt_keys = set(state_dict.keys())

        # ---- More detailed debug prints ----
        print(f"[load_pt] #model_keys={len(model_keys)}, #ckpt_keys={len(ckpt_keys)}")

        miss_keys = sorted(model_keys - ckpt_keys)
        extra_keys = sorted(ckpt_keys - model_keys)

        if miss_keys:
            print("[load_pt] Missing keys in checkpoint (expected by model but not found in ckpt):")
            for i, k in enumerate(miss_keys[:200]):  # print at most 200 entries
                print(f"  [{i}] {repr(k)}")
            if len(miss_keys) > 200:
                print(f"  ... and {len(miss_keys) - 200} more")

        if extra_keys:
            print("[load_pt] Unexpected keys in checkpoint (present in ckpt but not used by model):")
            for i, k in enumerate(extra_keys[:200]):
                print(f"  [{i}] {repr(k)}")
            if len(extra_keys) > 200:
                print(f"  ... and {len(extra_keys) - 200} more")

        # Load weights
        load_info = model.load_state_dict(state_dict, strict=False)

        # PyTorch missing/unexpected keys are authoritative; print them as well
        if load_info.missing_keys:
            print("[load_pt] torch.load_state_dict missing_keys:")
            for i, k in enumerate(load_info.missing_keys):
                print(f"  [{i}] {repr(k)}")

        if load_info.unexpected_keys:
            print("[load_pt] torch.load_state_dict unexpected_keys:")
            for i, k in enumerate(load_info.unexpected_keys):
                print(f"  [{i}] {repr(k)}")

        # dtype cast (e.g., bf16 for evaluation)
        if dtype is not None:
            model = model.to(dtype)

        return model


def compute_mlm_loss_with_positions(
    logits: torch.Tensor,             # (B, vocab_size, L_compressed)
    mask_positions: torch.Tensor,     # (B, L_original) - int64, -1 indicates invalid positions
    tokens_label: torch.Tensor,       # (B, L_original)
    token_mask: torch.Tensor,         # (B, L_original) bool
    label_smoothing: float = 0.0,
):
    """
    Fully vectorized MLM loss computation (zero Python for-loops, pure GPU, BF16-friendly).

    Args:
        logits: (B, vocab_size, L_compressed)
        mask_positions: (B, L_original) - compressed indices; -1 indicates non-masked positions
        tokens_label: (B, L_original) - target tokens
        token_mask: (B, L_original) - bool mask
        label_smoothing: float

    Returns:
        loss: scalar tensor
    """
    device = logits.device
    dtype = logits.dtype
    B, vocab_size, L_compressed = logits.shape
    _, L_original = tokens_label.shape

    # ========== Step 1: transpose logits ==========
    logits = logits.transpose(1, 2)  # (B, L_compressed, vocab_size)

    # ========== Step 2: build validity mask (fully vectorized) ==========
    # Valid if: token_mask=True AND mask_positions>=0 AND tokens_label>0
    valid_mask = (
        token_mask
        & (mask_positions >= 0)
        & (mask_positions < L_compressed)
        & (tokens_label > 0)
    )  # (B, L_original)

    num_valid = valid_mask.sum()
    if num_valid == 0:
        return torch.tensor(0.0, device=device, dtype=dtype)

    # ========== Step 3: batch gather logits (key optimization) ==========
    # Clamp mask_positions to a safe range to avoid out-of-bounds access
    safe_positions = mask_positions.clamp(0, L_compressed - 1)  # (B, L_original)

    # Expand to match vocab_size
    safe_positions_expanded = safe_positions.unsqueeze(-1).expand(
        B, L_original, vocab_size
    )  # (B, L_original, vocab_size)

    # Batch gather: extract all requested positions at once
    gathered_logits = torch.gather(
        logits,                 # (B, L_compressed, vocab_size)
        1,                      # dim=1
        safe_positions_expanded # (B, L_original, vocab_size)
    )  # (B, L_original, vocab_size)

    # ========== Step 4: flatten and filter (vectorized) ==========
    gathered_logits_flat = gathered_logits.reshape(-1, vocab_size)  # (B*L_original, vocab_size)
    tokens_label_flat = tokens_label.reshape(-1)                    # (B*L_original,)
    valid_mask_flat = valid_mask.reshape(-1)                        # (B*L_original,)

    # Filter via boolean indexing (pure GPU operation)
    pred_logits_valid = gathered_logits_flat[valid_mask_flat]  # (num_valid, vocab_size)
    true_tokens_valid = tokens_label_flat[valid_mask_flat]     # (num_valid,)

    # ========== Step 5: compute loss (fused op) ==========
    loss = F.cross_entropy(
        pred_logits_valid,
        true_tokens_valid,
        label_smoothing=label_smoothing,
        reduction="mean",
    )

    return loss
