from __future__ import annotations

import torch
from torch import Tensor
from torch import nn

from ..transformer.decoder import PeptideDecoder
from ..transformer.layers import MultiScalePeakEmbedding


def transform(m: nn.Module) -> nn.Module:
    gm: torch.fx.GraphModule = torch.fx.symbolic_trace(m)

    # Recompile the forward() method of `gm` from its Graph
    gm.recompile()

    return gm


class MaskedLanguageModel(nn.Module):
    """
    predicting origin token from masked input sequence
    n-class classification problem, n-class = vocab_size
    """

    def __init__(self, hidden, vocab_size):
        """
        :param hidden: output size of BERT model
        :param vocab_size: total vocab size
        """
        super().__init__()
        self.linear = nn.Linear(hidden, vocab_size)

    def forward(self, x):
        return self.linear(x)


class DDA_BERT(nn.Module):

    def __init__(
            self,
            dim_model: int = 768,
            n_head: int = 16,
            dim_feedforward: int = 2048,
            n_layers: int = 9,
            dropout: float = 0,
            max_length: int = 50,  # peptide max_length
            vocab_size: int = 29,  # vocab size of amino acid
            max_charge: int = 10,
    ) -> None:
        super().__init__()
        self.dim_model = dim_model
        self.max_length = max_length

        # The latent representations for the spectrum and each of its peaks.
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

        # spectrum sequence encoder, where no spectrum mask
        self.spectrum_sequence_encoder = PeptideDecoder(
            dim_model=dim_model,
            n_head=n_head,
            dim_feedforward=dim_feedforward,
            n_layers=n_layers,
            dropout=dropout,
            residues_length=vocab_size,
            max_charge=max_charge,
            hidden_size=max_length  # peptide max_length
        )

        self.psm_0 = nn.Linear(self.max_length + 1, 1)
        self.psm_1 = nn.Linear(self.dim_model, 64)
        self.psm_2 = nn.Linear(64, 1)

        self.mask_lm = MaskedLanguageModel(self.dim_model, vocab_size)

        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(p=dropout)

    @classmethod
    def load_ckpt(cls, path: str) -> nn.Module:
        """Load model from checkpoint."""
        ckpt = torch.load(path, map_location="cpu")
        config = ckpt["config"]

        # fix
        if 'loss_fn.weight' in ckpt["state_dict"]:
            ckpt["state_dict"].pop('loss_fn.weight')

        # check if PTL checkpoint
        if all([x.startswith("model") for x in ckpt["state_dict"].keys()]):
            ckpt["state_dict"] = {k.replace("model.", ""): v for k, v in ckpt["state_dict"].items() if
                                  k.startswith('model.')}

        model = cls(
            dim_model=config["dim_model"],
            n_head=config["n_head"],
            dim_feedforward=config["dim_feedforward"],
            n_layers=config["n_layers"],
            dropout=config["dropout"],
            max_length=config["max_length"],
            vocab_size=len(config["vocab"]),
            max_charge=config["max_charge"],
        )
        model.load_state_dict(ckpt["state_dict"], strict=False)
        return model, config

    @classmethod
    def load_ckpt_resume(cls, path: str) -> nn.Module:
        """Load model from checkpoint."""
        ckpt = torch.load(path, map_location="cpu")
        config = ckpt["config"]

        if 'loss_fn.weight' in ckpt["model"].keys():
            ckpt["model"].pop('loss_fn.weight')

        model = cls(
            dim_model=config["dim_model"],
            n_head=config["n_head"],
            dim_feedforward=config["dim_feedforward"],
            n_layers=config["n_layers"],
            dropout=config["dropout"],
            max_length=config["max_length"],
            vocab_size=len(config["vocab"]),
            max_charge=config["max_charge"],
        )
        model.load_state_dict(ckpt["model"], strict=False)
        return model, ckpt

    @classmethod
    def load_pt(cls, path: str, config: str) -> nn.Module:
        """Load model from checkpoint."""
        ckpt = torch.load(path, map_location="cpu")
        model = cls(
            dim_model=config["dim_model"],
            n_head=config["n_head"],
            dim_feedforward=config["dim_feedforward"],
            n_layers=config["n_layers"],
            dropout=config["dropout"],
            max_length=config["max_length"],
            vocab_size=len(config["vocab"]),
            max_charge=config["max_charge"],
        )

        miss_keys = set(model.state_dict().keys()) - set(ckpt['module'].keys())
        if len(miss_keys) > 0:
            print('miss keys： ', miss_keys)

        model.load_state_dict(ckpt['module'], strict=False)
        return model

    def forward(
            self,
            spectra: Tensor,
            spectra_mask: Tensor,
            precursors: Tensor,
            tokens: Tensor,
    ) -> Tensor:
        """Model forward pass.

        Args:
            spectra: float Tensor (batch, n_peaks, 2) . 2: [mz_array, int_array]
            spectra_mask: Spectra padding mask, True for padded indices, bool Tensor (batch, n_peaks)
            precursors: float Tensor (batch, 4) . 4: [precursor_masses, precursor_charges, deltaRT, predictedRT]
            tokens: float Tensor (batch, 50)
        Returns:
            # PSM
            pred: float Tensor (batch, 1)
        """
        spectra, spectra_mask = self._encoder(spectra, spectra_mask)
        return self._psm_encoder(spectra, spectra_mask, precursors, tokens)

    def pred(
            self,
            spectra: Tensor,
            spectra_mask: Tensor,
            precursors: Tensor,
            tokens: Tensor,
    ) -> Tensor:
        """Model forward pass.

        Args:
            spectra: float Tensor (batch, n_peaks, 2) . 2: [mz_array, int_array]
            spectra_mask: Spectra padding mask, True for padded indices, bool Tensor (batch, n_peaks)
            precursors: float Tensor (batch, 4) . 4: [precursor_masses, precursor_charges, deltaRT, predictedRT]
            tokens: float Tensor (batch, 50)
        Returns:
            # PSM
            dda_pred: float Tensor (batch)
            # MaskedLanguageModel 
            mask_pred: float Tensor (batch)
        """
        with torch.no_grad():
            dda_pred, mask_pred = self.forward(spectra, spectra_mask, precursors, tokens)
            sigmod = nn.Sigmoid()
            dda_pred = sigmod(dda_pred)
        return dda_pred, mask_pred

    # SpectrumEncoder
    def _encoder(self, spectra: Tensor, spectra_mask: Tensor) -> tuple[Tensor, Tensor]:
        # Peak encoding
        spectra = self.peak_encoder(spectra[:, :, [0]], spectra[:, :, [1]])

        # Self-attention on latent spectra AND peaks
        latent_spectra = self.latent_spectrum.expand(spectra.shape[0], -1, -1)
        spectra = torch.cat([latent_spectra, spectra], dim=1)
        latent_mask = torch.zeros((spectra_mask.shape[0], 1), dtype=bool, device=spectra_mask.device)
        spectra_mask = torch.cat([latent_mask, spectra_mask], dim=1).bool()

        spectra = self.encoder(spectra, src_key_padding_mask=spectra_mask)
        return spectra, spectra_mask

    def _psm_encoder(
            self,
            spectra: Tensor,
            spectra_mask: Tensor,
            precursors: Tensor,
            tokens: Tensor,
    ) -> Tensor:
        # pred： (batch , 51, 768)
        decoder_output = self.spectrum_sequence_encoder(spectra, spectra_mask, precursors, tokens)

        # pred： (batch, 51, 768) ==> (batch, 768)
        pred = decoder_output[:, 0, :]

        # mask pred：(batch, 51, 768) ==> (batch, 51, 29) ==> (batch, 29, 51) ==> (batch, 29, 50)
        mask_pred = self.mask_lm(decoder_output).transpose(1, 2)[:, :, 1:]

        # pred：(batch, 768) ==> (batch, 64)
        dda_pred = self.dropout(self.relu(self.psm_1(pred)))

        # preds：(batch, 64) ==> (batch, 1)
        dda_pred = self.psm_2(dda_pred).squeeze()

        return dda_pred, mask_pred
