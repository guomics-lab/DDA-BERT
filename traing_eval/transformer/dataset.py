from __future__ import annotations

import re
import os

import numpy as np
import pandas as pd
import polars as pl
import spectrum_utils.spectrum as sus
import torch
from torch import nn
from torch import Tensor
from torch.utils.data import Dataset
import torch.nn.functional as F

PROTON_MASS_AMU = 1.007276


class SpectrumDataset(Dataset):
    """Spectrum dataset class supporting `.ipc`."""

    def __init__(
        self,
        df: pd.DataFrame | pl.DataFrame,
        s2i: dict,
        n_peaks: int = 200,
        max_length=50,
        min_mz: float = 50.0,
        max_mz: float = 2500.0,
        min_intensity: float = 0.01,
        remove_precursor_tol: float = 2.0,
        reverse_peptide: bool = True,
        annotated: bool = True,
    ) -> None:
        super().__init__()
        self.df = df
        self.s2i = s2i
        self.n_peaks = n_peaks
        self.max_length = max_length
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.remove_precursor_tol = remove_precursor_tol
        self.min_intensity = min_intensity
        self.annotated = annotated

        assert isinstance(df, pl.DataFrame):

    def __len__(self) -> int:
        return int(self.df.shape[0])

    def __getitem__(self, idx: int) -> tuple[Tensor, float, int, Tensor | list[str]]:
        mz_array = torch.Tensor(self.df[idx, "mz_array"].to_list())
        int_array = torch.Tensor(self.df[idx, "intensity_array"].to_list())
        precursor_mz = self.df[idx, "precursor_mz"]
        precursor_charge = self.df[idx, "precursor_charge"]
        peptide = self.df[idx, "modified_sequence"]
        label = self.df[idx, "label"]
        weight = self.df[idx, "weight"]
        deltaRT = self.df[idx, "delta_rt_model"] or 0.0
        predictedRT = self.df[idx, "predicted_rt"] or 0.0

        spectrum = self._process_peaks(mz_array, int_array, precursor_mz, precursor_charge)
        tokens = self._tokenize(peptide)
        return spectrum, precursor_mz, precursor_charge, deltaRT, predictedRT, tokens, peptide, label, weight

    def _process_peaks(
        self,
        mz_array: Tensor,
        int_array: Tensor,
        precursor_mz: Tensor,
        precursor_charge: Tensor,
    ) -> Tensor:
        """Preprocess the spectrum by removing noise peaks and scaling the peak intensities.

        Parameters
        ----------
        mz_array : numpy.ndarray of shape (n_peaks,)
            The spectrum peak m/z values.
        int_array : numpy.ndarray of shape (n_peaks,)
            The spectrum peak intensity values.

        Returns
        -------
        torch.Tensor of shape (n_peaks, 2)
            A tensor of the spectrum with the m/z and intensity peak values.
        """
        
        spectrum = sus.MsmsSpectrum(
            "",
            precursor_mz,
            precursor_charge,
            np.array(mz_array).astype(np.float32),
            np.array(int_array).astype(np.float32),
        )
        try:
            spectrum.set_mz_range(self.min_mz, self.max_mz)
            if len(spectrum.mz) == 0:
                raise ValueError
            spectrum.remove_precursor_peak(self.remove_precursor_tol, "Da")
            if len(spectrum.mz) == 0:
                raise ValueError
            spectrum.filter_intensity(self.min_intensity, self.n_peaks)
            if len(spectrum.mz) == 0:
                raise ValueError
            spectrum.scale_intensity("root", 1)
            intensities = spectrum.intensity / np.linalg.norm(spectrum.intensity)
            return torch.tensor(np.array([spectrum.mz, intensities])).T.float()
        except ValueError:
            # Replace invalid spectra by a dummy spectrum.
            return torch.tensor([[0, 1]]).float()


    def _tokenize(self, sequence):
        """Transform a peptide sequence into tokens

        Parameters
        ----------
        sequence : str
            A peptide sequence.

        Returns
        -------
        torch.Tensor
            The token for each amino acid in the peptide sequence.
        """
        sequence = sequence.replace("I", "L").replace('n[42]', 'X')
        sequence = sequence.replace('cC', 'C[57.02]')\
                           .replace('oxM', 'M[15.99]')\
                           .replace('M(ox)', 'M[15.99]')\
                           .replace('deamN', 'N[.98]')\
                           .replace('deamQ', 'Q[.98]')\
                           .replace('a', 'X')
        
        sequence = re.split(r"(?<=.)(?=[A-Z])", sequence)

        tokens = torch.tensor([self.s2i[aa] for aa in sequence])
        tokens = F.pad(tokens, (0, self.max_length - tokens.shape[0]), 'constant', 0)# padding
        return tokens


def padding(data):
    ll = torch.tensor([x.shape[0] for x in data], dtype=torch.long)
    data = nn.utils.rnn.pad_sequence(data, batch_first=True)
    data_mask = torch.arange(data.shape[1], dtype=torch.long)[None, :] >= ll[:, None]
    return data, data_mask


def collate_batch_weight_deltaRT(
    batch: list[tuple[Tensor, float, int, Tensor, Tensor]]
) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor, Tensor]:
    """Collate batch of samples."""
    spectrum, precursor_mzs, precursor_charges, deltaRT, predictedRT, tokens, peptide, label, weight = zip(*batch)
    
    # Pad spectra
    spectra, spectra_mask = padding(spectrum)
    
    # stack tokens
    tokens = torch.stack(tokens, dim=0)

    precursor_mzs = torch.tensor(precursor_mzs)
    precursor_charges = torch.tensor(precursor_charges)
    precursor_masses = (precursor_mzs - PROTON_MASS_AMU) * precursor_charges

    deltaRT = torch.tensor(deltaRT)
    predictedRT = torch.tensor(predictedRT)
    precursors = torch.vstack([precursor_masses, precursor_charges, deltaRT, predictedRT]).T.float()
    
    label = torch.tensor(label).to(torch.float)
    weight = torch.tensor(weight).to(torch.float)
    return spectra, spectra_mask, precursors, tokens, peptide, label, weight
