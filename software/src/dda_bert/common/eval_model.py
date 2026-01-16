import os

import lightning.pytorch as ptl
import pandas as pd
import torch
from torch import Tensor

from ..transformer.model import DDA_BERT


class Evalute(ptl.LightningModule):
    """evaluate for model."""

    def __init__(
            self,
            out_path: str,
            file_name: str,
            model: DDA_BERT,
    ) -> None:
        super().__init__()
        self.out_path = out_path
        self.file_name = file_name
        self.model = model
        self._reset_metrics()

    def test_step(
            self,
            batch: tuple[Tensor, Tensor, Tensor, Tensor, list, list],
    ) -> torch.Tensor:
        """Single test step."""
        spectra, spectra_mask, precursors, tokens, peptide, label, weight, index = batch

        spectra = spectra.to(self.device).to(torch.bfloat16)
        spectra_mask = spectra_mask.to(self.device).to(torch.bfloat16)
        precursors = precursors.to(self.device).to(torch.bfloat16)
        tokens = tokens.to(self.device).to(torch.long)

        # Loss
        with torch.autocast(device_type='cuda', dtype=torch.bfloat16):
            pred, _ = self.model.pred(spectra, spectra_mask, precursors, tokens)

        label = label.to(torch.int32).cpu().detach().numpy()
        pred = pred.to(torch.float32).cpu().detach().numpy()
        index = index.to(torch.int32).cpu().detach().numpy()
        weight = weight.to(torch.float32).cpu().detach().numpy()

        if label.shape[0] > 1 and (not isinstance(pred, float)) and pred is not None:
            self.pred_list.extend(pred.tolist())
            self.label_list.extend(label.tolist())
            self.index_list.extend(index.tolist())
            self.weight_list.extend(weight.tolist())

    def on_test_end(self) -> None:
        df = pd.DataFrame({"index": self.index_list,
                           "score": self.pred_list,
                           "label": self.label_list,
                           "weight": self.weight_list})
        df.to_csv(os.path.join(self.out_path, f"%s_batch_pred.csv" % self.file_name), mode='a+', header=False, index=None)

    def _reset_metrics(self) -> None:
        self.pred_list = []
        self.label_list = []
        self.index_list = []
        self.weight_list = []
