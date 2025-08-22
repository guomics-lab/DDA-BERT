from __future__ import annotations

import os
import polars as pl
import pandas as pd
import numpy as np
import yaml
from tqdm import tqdm
from sklearn.metrics import roc_auc_score, accuracy_score
import logging
import argparse
import random
import glob
import datetime

import torch
from torch import nn
from torch import Tensor
from torch.utils.data import DataLoader
import lightning.pytorch as ptl
from lightning.pytorch.strategies import DDPStrategy
from torch.utils.tensorboard import SummaryWriter

from transformer.dataset import SpectrumDataset, collate_batch_weight_deltaRT, mkdir_p
from transformer.model import DDABert
from transformer.iterable_dataset import create_iterable_dataset

from utils import set_seeds

logger = logging.getLogger()
logger.setLevel(logging.INFO)

import warnings
warnings.filterwarnings("ignore")

class Evalute(ptl.LightningModule):
    """evaluate for model."""

    def __init__(
            self,
            config: dict[str, Any],
            model: DDABert,
            sw: SummaryWriter,
            epoch: int,
    ) -> None:
        super().__init__()
        self.config = config
        self.model = model
        self.sw = sw
        self.epoch = epoch
        
        self._reset_metrics()
        self.loss_fn = nn.BCEWithLogitsLoss()
        self.sigmod = nn.Sigmoid()

        self.val_running_loss = None
        self.val_weight_running_loss = None
        self.val_steps = 0
        self.val_step_scale = config["step_scale"]

    def test_step(
            self,
            batch: tuple[Tensor, Tensor, Tensor, Tensor, list, list],
    ) -> torch.Tensor:
        """Single test step."""
        spectra, spectra_mask, precursors, tokens, peptide, label, weight = batch
        
        spectra = spectra.to(self.device).to(torch.bfloat16)
        spectra_mask = spectra_mask.to(self.device).to(torch.bfloat16)
        precursors = precursors.to(self.device).to(torch.bfloat16)
        tokens = tokens.to(self.device).to(torch.long)
        label = label.to(self.device).to(torch.bfloat16)
        weight = weight.to(self.device).to(torch.bfloat16)
        
        # Set the weight of false targets to 0.3
        weight[weight < 0.6] = 0.3

        loss_fn_weight = nn.BCEWithLogitsLoss(weight=weight)

        # Loss
        with torch.autocast(device_type='cuda', dtype=torch.bfloat16):
            pred, _ = self.model(spectra, spectra_mask, precursors, tokens)

            # calc loss with logit
            loss_weight = loss_fn_weight(pred, label.flatten())
            loss = self.loss_fn(pred, label.flatten())

        pred = self.sigmod(pred)
        
        if self.val_running_loss is None:
            self.val_running_loss = loss.item()
            self.val_weight_running_loss = loss_weight.item()
        else:
            self.val_running_loss = 0.99 * self.val_running_loss + (1 - 0.99) * loss.item()
            self.val_weight_running_loss = 0.99 * self.val_weight_running_loss + (1 - 0.99) * loss_weight.item()

        if ((self.val_steps + 1) % int(self.val_step_scale)) == 0:
            self.sw.add_scalar("eval_step/dda_loss_raw", loss.item(), self.val_steps - 1)
            self.sw.add_scalar("eval_step/dda_loss_smooth", self.val_running_loss, self.val_steps - 1)
            
            self.sw.add_scalar("eval_step/dda_loss_weight_raw", loss_weight.item(), self.val_steps - 1)
            self.sw.add_scalar("eval_step/dda_loss_weight_smooth", self.val_weight_running_loss, self.val_steps - 1)
        self.val_steps += 1

        label = label.to(torch.float32).cpu().detach().numpy()
        pred = pred.to(torch.float32).cpu().detach().numpy()
        weight = weight.to(torch.float32).cpu().detach().numpy()
        
        self.pred_list.extend(pred.tolist())
        self.label_list.extend(label.tolist())        
        self.weight_list.extend(weight.tolist())

        self.loss_list.extend([loss.item()])
        self.loss_weight_list.extend([loss_weight.item()])

    def on_test_end(self) -> None:
        df = pd.DataFrame({"score": self.pred_list,
                           "label": self.label_list,
                           "weight": self.weight_list})
        df.to_csv(os.path.join(self.config['out_path'], f"epoch{self.epoch}_pred.csv"), mode='a+', header=False, index=None)
        
        df['type'] = np.where(df['label'] == 1.0, 'target', \
                     np.where(df['weight'] < 0.9, 'false_target',  'decoy'))
        
        df = df[df['type'] != 'false_target']
        label_list = torch.tensor(df['label'].to_numpy())
        pred_list = torch.tensor(df['score'].to_numpy())
        auc = roc_auc_score(label_list, pred_list)
        
        pred = torch.tensor(np.array([1 if v >= 0.5 else 0 for v in pred_list]))
        acc = accuracy_score(label_list, pred)
        
        # calc loss
        loss_fn = nn.BCELoss()
        loss = loss_fn(torch.tensor(df['score'].to_numpy()),\
                       torch.tensor(df['label'].to_numpy())).item()
        
        loss_fn_weight = nn.BCEWithLogitsLoss(weight=torch.tensor(df['weight'].to_numpy()))
        loss_weight = loss_fn_weight(torch.tensor(df['score'].to_numpy()),\
                                     torch.tensor(df['label'].to_numpy())).item()
        logging.info(f"epoch: {self.epoch} - auc: {auc:.3f} - acc: {acc:.3f} - loss: {loss:.3f} - loss_weight: {loss_weight:.3f}")
        
        # calc loss by source
        df['label'] = df['label'].astype(float)
        df['score'] = df['score'].astype(float)
        df['weight'] = df['weight'].astype(float)

        loss_fn = nn.BCELoss()
        target = df[df['type'] == 'target']
        target_loss = loss_fn(torch.tensor(target['score'].to_numpy()),\
                              torch.tensor(target['label'].to_numpy())).item()

        decoy = df[df['type'] == 'decoy']
        decoy_loss = loss_fn(torch.tensor(decoy['score'].to_numpy()),\
                             torch.tensor(decoy['label'].to_numpy())).item()

        self.sw.add_scalar("eval_epoch/auc", auc, int(self.epoch))
        self.sw.add_scalar("eval_epoch/acc", acc, int(self.epoch))
        self.sw.add_scalar("eval_epoch/dda_loss", loss, int(self.epoch))
        self.sw.add_scalar("eval_epoch/dda_loss_weight", loss_weight, int(self.epoch))
        
        # add monitor
        self.sw.add_scalar("eval_epoch/target_loss", target_loss, int(self.epoch))
        self.sw.add_scalar("eval_epoch/decoy_loss", decoy_loss, int(self.epoch))

    def _reset_metrics(self) -> None:
        self.pred_list = []
        self.label_list = []
        self.weight_list = []
        self.loss_list = []
        self.loss_weight_list = []

    
def combine_data(file_list):
    df_ret = []
    use_cols = ['mz_array', 'intensity_array', 'precursor_mz',
                'precursor_charge', 'modified_sequence', 'label', 'weight', 'delta_rt_model', 'predicted_rt']
    
    for file in file_list:
        df = pl.read_ipc(file, columns=use_cols)
        df = df.with_columns(pl.col("precursor_charge").cast(pl.Int8))
        df = df.with_columns(pl.col("label").cast(pl.Int8))

        for col in ['delta_rt_model', 'predicted_rt', 'precursor_mz', 'weight']:
            df = df.with_columns(pl.col(col).cast(pl.Float32))
            
        df_ret.append(df)

    df = pl.concat(df_ret)
    df = df.select(pl.col(use_cols))
    return df

def gen_dl(df, config):
    s2i = {v: i for i, v in enumerate(config["vocab"])}
    logging.info(f"gen Vocab: {s2i}")
    
    if "weight" not in df.columns:
        df = df.with_columns(pl.lit(1.0).alias("weight"))
        
    if df['predicted_rt'].max() > 1.2:
        # predicted_rt normalization
        df = df.with_columns((pl.col("predicted_rt") / pl.col("predicted_rt").max()).alias("predicted_rt"))
     
    ds = SpectrumDataset(df, s2i, config["n_peaks"], need_label=True, need_weight=True, need_deltaRT=True)
    dl = DataLoader(ds,
                    batch_size=config["predict_batch_size"],
                    shuffle=False,
                    num_workers=0,
                    collate_fn=collate_batch_weight_deltaRT,
                    timeout=10)
    logging.info(f"Data: {len(ds):,} samples, DataLoader: {len(dl):,}")
    return dl


def main() -> None:
    """Predict with the model."""
    logging.info("Initializing inference.")

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="/DDA_BERT/yaml/eval_model.yaml")
    args = parser.parse_args()
    
    # load config
    with open(args.config) as f_in:
        config = yaml.safe_load(f_in)
    vocab = ['<pad>', '<mask>'] + list(config["residues"].keys()) + ['<unk>']
    config["vocab"] = vocab
    s2i = {v: i for i, v in enumerate(vocab)}
    logging.info(f"Vocab: {s2i}, n_peaks: {config['n_peaks']}")
    
    mkdir_p(config["out_path"])
    logging.info(f"config:  {args.config}")

    set_seeds(config['seed'])
    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    config["tb_summarywriter"] = os.path.join(config["tb_summarywriter"], datetime.datetime.now().strftime("EVAL_%y_%m_%d_%H_%M"))
    sw = SummaryWriter(config["tb_summarywriter"])

    logging.info(f"valid_path:  {config['valid_path']}")
    if ';' in config['valid_path']:
        valid_path_list = config['valid_path'].split(';')
        
        ipc_file_path = []
        for valid_path in valid_path_list:
            ipc_file_path.extend(glob.glob(f'{valid_path}/*ipc'))
        df = combine_data(ipc_file_path)
        dl = gen_dl(df, config)
        
    elif os.path.isdir(config['valid_path']):
        ipc_file_path = [os.path.join(config['valid_path'], f) for f in os.listdir(config['valid_path']) if f.endswith('.ipc')][:1]
        print('load file: ', ipc_file_path)
        df = combine_data(ipc_file_path)
        dl = gen_dl(df, config)
        
    else:
        df = pl.read_ipc(config['valid_path'])
        dl = gen_dl(df, config)

    # evaluate
    gpu_num = torch.cuda.device_count()
    
    # Update parameters
    one_epoch_iters = int(len(dl))
    config["step_scale"] = int(one_epoch_iters * float(config["step_ratio"]))

    strategy = DDPStrategy(gradient_as_bucket_view=True, find_unused_parameters=True)
    trainer = ptl.Trainer(
        accelerator="auto",
        devices="auto",
        strategy=strategy,
    )
    
    logging.info(f"model_path:  {config['model_path']}")
    if os.path.isdir(config['model_path']):
        model_path_list = glob.glob(f"{config['model_path']}/*/*.pt")

        if config['model_type'] == 'descend':
            model_step_id = sorted(range(len(model_path_list)), key=lambda x: int(model_path_list[x].split('/')[-2]), reverse=True)
        elif config['model_type'] == 'ascend':
            model_step_id = sorted(range(len(model_path_list)), key=lambda x: int(model_path_list[x].split('/')[-2]))

        model_path_list = [model_path_list[i] for i in model_step_id]

    else:
        # load single model
        model_path_list = [config['model_path']]
    
    print('model_path_list: ', model_path_list)
    for model_path in model_path_list:
        # load model
        model = DDABert.load_pt(model_path, config)
        model.eval()
        model.to(torch.bfloat16).to(device)

        # parse epoch
        epoch = model_path.split('/')[-2]
        evaluate = Evalute(config, model, sw, epoch)
        trainer.test(evaluate, dataloaders=dl)


if __name__ == "__main__":
    main()