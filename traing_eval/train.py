import torch
import argparse
import datetime
import logging
import os
os.environ['CUDA_LAUNCH_BLOCKING'] = '1'

import deepspeed
from deepspeed.accelerator import get_accelerator
import matplotlib.pyplot as plt
import numpy as np
import json
import yaml
import random
import pandas as pd
import polars as pl
import math
from tqdm import tqdm

from torch import nn
from torch.utils.tensorboard import SummaryWriter
from torch.optim import Optimizer

from deepspeed.runtime.checkpoint_engine.torch_checkpoint_engine import TorchCheckpointEngine

from transformer.utils import mkdir_p, set_seeds
from transformer.iterable_dataset_online import create_iterable_dataset
from transformer.model import DDABert, compute_mlm_loss_with_positions
from transformer.data_augment import mask_batch_with_augmentation

import warnings
warnings.filterwarnings("ignore")

logger = logging.getLogger()
logger.setLevel(logging.INFO)


def add_argument():
    parser = argparse.ArgumentParser(description='DDABert')

    # cuda / distributed
    parser.add_argument('--with_cuda',
                        default=True,
                        action='store_true',
                        help='use CPU in case there\'s no GPU support')
    parser.add_argument('--node_num',
                        default=1,
                        type=int,
                        help='the num of mpi node')
    parser.add_argument('--gpu_num',
                        default=2,
                        type=int,
                        help='the total num of gpu')
    parser.add_argument('--local_rank',
                        type=int,
                        default=-1,
                        help='local rank passed from distributed launcher')

    # train
    parser.add_argument("--config", default="/DDA_BERT/yaml/model.yaml")
    parser.add_argument('-f',
                        '--file_size',
                        default=128 * 1024,
                        type=int,
                        help='the size of file')

    # resume
    parser.add_argument('--resume',
                        action='store_true',
                        default=False,
                        help='resume training from checkpoint')
    parser.add_argument('--resume_tag',
                        type=str,
                        default=None,
                        help='tag (usually epoch) to resume from; if None, auto-detect latest')

    # deepspeed
    parser.add_argument(
        '--dtype',
        default='bf16',
        type=str,
        choices=['bf16', 'fp16', 'fp32'],
        help='Datatype used for training'
    )
    parser.add_argument(
        '--stage',
        default=0,
        type=int,
        choices=[0, 1, 2, 3],
        help='ZeRO stage'
    )

    # Include DeepSpeed configuration arguments
    parser = deepspeed.add_config_arguments(parser)
    args = parser.parse_args()
    return args


args = add_argument()
print("------------")
print(args.deepspeed)
print(args.deepspeed_config)


def get_rank():
    if torch.distributed.is_available() and torch.distributed.is_initialized():
        return torch.distributed.get_rank()
    return 0

def is_rank0():
    return get_rank() == 0

def get_torch_optimizer(optimizer):
    if isinstance(optimizer, Optimizer):
        return optimizer

    if hasattr(optimizer, 'optimizer') and isinstance(optimizer.optimizer, Optimizer):
        return optimizer.optimizer

    raise TypeError('{} is not a subclass of torch.optim.Optimizer'.format(type(optimizer).__name__))


class TwoPhaseScheduler(object):
    """
    Two-phase learning rate schedule:
      - Phase 1 (first half of training, step-based): warmup + decay (exp or cosine)
        * Linearly warm up from min_lr to max_lr over (warmup_ratio * phase1_max_iter) steps
        * Then decay from max_lr down to min_lr according to the chosen strategy
      - Phase 2 (second half of training, step-based): keep a constant second_min_lr

    Driven by the global step. Call step() every iteration (typically driven by DeepSpeed engine.step()).
    """

    def __init__(self,
                 optimizer: Optimizer,
                 phase1_max_iter: int,
                 warmup_ratio: float,
                 max_lr: float,
                 min_lr: float,
                 second_min_lr: float,
                 warmup_type: str = 'cos',
                 last_batch_iteration: int = -1):
        self.optimizer = get_torch_optimizer(optimizer)

        self.phase1_max_iter = max(1, int(phase1_max_iter))
        self.warmup_ratio = float(warmup_ratio)
        self.warmup_iters = max(1, int(self.phase1_max_iter * self.warmup_ratio))
        self.warmup_type = warmup_type

        self.max_lr = float(max_lr)
        self.min_lr = float(min_lr)
        self.second_min_lr = float(second_min_lr)

        self.last_batch_iteration = last_batch_iteration

        # Record the initial lrs (used for ratio-based scaling in some schedulers).
        # Here we use an absolute LR schedule to keep max_lr/min_lr semantics consistent.
        self.org_lrs = [group['lr'] for group in self.optimizer.param_groups]

        # Flag for the current phase
        self.in_phase1 = True

    def _phase1_lr(self, step: int) -> float:
        """Compute the absolute learning rate for phase 1 given the global step (starting from 0)."""
        if step <= self.warmup_iters:
            # Linear warmup: min_lr -> max_lr
            ratio = step / float(self.warmup_iters)
            lr = self.min_lr + ratio * (self.max_lr - self.min_lr)
            return lr

        # Decay part: max_lr -> min_lr
        t = step - self.warmup_iters
        T = max(1, self.phase1_max_iter - self.warmup_iters)

        if self.warmup_type == 'exp':
            # Exponential-like decay with exponent 0.9 (tunable)
            decay_ratio = (1 - t / float(T)) ** 0.9
            lr = self.min_lr + decay_ratio * (self.max_lr - self.min_lr)
        else:
            # Cosine decay
            lr = self.min_lr + 0.5 * (self.max_lr - self.min_lr) * (1 + np.cos(t / float(T) * np.pi))

        return float(lr)

    def get_lr(self):
        """Return the list of absolute lrs for each param_group."""
        if self.last_batch_iteration < 0:
            return [self.min_lr for _ in self.optimizer.param_groups]

        if self.last_batch_iteration <= self.phase1_max_iter:
            lr = self._phase1_lr(self.last_batch_iteration)
        else:
            lr = self.second_min_lr

        return [float(lr) for _ in self.optimizer.param_groups]

    def step(self, last_batch_iteration=None):
        if last_batch_iteration is None:
            last_batch_iteration = self.last_batch_iteration + 1
        self.last_batch_iteration = last_batch_iteration

        lrs = self.get_lr()
        for param_group, lr in zip(self.optimizer.param_groups, lrs):
            param_group['lr'] = lr
        self._last_lr = [group['lr'] for group in self.optimizer.param_groups]

    def get_last_lr(self):
        assert getattr(self, '_last_lr', None) is not None, "need to call step() first"
        return self._last_lr

    def state_dict(self):
        return {
            'last_batch_iteration': self.last_batch_iteration,
            'phase1_max_iter': self.phase1_max_iter,
            'warmup_ratio': self.warmup_ratio,
            'max_lr': self.max_lr,
            'min_lr': self.min_lr,
            'second_min_lr': self.second_min_lr,
            'warmup_type': self.warmup_type,
        }

    def load_state_dict(self, sd):
        self.last_batch_iteration = sd.get('last_batch_iteration', -1)
        self.phase1_max_iter = sd.get('phase1_max_iter', self.phase1_max_iter)
        self.warmup_ratio = sd.get('warmup_ratio', self.warmup_ratio)
        self.warmup_iters = max(1, int(self.phase1_max_iter * self.warmup_ratio))
        self.max_lr = sd.get('max_lr', self.max_lr)
        self.min_lr = sd.get('min_lr', self.min_lr)
        self.second_min_lr = sd.get('second_min_lr', self.second_min_lr)
        self.warmup_type = sd.get('warmup_type', self.warmup_type)


def save_checkpoint(model_engine, save_dir, tag, client_state=None):
    if client_state is None:
        client_state = {}
    tag = str(tag)
    model_engine.save_checkpoint(save_dir, tag=tag, client_state=client_state)


def calc_loss(mask, weight, label, pred):
    mask_weight = weight[mask]
    mask_label = label[mask]
    mask_pred = pred[mask]

    if len(mask_weight) > 0:
        dda_criterion = nn.BCEWithLogitsLoss(weight=mask_weight)
        dda_loss = dda_criterion(mask_pred, mask_label.flatten())
        return dda_loss.item()
    else:
        return 0


def train(model_engine, trainloader, sw, optim, scheduler, local_device,
          target_dtype, local_rank, config, epoch, sw_step):
    running_loss, dda_running_loss, mask_running_loss = None, None, None
    target_running_loss, ft_part1_running_loss, ft_part2_running_loss = None, None, None

    trainloader.set_epoch(epoch)
    train_bar = tqdm(
        trainloader, 
        total=len(trainloader),
        disable=True,
    )

    for batch_idx, batch in enumerate(train_bar):
        if batch_idx > config["one_epoch_iters"]:
            break
        (
            spectra,
            spectra_mask,
            precursors,
            tokens_input,
            tokens_label,
            token_mask,
            labels,
            weights,
        ) = mask_batch_with_augmentation(
            batch,
            pad_index=config['pad_token_id'],
            mask_index=config['mask_token_id'],
            unk_index=config['unk_token_id'],
            mlm_sample_ratio=config['mlm_sample_ratio'],
            token_mask_ratio=config['token_mask_ratio'],
            spectrum_topk=config['spectrum_topk'],
            spectrum_frac_shuffle=config['spectrum_frac_shuffle'],
            spectrum_num_shuffle=config['spectrum_num_shuffle'],
            spectrum_mask_ratio=config['spectrum_mask_ratio'],
            spectrum_drop_ratio=config['spectrum_drop_ratio'],
            spectrum_augment=config.get('spectrum_augment', True),
            device=local_device,
        )

        if target_dtype is not None:
            spectra = spectra.to(target_dtype)
            precursors = precursors.to(target_dtype)
            labels = labels.to(target_dtype)
            weights = weights.to(target_dtype)

        tokens_input = tokens_input.long()
        tokens_label = tokens_label.long()
        token_mask = token_mask.bool()
        spectra_mask = spectra_mask.bool()

        dda_pred, mask_pred, mask_positions = model_engine(
            spectra,
            spectra_mask,
            precursors,
            tokens_input,
            token_mask,
        )

        dda_criterion = nn.BCEWithLogitsLoss(weight=weights)
        dda_loss = dda_criterion(dda_pred, labels.flatten())

        target_mask = (labels > 0.5)
        target_dda_loss = calc_loss(target_mask, weights, labels, dda_pred)

        ft_part1_mask = (labels < 0.5) & (weights < 0.9)
        ft_part1_dda_loss = calc_loss(ft_part1_mask, weights, labels, dda_pred)

        ft_part2_mask = (labels < 0.5) & (weights >= 0.9)
        ft_part2_dda_loss = calc_loss(ft_part2_mask, weights, labels, dda_pred)

        mlm_loss = compute_mlm_loss_with_positions(
            mask_pred,
            mask_positions,
            tokens_label,
            token_mask,
            label_smoothing=config.get('label_smoothing', 0.0),
        )

        loss = (
            config['dda_loss_weight'] * dda_loss +
            config['mlm_loss_weight'] * mlm_loss
        )

        try:
            model_engine.backward(loss)
            model_engine.step()
        except Exception as e:
            logging.info(f"epoch: {epoch}, rank: {int(os.environ['RANK'])}, error: {e}!!!")
            raise

        if running_loss is None:
            running_loss = loss.item()
            dda_running_loss = dda_loss.item()
            mask_running_loss = mlm_loss.item()

            target_running_loss = target_dda_loss
            ft_part1_running_loss = ft_part1_dda_loss
            ft_part2_running_loss = ft_part2_dda_loss
        else:
            running_loss = 0.99 * running_loss + (1 - 0.99) * loss.item()
            dda_running_loss = 0.99 * dda_running_loss + (1 - 0.99) * dda_loss.item()
            mask_running_loss = 0.99 * mask_running_loss + (1 - 0.99) * mlm_loss.item()

            target_running_loss = 0.99 * target_running_loss + (1 - 0.99) * target_dda_loss
            ft_part1_running_loss = 0.99 * ft_part1_running_loss + (1 - 0.99) * ft_part1_dda_loss
            ft_part2_running_loss = 0.99 * ft_part2_running_loss + (1 - 0.99) * ft_part2_dda_loss

        if is_rank0():
            sw_step += 1

            if sw_step % int(config["one_epoch_iters"] * 0.5) == 0:
                if hasattr(model_engine, 'lr_scheduler'):
                    scheduler_lr = float(model_engine.lr_scheduler.get_last_lr()[0])
                else:
                    scheduler_lr = float(optim.param_groups[0]['lr'])

                logging.info(f"[Epoch={epoch}, sw_step={sw_step}, rank={int(os.environ['RANK'])}, loss={running_loss}, lr={scheduler_lr}]")

                sw.add_scalar("train/loss", loss.item(), sw_step)
                sw.add_scalar("train/loss_smooth", running_loss, sw_step)

                sw.add_scalar("train/dda_loss", dda_loss.item(), sw_step)
                sw.add_scalar("train/dda_loss_smooth", dda_running_loss, sw_step)

                sw.add_scalar("train/target_loss", target_running_loss, sw_step)
                sw.add_scalar("train/ft1_loss", ft_part1_running_loss, sw_step)
                sw.add_scalar("train/ft2_loss", ft_part2_running_loss, sw_step)

                sw.add_scalar("train/mlm_loss", mlm_loss.item(), sw_step)
                sw.add_scalar("train/mlm_loss_smooth", mask_running_loss, sw_step)

                sw.add_scalar("optim/scheduler_lr", scheduler_lr, sw_step)
                sw.add_scalar("optim/epoch", epoch, sw_step)
                sw.flush()

        global_step = model_engine.global_steps

        if (global_step > 0) and (global_step % config["ckpt_interval"] == 0):
            logging.info(
                f"[rank={get_rank()}] saving ckpt at global_step={global_step}, sw_step={sw_step}"
            )
            save_checkpoint(
                model_engine,
                save_dir=config["model_save_path"],
                tag=epoch,
                client_state={'epoch': epoch, 'sw_step': int(global_step)},
            )
            logging.info(
                f"[rank={get_rank()}] finished ckpt at global_step={global_step}, sw_step={sw_step}"
            )

    if is_rank0() and running_loss is not None:
        sw.add_scalar("train/train_loss", running_loss, epoch)
    return sw_step


def auto_find_latest_tag(ckpt_dir: str):
    if not os.path.exists(ckpt_dir):
        return None
    tags = []
    for name in os.listdir(ckpt_dir):
        full = os.path.join(ckpt_dir, name)
        if os.path.isdir(full):
            try:
                tags.append(int(name))
            except ValueError:
                continue
    if not tags:
        return None
    return str(max(tags))


def main():
    deepspeed.init_distributed(timeout=datetime.timedelta(seconds=5400))
    torch.cuda.set_device(args.local_rank)

    config_path = args.config
    with open(config_path) as f_in:
        config = yaml.safe_load(f_in)

    if torch.distributed.get_rank() != 0:
        torch.distributed.barrier()

    vocab = ['<pad>', '<mask>'] + list(config["residues"].keys()) + ['<unk>']
    config["vocab"] = vocab
    config['pad_token_id'] = 0
    config['mask_token_id'] = 1
    config['unk_token_id'] = 28
    config["node_num"] = args.node_num
    config["gpu_num"] = args.gpu_num

    set_seeds(config['seed'])

    s2i = {v: i for i, v in enumerate(vocab)}
    logging.info(f"Vocab: {s2i}")

    mkdir_p(config["tb_summarywriter"])
    mkdir_p(config["model_save_path"])

    if torch.distributed.get_rank() == 0:
        log_dir = config["tb_summarywriter"] + datetime.datetime.now().strftime("DDABert_%y_%m_%d_%H_%M")
        sw = SummaryWriter(log_dir)
    else:
        sw = SummaryWriter(log_dir="/tmp/ddabert_dummy")

    if torch.distributed.get_rank() == 0:
        torch.distributed.barrier()

    multi_node = args.node_num > 1
    train_dl = create_iterable_dataset(logging, config, s2i, parse='train', multi_node=multi_node)

    model = DDABert(
        dim_model=config['dim_model'],
        n_head=config['n_head'],
        dim_feedforward=config['dim_feedforward'],
        n_layers=config['n_layers'],
        dropout=config['dropout'],
        max_length=config['max_length'],
        vocab=vocab,
        max_charge=config['max_charge'],
    )

    if (config.get('init_model_path') is not None) and config['init_model_path'] != '':
        logging.info(f"Loading model checkpoint from '{config['init_model_path']}'")
        model = DDABert.load_pt(config['init_model_path'], config)
        model = model.to(torch.bfloat16)

    parameters = filter(lambda p: p.requires_grad, model.parameters())

    base_optimizer = torch.optim.Adam(
        model.parameters(),
        lr=float(config["learning_rate"]),
        weight_decay=float(config["weight_decay"]),
    )

    one_epoch_iters = max(len(train_dl), 1)
    train_step_ratio = float(config.get("train_step_ratio", 1.0))
    one_epoch_iters = max(int(one_epoch_iters * train_step_ratio), 1)

    total_iters = config["epochs"] * one_epoch_iters
    phase1_max_iter = max(1, int(total_iters * 0.5))

    warmup_ratio = float(config.get("warmup_ratio", 0.2))
    config["one_epoch_iters"] = one_epoch_iters
    config["ckpt_interval"] = int(one_epoch_iters * float(config.get("ckpt_step_ratio", 0.9)))

    logging.info(
        f"one_epoch_iters={one_epoch_iters}, total_iters={total_iters}, "
        f"phase1_max_iter={phase1_max_iter}, warmup_ratio={warmup_ratio}, ckpt_interval={config['ckpt_interval']}"
    )

    model_engine, ds_optimizer, _, _ = deepspeed.initialize(
        args=args,
        model=model,
        model_parameters=parameters,
        optimizer=base_optimizer,
        lr_scheduler=None,
    )

    scheduler = TwoPhaseScheduler(
        optimizer=ds_optimizer,
        phase1_max_iter=phase1_max_iter,
        warmup_ratio=warmup_ratio,
        max_lr=float(config['learning_rate']),
        min_lr=float(config['min_lr']),
        second_min_lr=float(config['second_min_lr']),
        warmup_type=config.get('warmup_strategy', 'cos'),
    )
    model_engine.lr_scheduler = scheduler

    local_device = get_accelerator().device_name(model_engine.local_rank)
    local_rank = model_engine.local_rank

    target_dtype = None
    if model_engine.bfloat16_enabled():
        target_dtype = torch.bfloat16
    elif model_engine.fp16_enabled():
        target_dtype = torch.half
    logging.info(f"target_dtype: {target_dtype}")

    # ===========================
    # Resume logic
    # ===========================
    start_epoch = 0
    sw_step = 0

    if args.resume:
        ckpt_dir = config["resume_model_path"]
        tag = args.resume_tag
        if tag is None:
            tag = auto_find_latest_tag(ckpt_dir)
        if tag is None:
            logging.warning(f"No checkpoint found in {ckpt_dir}, start from scratch.")
        else:
            logging.info(f"Resuming from checkpoint tag={tag} in {ckpt_dir}")
            load_path, client_state = model_engine.load_checkpoint(ckpt_dir, tag=str(tag))
            if load_path is None:
                logging.warning("DeepSpeed load_checkpoint returned None, start from scratch.")
            else:
                start_epoch = client_state.get('epoch', 0) + 1
                sw_step = client_state.get('sw_step', 0)
                logging.info(f"Resume from epoch={start_epoch}, sw_step={sw_step}")

    # ===========================
    # Train loop
    # ===========================
    for epoch in range(start_epoch, config['epochs']):
        try:
            sw_step = train(
                model_engine,
                train_dl,
                sw,
                ds_optimizer,
                scheduler,
                local_device,
                target_dtype,
                local_rank,
                config,
                epoch,
                sw_step,
            )
        except Exception as e:
            logging.info(f"epoch: {epoch}, rank: {int(os.environ['RANK'])}, error: {e}!!!")
            continue

    logging.info('Finished Training')


if __name__ == '__main__':
    main()
