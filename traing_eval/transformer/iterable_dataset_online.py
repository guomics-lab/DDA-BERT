from __future__ import annotations

import torch
import torch.utils.data as pt_data
import torch.nn.functional as F
from torch.utils.data import IterableDataset, DataLoader, ConcatDataset, random_split
from torch import Tensor
from torch.utils.data import random_split
from typing import List, Dict, Any, Iterator

import os
import numpy as np
import polars as pl
import pickle
import random
from random import sample
import math
import re
from collections import defaultdict
import glob
import time


class Dataset(pt_data.Dataset):
    def __init__(self):
        self.spectra = None
        self.spectra_mask = None
        self.precursors = None
        self.tokens = None
        self.label = None
        self.weight = None

    def __getitem__(self, idx):
        return {
            "spectra": self.spectra[idx],
            "spectra_mask": self.spectra_mask[idx],
            "precursors": self.precursors[idx],
            "tokens": self.tokens[idx],
            "label": self.label[idx],
            "weight": self.weight[idx],
        }

    def __len__(self):
        return len(self.label)

    def fit_scale(self):
        pass


def dict2pt(batch: Dict[str, Any]) -> Dataset:
    """Convert a dict loaded from a pkl file into a PyTorch Dataset."""
    data_set = Dataset()
    data_set.spectra = np.nan_to_num(batch["spectra"])
    data_set.spectra_mask = np.nan_to_num(batch["spectra_mask"])
    data_set.precursors = np.nan_to_num(batch["precursors"])
    data_set.tokens = np.nan_to_num(batch["tokens"])
    data_set.weight = np.nan_to_num(batch["weight"])
    data_set.label = np.nan_to_num(batch["label"])
    return data_set


def collate_batch(batch_data):
    """Collate a batch of samples."""
    one_batch_spectra = torch.tensor(np.array([batch["spectra"] for batch in batch_data]), dtype=torch.float)
    one_batch_spectra_mask = torch.tensor(np.array([batch["spectra_mask"] for batch in batch_data]), dtype=torch.float)
    one_batch_precursors = torch.tensor(np.array([batch["precursors"] for batch in batch_data]), dtype=torch.float)
    one_batch_tokens = torch.tensor(np.array([batch["tokens"] for batch in batch_data]), dtype=torch.float)
    one_batch_label = torch.tensor(np.array([batch["label"] for batch in batch_data]), dtype=torch.float)
    one_batch_weight = torch.tensor(np.array([batch["weight"] for batch in batch_data]), dtype=torch.float)

    return (
        one_batch_spectra,
        one_batch_spectra_mask,
        one_batch_precursors,
        one_batch_tokens,
        one_batch_label,
        one_batch_weight,
    )


def shuffle_file_list(file_list: List[str], seed: int) -> List[str]:
    generator = torch.Generator()
    generator.manual_seed(seed)
    idx = torch.randperm(len(file_list), generator=generator).numpy()
    file_list = (np.array(file_list)[idx]).tolist()
    return file_list


def create_iterable_dataset(
    logging,
    config,
    s2i,
    parse: str = "train",
    multi_node: bool = False,
    seed: int = 123,
):
    """
    Streaming IterableDataset based on a global-batch stream sharded by global_step.

    - Iterate over all pkl files and collect the number of samples per pkl
    - Remove pkls with too few samples
    - Remove pkls whose sample counts deviate too much from the "mainstream" size
    - Truncate the file list to a multiple of the GPU count (world_size or gpu_num) to keep training
      as balanced as possible across GPUs
    """
    # -------- 1. Load file list --------
    if parse == "train":
        if ";" in config["train_path"]:
            total_train_path = config["train_path"].split(";")
            data_file_list = []
            for train_path in total_train_path:
                train_part_file_list = glob.glob(f"{train_path}/*.pkl")
                if len(train_part_file_list) > 0:
                    data_file_list.extend(train_part_file_list)
            logging.info(
                f"******************{parse} {config['train_path']}, total loaded: {len(data_file_list)};**********"
            )
        else:
            data_file_list = glob.glob(f"{config['train_path']}/*pkl")
            logging.info(
                f"******************{parse} {config['train_path']}, loaded: {len(data_file_list)};**********"
            )
    else:
        data_file_list = glob.glob(f"{config['valid_path']}/*pkl")
        logging.info(f"******************{parse} loaded: {len(data_file_list)};**********")

    if len(data_file_list) == 0:
        logging.error(f"No pkl files found for parse={parse}")
        raise RuntimeError("Empty dataset")

    # -------- 2. Count the number of samples in each pkl --------
    file_sizes = []
    for fpath in data_file_list:
        try:
            with open(fpath, "rb") as f:
                ds = pickle.loads(f.read())
            # Assume 'label' is the sample dimension
            n = len(ds["label"])
            file_sizes.append((fpath, n))
        except Exception as e:
            logging.warning(f"Failed to read {fpath}: {e}")
            continue

    if not file_sizes:
        logging.error("All pkl files failed to load.")
        raise RuntimeError("All pkl files failed to load")

    # -------- 3. Remove pkls that are too small or inconsistent with the mainstream size --------
    sizes = np.array([n for _, n in file_sizes], dtype=np.int64)

    # 3.1 Filter out files with too few samples (e.g., < min_size)
    min_size = int(config.get("min_pkl_size", 131072))  # configurable in yaml; default is 128k samples
    filtered = [(f, n) for (f, n) in file_sizes if n >= min_size]
    if not filtered:
        logging.warning(
            f"All pkls smaller than min_pkl_size={min_size}, fallback to all files."
        )
        filtered = file_sizes

    sizes = np.array([n for _, n in filtered], dtype=np.int64)

    # 3.2 Use the median as the "mainstream size" and drop files that deviate too much
    median_size = float(np.median(sizes))
    # Allowed relative deviation, e.g., 50%
    size_tol = float(config.get("pkl_size_tolerance", 0.6))  # 0.6 -> ±60%
    low_bound = median_size * (1.0 - size_tol)
    high_bound = median_size * (1.0 + size_tol)

    balanced = [
        (f, n) for (f, n) in filtered if (n >= low_bound and n <= high_bound)
    ]
    if not balanced:
        logging.warning(
            f"No pkls within [{low_bound:.1f}, {high_bound:.1f}], "
            f"fallback to min_size-filtered list."
        )
        balanced = filtered

    # Final file list for training/validation
    data_file_list = [f for (f, _) in balanced]
    logging.info(
        f"After filtering: {len(data_file_list)} pkls kept "
        f"(median_size={median_size:.1f}, tol={size_tol})"
    )

    sizes = np.array([n for _, n in balanced], dtype=np.int64)
    total_samples = sizes.sum()
    batch_size = config["train_batch_size"] if parse == "train" else config["predict_batch_size"]
    approx_total_batches = int(total_samples // batch_size)

    logging.info(
        f"approx_total_batches: {approx_total_batches} "
        f"(batch_size={batch_size}!!!!!!!!!!)"
    )

    # -------- 4. Truncate by GPU count to balance data across GPUs --------
    # For global-batch streaming + sharding, we prefer:
    # - total number of files to be a multiple of world_size so that each rank sees a similar
    #   number of files per epoch (sharding is actually done on the batch dimension; file-level
    #   balancing is auxiliary)
    if multi_node and torch.distributed.is_available() and torch.distributed.is_initialized():
        world_size = torch.distributed.get_world_size()
    else:
        # Single-node case: use config['gpu_num'] as an approximation of world_size
        world_size = int(config.get("gpu_num", 1))
        world_size = max(world_size, 1)

    if len(data_file_list) < world_size:
        logging.warning(
            f"Number of pkls ({len(data_file_list)}) < world_size ({world_size}), "
            f"some ranks may see fewer batches."
        )
    else:
        # Truncate to a multiple of world_size
        trunc_len = (len(data_file_list) // world_size) * world_size
        if trunc_len <= 0:
            trunc_len = world_size  # ensure at least one file per rank (if enough files exist)
        if trunc_len < len(data_file_list):
            logging.info(
                f"Truncating pkls from {len(data_file_list)} to {trunc_len} "
                f"to better balance across {world_size} GPUs."
            )
            data_file_list = data_file_list[:trunc_len]

    # -------- 5. Shuffle file order --------
    random.shuffle(data_file_list)
    data_file_list = shuffle_file_list(data_file_list, seed)

    # -------- 6. Build GlobalShardIterableDataset --------
    if parse == "train":
        dl = GlobalShardIterableDataset(
            file_list=data_file_list,
            batch_size=config["train_batch_size"],
            shuffle=True,
            seed=seed,
            approx_total_batches=approx_total_batches,
        )
        logging.info(
            f"Data loaded: ~{len(dl) * config['train_batch_size']:,} training samples (approx)"
        )
        return dl
    else:
        dl = GlobalShardIterableDataset(
            file_list=data_file_list,
            batch_size=config["predict_batch_size"],
            shuffle=False,
            seed=seed,
            approx_total_batches=approx_total_batches,
        )
        logging.info(
            f"~{len(dl) * config['predict_batch_size']:,} validation samples (approx)"
        )
        return dl


class GlobalShardIterableDataset(IterableDataset):
    """
    Global-batch stream sharded by global_step:
    - In both single-node and multi-node settings, every rank iterates over the same file list
      and the same global batch stream order.
    - Batches are assigned by: global_step % world_size == rank
    """

    def __init__(
        self,
        file_list: List[str],
        batch_size: int,
        shuffle: bool = False,
        seed: int = 0,
        approx_total_batches: int | None = None,
        **kwargs,
    ):
        super(GlobalShardIterableDataset, self).__init__()
        self.file_list = file_list
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.seed = seed
        self.epoch = 0
        self.approx_total_batches = approx_total_batches

    def set_epoch(self, epoch: int):
        self.epoch = epoch

    def _get_rank_and_world_size(self) -> Tuple[int, int]:
        if torch.distributed.is_available() and torch.distributed.is_initialized():
            rank = torch.distributed.get_rank()
            world_size = torch.distributed.get_world_size()
        else:
            rank = 0
            world_size = 1
        return rank, world_size

    def _iter_one_file(self, file_name: str) -> Iterator:
        with open(file_name, "rb") as f:
            ds = pickle.loads(f.read())

        dpt = dict2pt(ds)

        loader = DataLoader(
            dpt,
            batch_size=self.batch_size,
            shuffle=self.shuffle,
            collate_fn=collate_batch,
            num_workers=0,
            drop_last=True,
            pin_memory=True,
        )

        for batch in loader:
            yield batch

    def _global_batch_stream(self) -> Iterator:
        """Global batch stream from a single-process perspective: iterate over all files and yield batches sequentially."""
        file_list = list(self.file_list)
        if self.shuffle:
            rng = random.Random(self.seed + self.epoch)
            rng.shuffle(file_list)

        for fname in file_list:
            for batch in self._iter_one_file(fname):
                yield batch

    def __iter__(self):
        rank, world_size = self._get_rank_and_world_size()
        global_step = 0

        for batch in self._global_batch_stream():
            if global_step % world_size == rank:
                yield batch
            global_step += 1

    def __len__(self):
        _, world_size = self._get_rank_and_world_size()
        if self.approx_total_batches is not None:
            return max(1, self.approx_total_batches // max(world_size, 1))
        # fallback: always return a positive value
        return max(1, len(self.file_list) // max(world_size, 1))
