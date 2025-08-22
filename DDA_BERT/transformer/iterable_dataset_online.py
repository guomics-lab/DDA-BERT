from __future__ import annotations

import torch
import torch.utils.data as pt_data
import torch.nn.functional as F
from torch.utils.data import IterableDataset, DataLoader, ConcatDataset
from torch import Tensor
from torch.utils.data import random_split

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


def shuffle_file_list(file_list, seed):
    generator = torch.Generator()
    generator.manual_seed(seed)
    idx = torch.randperm(len(file_list), generator=generator).numpy()
    file_list = (np.array(file_list)[idx]).tolist()
    return file_list


class Dataset_weight(pt_data.Dataset):
    def __init__(self):
        self.spectra = None
        self.spectra_mask = None
        self.precursors = None
        self.tokens = None
        self.label = None
        self.weight = None

    def __getitem__(self, idx):
        return_dict = {"spectra": self.spectra[idx],
                       "spectra_mask": self.spectra_mask[idx],
                       "precursors": self.precursors[idx],
                       "tokens": self.tokens[idx],
                       "label": self.label[idx],
                       "weight": self.weight[idx]}
        return return_dict

    def __len__(self):
        return len(self.label)

    def fit_scale(self):
        pass


def dict2pt_weight(
    batch: dict,
    s2i: dict,
    max_length: int,
) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor, Tensor]:
    """dict to pt_data.Dataset."""
    # gen dataset
    data_set = Dataset_weight()
    data_set.spectra = np.nan_to_num(batch['spectra'])
    data_set.spectra_mask = np.nan_to_num(batch['spectra_mask'])
    data_set.precursors = np.nan_to_num(batch['precursors'])
    data_set.tokens = np.nan_to_num(batch['tokens'])
    data_set.weight = np.nan_to_num(batch['weight'])
    data_set.label = np.nan_to_num(batch['label'])
    return data_set


def collate_numpy_batch_weight(batch_data):
    """Collate batch of samples."""
    one_batch_spectra = torch.tensor(np.array([batch["spectra"] for batch in batch_data]), dtype=torch.float)
    one_batch_spectra_mask = torch.tensor(np.array([batch["spectra_mask"] for batch in batch_data]), dtype=torch.float)
    one_batch_precursors = torch.tensor(np.array([batch["precursors"] for batch in batch_data]), dtype=torch.float)
    one_batch_tokens = torch.tensor(np.array([batch["tokens"] for batch in batch_data]), dtype=torch.float)
    one_batch_label = torch.tensor(np.array([batch["label"] for batch in batch_data]), dtype=torch.float)
    one_batch_weight = torch.tensor(np.array([batch["weight"] for batch in batch_data]), dtype=torch.float)

    return one_batch_spectra, one_batch_spectra_mask, one_batch_precursors, one_batch_tokens, one_batch_label, one_batch_weight


def create_iterable_dataset(logging,
                            config,
                            s2i,
                            parse='train',
                            multi_node=False,
                            seed=123):
    """
    Note: If you want to load all data in the memory, please set "read_part" to False.
    Args:
        :param logging: out logging.
        :param config: data from the yaml file.
        :param s2i: vocab.
        :param buffer_size: An integer. the size of file_name buffer.
    :return:
    """
    # update gpu_num
    if multi_node:
        gpu_num = int(config['gpu_num'])
    else:
        gpu_num = torch.cuda.device_count() if torch.cuda.is_available() else 1
    logging.info(f"******************multi_node: {multi_node}, need_augment: {need_augment}, gpu_num: {gpu_num};**********")

    if parse == 'train':
        if ';' in config['train_path']:
            total_train_path = config['train_path'].split(';')
            data_file_list = []
            for train_path in total_train_path:
                train_part_file_list = glob.glob(f'{train_path}/*.pkl')
                if len(train_part_file_list) > 0:
                    data_file_list.extend(train_part_file_list)
            logging.info(f"******************{parse} {config['train_path']}, total loaded: {len(data_file_list)};**********")
        else:
            data_file_list = glob.glob(f"{config['train_path']}/*pkl")
            logging.info(f"******************{parse} {config['train_path']}, loaded: {len(data_file_list)};**********")
        
        random.shuffle(data_file_list)
        data_file_list = shuffle_file_list(data_file_list, config['seed'])
        
        # Truncate the dataset according to the number of GPUs so that each GPU processes the same number of samples
        file_bin_num = len(data_file_list) // gpu_num
        file_truncation_num = file_bin_num * gpu_num
        data_file_list = data_file_list[:file_truncation_num]

        train_dl = IterableDiartDataset(
            data_file_list,
            config["train_batch_size"],
            s2i,
            max_length=config["max_length"],
            buffer_size=config["buffer_size"],
            gpu_num=gpu_num,
            shuffle=True,
            multi_node=multi_node,
            seed=config['seed']
        )
        logging.info(f"Data loaded: {len(train_dl) * config['train_batch_size']:,} training samples")
        return train_dl
    else:
        data_file_list = glob.glob(f"{config['valid_path']}/*pkl")
        logging.info(f"******************{parse} loaded: {len(data_file_list)};**********")
        
        val_dl = IterableDiartDataset(
            data_file_list,
            config["predict_batch_size"],
            s2i, # vocab
            max_length=config["max_length"],
            gpu_num=gpu_num,
            shuffle=False,
            multi_node=multi_node
        )
        logging.info(f"{len(val_dl) * config['predict_batch_size']:,} validation samples")
        return val_dl


class IterableDiartDataset(IterableDataset):
    """
    Custom dataset class for dataset in order to use efficient
    dataloader tool provided by PyTorch.
    """
    def __init__(self,
                 file_list: list,
                 batch_size,
                 s2i, # vocab
                 max_length=50,
                 buffer_size=1,
                 gpu_num=1,
                 shuffle=False,
                 multi_node=False,
                 seed=0,
                 bath_file_size=1,
                 **kwargs):
        super(IterableDiartDataset).__init__()
        self.file_list = file_list
        self.batch_size = batch_size
        self.s2i = s2i
        self.max_length = max_length

        self.shuffle = shuffle
        self.seed = seed
        self.epoch = 0
        
        self.bath_file_size = bath_file_size
        self.buffer_size = buffer_size

        self.gpu_num = gpu_num
        self.multi_node = multi_node
        
    def parse_file(self, file_name):
        ds = pd.read_pickle(file_name)
        dpt = dict2pt_weight(ds, self.s2i, self.max_length)
        return DataLoader(dpt,
                          batch_size=self.batch_size,
                          shuffle=self.shuffle,
                          collate_fn=collate_numpy_batch_weight,
                          num_workers=2,
                          drop_last=True,
                          pin_memory=True)

    def file_mapper(self, file_list):
        idx = 0
        file_num = len(file_list)
        while idx < file_num:
            yield self.parse_file(file_list[idx])
            idx += 1

    def __iter__(self):
        if self.gpu_num > 1:
            if self.multi_node:
                if 'RANK' in os.environ:
                    rank = int(os.environ['RANK'])
                else:
                    rank = 0
            else:
                if 'LOCAL_RANK' in os.environ:
                    rank = int(os.environ['LOCAL_RANK'])
                else:
                    rank = 0
            
            file_itr = self.file_list[rank::self.gpu_num]
        else:
            file_itr = self.file_list

        file_mapped_itr = self.file_mapper(file_itr)

        if self.shuffle:
            return self._shuffle(file_mapped_itr)
        else:
            return file_mapped_itr

    def __len__(self):
        if self.gpu_num > 1:
            return math.ceil(len(self.file_list) / self.gpu_num)
        else:
            return len(self.file_list)
        
    def set_epoch(self, epoch):
        self.epoch = epoch

    def generate_random_num(self):
        while True:
            random.seed(self.epoch + self.seed)
            random_nums = random.sample(range(self.buffer_size), self.bath_file_size)
            yield from random_nums

    def _shuffle(self, mapped_itr):
        buffer = []
        for dt in mapped_itr:
            if len(buffer) < self.buffer_size:
                buffer.append(dt)
            else:
                i = next(self.generate_random_num())
                yield buffer[i]
                buffer[i] = dt
        random.shuffle(buffer)
        yield from buffer


def aa_tokenize(sequence, s2i, max_length):
    """Transform a peptide sequence into tokens

    Parameters
    ----------
    sequence : str
        A peptide sequence.
    s2i : str
        amino acid s2i.
    max_length : int
        default is 50.

    Returns
    -------
    torch.Tensor
        The token for each amino acid in the peptide sequence.
    """
    sequence = sequence.replace("I", "L").replace('n[42]', 'X')
    sequence = re.split(r"(?<=.)(?=[A-Z])", sequence)

    tokens = torch.tensor([s2i[aa] for aa in sequence])
    tokens = F.pad(tokens, (0, max_length - tokens.shape[0]), 'constant', 0)
    return tokens


def mask_target_decoy_1unk_tokens(tokens, label, mask_index, unk_index, token_mask_ratio=0.15, device='cpu'):
    """Args:
            tokens: float Tensor (batch, 50)
            label: float Tensor (batch)
            mask_index: int=1
            unk_index: int=28
            token_mask_ratio: float=0.15
            device: cpu or cuda
        Returns:
            tokens_mask: float Tensor (batch, 50)
            tokens_label: float Tensor (batch)
    """
    mask = torch.rand(tokens.shape[0], tokens.shape[1]).to(device)
    mask = torch.where(mask < token_mask_ratio, 1, 0)
    # Set mask of first column (index 0) to 0, i.e., no masking.
    mask[:, 0] = 0
    
    tokens = tokens.to(device)
    tokens_mask = torch.where((mask==1) & (tokens > mask_index), mask_index, tokens).to(device)

    label_adjust = torch.where(label==0, -1, label).to(device)
    label_mask = torch.mul(mask, label_adjust.reshape(label_adjust.shape[0], -1)).to(torch.int8).to(device)

    # Keep one per row in decoy, set others to unk
    decoy_tokens = mask_decoy_1unk(tokens, unk_index, device)
    
    # target: predict masked amino acids
    # decoy: predict fixed placeholders
    tokens_label = torch.where((label_mask==1) & (tokens > mask_index), tokens, \
                   torch.where((label_mask==-1) & (tokens > mask_index), decoy_tokens, 0)).to(device)
    return tokens_mask, tokens_label


def mask_decoy_1unk(tokens, unk_index=28, device='cpu'):
    num_rows, num_cols = tokens.shape
    random_indices = torch.randint(0, num_cols, (num_rows,)).to(device)

    mask = torch.zeros_like(tokens, dtype=torch.bool).to(device)
    mask[torch.arange(num_rows), random_indices] = True

    result = torch.where(mask, tokens, torch.tensor(unk_index)).to(device)
    return result


def mask_batch_decoy_1unk_data(batch, pad_index=0, mask_index=1, unk_index=28, mask_ratio=0.9, token_mask_ratio=0.15, device='cpu'):
    """Args:
            batch: spectra, spectra_mask, precursors, tokens, label, weight. the dim is 
            float Tensor (batch, 300, 2), float Tensor (batch, 300), float Tensor (batch, 3), float Tensor (batch, 50), float Tensor (batch)
            pad_index: int=0
            mask_index: int=1
            unk_index: int=28
            mask_ratio: float=0.9
            token_mask_ratio: float=0.15
            device: cpu or cuda
        Returns:
            batch: spectra, spectra_mask, precursors, tokens, label. 
    """
    spectra, spectra_mask, precursors, tokens, label, weight = batch

    full_indices = list(range(label.shape[0]))
    mask_size = int(mask_ratio * len(full_indices))
    remain_size = len(full_indices) - mask_size
    mask_index_list, remain_index_list = random_split(full_indices, [mask_size, remain_size])

    mask_spectra = spectra[mask_index_list].to(device)
    mask_spectra_mask = spectra_mask[mask_index_list].to(device)
    mask_precursors = precursors[mask_index_list].to(device)
    mask_tokens = tokens[mask_index_list].to(device)
    mask_label = label[mask_index_list].to(device)
    mask_weight = weight[mask_index_list].to(device)

    remain_spectra = spectra[remain_index_list].to(device)
    remain_spectra_mask = spectra_mask[remain_index_list].to(device)
    remain_precursors = precursors[remain_index_list].to(device)
    remain_tokens = tokens[remain_index_list].to(device)
    remain_label = label[remain_index_list].to(device)
    remain_weight = weight[remain_index_list].to(device)
    
    tokens_masked, tokens_label_masked = mask_target_decoy_1unk_tokens(mask_tokens,
                                                                       mask_label,
                                                                       mask_index,
                                                                       unk_index,
                                                                       token_mask_ratio=token_mask_ratio,
                                                                       device=device)
    remain_tokens_label = torch.zeros(remain_tokens.shape[0], remain_tokens.shape[1]).to(device)

    spectra = torch.cat((mask_spectra, remain_spectra))
    spectra_mask = torch.cat((mask_spectra_mask, remain_spectra_mask))
    precursors = torch.cat((mask_precursors, remain_precursors))
    tokens = torch.cat((tokens_masked, remain_tokens))
    tokens_label = torch.cat((tokens_label_masked, remain_tokens_label))
    label = torch.cat((mask_label, remain_label))
    weight = torch.cat((mask_weight, remain_weight))
    return spectra, spectra_mask, precursors, tokens, tokens_label, label, weight


def mask_spectra_data(spectra, spectra_mask, remain_ratio=0.1, spectra_zero_ratio=0.1, device='cpu'):
    """Args:
            spectra: float Tensor (batch, 300, 2)
            spectra_mask: float Tensor (batch, 300)
            remain_ratio: float=0.1
            spectra_zero_ratio: float=0.1
            device: cpu or cuda
        Returns:
            spectra: float Tensor (batch, 300, 2)
            spectra_mask: float Tensor (batch, 300)
    """
    spectra_mask_backup = torch.rand(spectra.shape[0], spectra.shape[1]).to(device)
    spectra_mask_backup = torch.where(spectra_mask_backup < spectra_zero_ratio, 0, 1)

    full_indices = list(range(spectra.shape[0]))
    remain_size = int(remain_ratio * len(full_indices))
    mask_size = len(spectra) - remain_size
    remain_spectra_index_list, mask_spectra_index_list = random_split(full_indices, [remain_size, mask_size])

    for i in remain_spectra_index_list:
        spectra_mask_backup[i] = 1

    spectra = torch.mul(spectra, torch.unsqueeze(spectra_mask_backup, dim=-1)).to(device)
    spectra_mask = torch.mul(spectra_mask, spectra_mask_backup).to(device)
    return spectra, spectra_mask