import logging
import os
import yaml
import random
import shutil
import glob
import numpy as np
import pandas as pd

import torch
from torch import nn
from sklearn import metrics
from collections import defaultdict


def set_seeds(seed):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        # if you are using multi-GPU.
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    # Numpy module.
    np.random.seed(seed)
    # Python random module.
    random.seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


def mkdir_p(dirs):
    """
    make a directory (dir) if it doesn't exist
    """
    if not os.path.exists(dirs):
        try:
            os.makedirs(dirs)
        except:
            pass

    return True, 'OK'


def eval_predict(predictions, targets):
    fpr, tpr, threshold = metrics.roc_curve(targets, predictions, pos_label=1)
    auc_results = metrics.auc(fpr, tpr)
    round_pred = np.round(predictions)
    correct_count = 0
    for i in range(len(round_pred)):
        if round_pred[i] == targets[i]:
            correct_count += 1
    correctness = float(correct_count) / len(round_pred)
    return auc_results, correctness