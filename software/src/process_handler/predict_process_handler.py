from __future__ import annotations

import json
import os
import shutil
import time
from copy import deepcopy
from multiprocessing import Process

import lightning.pytorch as ptl
import pandas as pd
import polars as pl
import torch
import yaml
from lightning.pytorch.strategies import DDPStrategy
from torch.utils.data import DataLoader

from src.common import constant
from src.common.eval_model import Evalute
from src.common.runtime_data_info import runtime_data
from src.obj.Entity import MzMLFileInfo, ProcessMsg
from src.process_handler.common_process_handler import CommonProcessHandler
from src.transformer.dataset import SpectrumDataset
from src.transformer.dataset import collate_batch_weight_deltaRT_index
from src.transformer.model import DDA_BERT


class PredictProcessHandler(CommonProcessHandler):

    def __init__(self, f_info: MzMLFileInfo, base_output_path, model_path, logger, env, input_param, start_time=0):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.PREDICT)
        self.model_path = model_path
        self.base_config = 'config/model.yaml'

    def deal_process(self):
        self.send_msg('Processing predict data', status=constant.ProgressStepStatusEnum.RUNNING)
        pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.PREDICT,
                        constant.ProgressStepStatusEnum.RUNNING)
        print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
        mzml_info = self.f_info
        try:
            ss = time.time()
            if not runtime_data.current_is_success:
                self.send_msg(status=constant.ProgressStepStatusEnum.ERROR)
                return
            if not self.input_param.open_pred:
                ee = time.time()
                self.send_msg('Finished predict data, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.SUCCESS)
                return
            ss = time.time()
            self.send_msg('Load config')
            config = self.load_config()
            config['model_path'] = self.model_path

            self.send_msg(
                'Processing predict file {}'.format(mzml_info.base_file_name))
            self.deal_one(config, mzml_info)
            ee = time.time()
            self.send_msg('Finished predict data, spend time: {}'.format(self.get_use_time(ss, ee)),
                          status=constant.ProgressStepStatusEnum.SUCCESS)
        except Exception as e:
            runtime_data.current_is_success = False
            self.logger.exception('Predict data exception')
            self.send_msg('Predict data exception: {}'.format(e), status=constant.ProgressStepStatusEnum.ERROR)
        finally:
            if runtime_data.current_is_success:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.PREDICT,
                                constant.ProgressStepStatusEnum.SUCCESS)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
            else:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.PREDICT,
                                constant.ProgressStepStatusEnum.ERROR)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))

    def load_model(self, model_path, config, gpu_device_index):
        model_type = model_path.split('.')[-1]
        if model_type in ('pth', 'bin', 'pt'):
            vocab = ['<pad>', '<mask>', 'A', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R',
                     'S', 'T', 'U', 'V', 'W', 'Y', 'C[57.02]', 'M[15.99]', 'N[.98]', 'Q[.98]', 'X', '<unk>']
            config["vocab"] = vocab
            # s2i = {v: i for i, v in enumerate(vocab)}
            # self.logger.info(f"Vocab: {s2i}, n_peaks: {config['n_peaks']}")
            if model_type == 'pth':
                model = torch.load(model_path)
            elif model_type == 'pt':
                model = DDA_BERT.load_pt(model_path, config)
                model.to(torch.bfloat16)
            else:
                model = torch.load_bin(model_path, config)
        elif model_type == 'ckpt':
            model, config = DDA_BERT.load_ckpt(model_path)
        else:
            raise Exception("Invalid model file type. Only .pth, .bin and .ckpt formats are supported.")

        model.eval()
        model = model.to(torch.bfloat16)
        model.to(f'cuda:{gpu_device_index}')
        return model

    def deal_one(self, config, mzml_info):
        self.deal_sage_pred(config, mzml_info)
        if not self.is_d:
            self.deal_alphapept_pred(config, mzml_info)
            self.deal_fragpipe_pred(config, mzml_info)

    def deal_sage_pred(self, config, mzml_info):
        self.logger.info('Process sage pred')
        ipc_path = os.path.join(self.sage_gen_pred_ipc_result_dir, '{}.ipc'.format(mzml_info.base_file_name))
        sage_pred_dir = os.path.join(self.pred_result_dir, 'sage')
        shutil.rmtree(sage_pred_dir, ignore_errors=True)
        os.makedirs(sage_pred_dir, exist_ok=True)

        self.evaluate(ipc_path, config, mzml_info.base_file_name, sage_pred_dir)
        self.logger.info('Finished sage pred')

    def deal_alphapept_pred(self, config, mzml_info):
        self.logger.info('Process alphapept pred')
        ipc_path = os.path.join(self.alphapept_gen_pred_output_dir, '{}.ipc'.format(mzml_info.base_file_name))

        alphapept_pred_dir = os.path.join(self.pred_result_dir, 'alphapept')
        shutil.rmtree(alphapept_pred_dir, ignore_errors=True)
        os.makedirs(alphapept_pred_dir, exist_ok=True)

        self.evaluate(ipc_path, config, mzml_info.base_file_name, alphapept_pred_dir)
        self.logger.info('Finished alphapept pred')

    def deal_fragpipe_pred(self, config, mzml_info):
        self.logger.info('Process fragpipe pred')
        ipc_path = os.path.join(self.fragpipe_fp_first_clean_dir, '{}.ipc'.format(mzml_info.base_file_name))
        fragpipe_pred_dir = os.path.join(self.pred_result_dir, 'fragpipe')
        shutil.rmtree(fragpipe_pred_dir, ignore_errors=True)
        os.makedirs(fragpipe_pred_dir, exist_ok=True)

        self.evaluate(ipc_path, config, mzml_info.base_file_name, fragpipe_pred_dir)
        self.logger.info('Finished fragpipe pred')

    '''
    load config
    '''

    def load_config(self):
        with open(self.base_config, 'r', encoding='utf-8') as f:
            config = yaml.load(f.read(), Loader=yaml.FullLoader)
        return config

    def load_pred_df(self, file_path):
        self.send_msg('Load ipc data, {}'.format(file_path))
        # load data
        df = pl.read_ipc(file_path)
        if 'modified_sequence' not in df.columns:
            df = df.with_columns(df['precursor_sequence'].alias('modified_sequence'))

        if 'index' not in df.columns:
            df = df.with_columns((pl.arange(0, pl.count()).alias("index")))
        if "weight" not in df.columns:
            df = df.with_columns(pl.lit(1.0).alias("weight"))
        return df

    '''
    gen dl
    '''

    def load_gen_dl(self, df, config, gpu_index):

        vocab = config["vocab"]
        s2i = {v: i for i, v in enumerate(vocab)}

        # add index when dataloader
        ds = SpectrumDataset(df, s2i, config["n_peaks"], need_label=True, need_weight=True, need_deltaRT=True,
                             need_index=True)

        # update batch_size with free_memory
        current_gpu_index = torch.device(f'cuda:{gpu_index}')
        total_memory = torch.cuda.get_device_properties(current_gpu_index).total_memory / (1024 ** 3)
        used_memory = torch.cuda.memory_allocated(current_gpu_index) / (1024 ** 3)
        free_memory = total_memory - used_memory
        predict_batch_size = int(free_memory * 100)
        self.logger.info(f'init predict_batch_size: {predict_batch_size}')

        predict_batch_size = min(predict_batch_size, config['predict_batch_size'])
        self.send_msg(
            "Eval, predict_bath_size: {}, config: {}".format(predict_batch_size, config['predict_batch_size']))

        dl = DataLoader(ds,
                        batch_size=config["predict_batch_size"],
                        num_workers=0,
                        shuffle=False,
                        collate_fn=collate_batch_weight_deltaRT_index,
                        pin_memory=True)

        self.send_msg('Data {} samples, DataLoader: {}'.format(len(ds), len(dl)))
        return dl

    def evaluate(self, ipc_path, config, base_file_name, result_dir):
        ipc_df = self.load_pred_df(ipc_path)
        process_num = len(self.input_param.gpu_devices)
        if process_num == 1:

            self.eval_one_thread(base_file_name, self.model_path, ipc_df, config, result_dir,
                                 self.input_param.gpu_devices[0])
        else:
            #
            chunk_df_list = self.split_df(ipc_df, process_num)
            processes = []
            for thread_index, gpu_device in enumerate(self.input_param.gpu_devices):
                process = Process(
                    target=self.eval_one_thread,
                    args=(base_file_name, self.model_path, deepcopy(chunk_df_list[thread_index]), config, result_dir,
                          gpu_device)
                )
                processes.append(process)
                process.start()

            for process in processes:
                process.join()

        all_df = []
        for root, dirs, files in os.walk(result_dir):
            for file_name in files:
                if file_name == f'{base_file_name}_pred.csv':
                    each_pred_df = pd.read_csv(os.path.join(root, file_name), header=None)
                    all_df.append(each_pred_df)
        save_df = pd.concat(all_df, axis=0)
        save_df.to_csv(os.path.join(result_dir, f'{base_file_name}_pred.csv'), index=False, header=False)

    def split_df(self, df, thread_num):
        num_rows = len(df)
        base_size = num_rows // thread_num
        remainder = num_rows % thread_num
        chunk_df_list = []
        start = 0
        for i in range(thread_num):
            size = base_size + (1 if i < remainder else 0)
            end = start + size
            #
            chunk = df[start:end]
            chunk_df_list.append(chunk)
            start = end
        return chunk_df_list

    def eval_one_thread(self, base_file_name, model_path, ipc_df, config, result_dir, gpu_index):
        model = self.load_model(model_path, config, gpu_index)
        thread_result_dir = os.path.join(result_dir, f'g{gpu_index}')
        os.makedirs(thread_result_dir, exist_ok=True)
        dl = self.load_gen_dl(ipc_df, config, gpu_index)
        self.send_msg('Processing pred')
        if self.env == constant.env_linux:
            strategy = DDPStrategy(gradient_as_bucket_view=True, find_unused_parameters=True, start_method='spawn')
        else:
            strategy = DDPStrategy(gradient_as_bucket_view=True, find_unused_parameters=True,
                                   process_group_backend="gloo", start_method='spawn')
        trainer = ptl.Trainer(
            accelerator="auto",
            devices=[int(gpu_index)],
            strategy=strategy,
        )
        evaluate = Evalute(thread_result_dir, base_file_name, model)
        trainer.test(evaluate, dataloaders=dl)
        return
