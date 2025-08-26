from __future__ import annotations

import json
import os
import time

import numpy as np

from src.common import constant
from src.common.runtime_data_info import runtime_data
from src.obj.Entity import MzMLFileInfo, ProcessMsg
from src.process_handler.common_process_handler import CommonProcessHandler
from src.services import data_combine


class FdrCalcNoFinetuneProcessHandler(CommonProcessHandler):

    def __init__(self, f_info: MzMLFileInfo, base_output_path, logger, env, input_param, start_time=0):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.FDR_CALC)

    def deal_process(self):
        self.send_msg('Processing calc fdr no finetune', status=constant.ProgressStepStatusEnum.RUNNING)
        pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.FDR_CALC,
                        constant.ProgressStepStatusEnum.RUNNING)
        print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))

        try:
            ss = time.time()
            if not runtime_data.current_is_success:
                self.send_msg(status=constant.ProgressStepStatusEnum.ERROR)
                return
            if not self.input_param.open_fdr_calc:
                ee = time.time()
                self.send_msg('Finished calc fdr, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.SUCCESS)
                return
            ss = time.time()
            mzml_info = self.f_info
            self.send_msg('Processing calc file {}'.format(mzml_info.base_file_name))

            sage_ipc = os.path.join(self.sage_gen_pred_ipc_result_dir, '{}.ipc'.format(mzml_info.base_file_name))
            ap_ipc = os.path.join(self.alphapept_gen_pred_output_dir, '{}.ipc'.format(mzml_info.base_file_name))
            fp_ipc = os.path.join(self.fragpipe_fp_first_clean_dir, '{}.ipc'.format(mzml_info.base_file_name))
            sage_pred = os.path.join(self.pred_result_dir, 'sage', f'{mzml_info.base_file_name}_pred.csv')
            ap_pred = os.path.join(self.pred_result_dir, 'alphapept', f'{mzml_info.base_file_name}_pred.csv')
            fp_pred = os.path.join(self.pred_result_dir, 'fragpipe', f'{mzml_info.base_file_name}_pred.csv')
            output_dir = self.pred_result_dir

            data_combine.deal_process(self.f_info.base_file_name, sage_ipc, ap_ipc, fp_ipc, sage_pred, ap_pred, fp_pred,
                                      output_dir)

            ee = time.time()
            self.send_msg('Finished calc fdr, spend time: {}'.format(self.get_use_time(ss, ee)),
                          status=constant.ProgressStepStatusEnum.SUCCESS)
        except Exception as e:
            runtime_data.current_is_success = False
            self.logger.exception('Calc fdr exception')
            self.send_msg('Calc fdr exception: {}'.format(e), status=constant.ProgressStepStatusEnum.ERROR)
        finally:
            if runtime_data.current_is_success:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.FDR_CALC,
                                constant.ProgressStepStatusEnum.SUCCESS)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
            else:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.FDR_CALC,
                                constant.ProgressStepStatusEnum.ERROR)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))

    def get_fdr_result(self, df):
        df['score'] = df['score'].astype(float)
        df['label'] = df['label'].astype(int)
        df = df.sort_values(by='score', ascending=False, ignore_index=True)
        df['decoy'] = np.where(df['label'] == 1, 0, 1)
        target_num = (df.decoy == 0).cumsum()
        decoy_num = (df.decoy == 1).cumsum()
        target_num[target_num == 0] = 1
        decoy_num[decoy_num == 0] = 1
        df['q_value'] = decoy_num / target_num
        df['q_value'] = df['q_value'][::-1].cummin()
        return df

    def calc_fdr(self, df):
        print(f'org calc_fdr df len: {df.shape}')
        df['precursor_id'] = df['precursor_charge'].astype(str) + '_' + df['modified_sequence']
        df = df.sort_values(by='score', ascending=False, ignore_index=True)
        df = df.drop_duplicates(subset=['scan_number'], keep='first')
        df['pred_label'] = df['score'].apply(lambda x: 1 if x > 0.5 else 0)
        print(f'calc_fdr df len: {df.shape}')

        # finetune_eval_df_temp_dir = os.path.join(self.finetune_eval_dir, 'temp')
        pred_eval_df_temp_dir = os.path.join(self.pred_result_dir, 'temp')
        os.makedirs(pred_eval_df_temp_dir, exist_ok=True)
        df.to_csv(os.path.join(pred_eval_df_temp_dir, f'df_all.csv'), index=0)
        df = self.get_fdr_result(df)

        stat = {}
        for threshold in [0.001 * i for i in range(1, 11)]:
            stat[f'{threshold:.3f}_fdr'] = len(df[(df['q_value'] <= threshold) & (df['label'] == 1)])

        #
        fdr_df = df[(df['q_value'] <= 0.01) & (df['label'] == 1)]
        if 'protein_new' in df.columns:
            fdr_df = fdr_df[['precursor_mz', 'precursor_charge', 'modified_sequence', 'protein_new', 'cleaned_sequence',
                             'label', 'score', 'q_value', 'scan_number', 'precursor_id']]
        else:
            fdr_df = fdr_df[['precursor_mz', 'precursor_charge', 'modified_sequence',
                             'label', 'score', 'q_value', 'scan_number', 'precursor_id']]
            fdr_df['cleaned_sequence'] = fdr_df['modified_sequence']
            fdr_df['cleaned_sequence'] = fdr_df['cleaned_sequence'].str.replace('n\[42\]', '', regex=True)
            fdr_df['cleaned_sequence'] = fdr_df['cleaned_sequence'].str.replace('N\[.98\]', 'N', regex=True)
            fdr_df['cleaned_sequence'] = fdr_df['cleaned_sequence'].str.replace('Q\[.98\]', 'Q', regex=True)
            fdr_df['cleaned_sequence'] = fdr_df['cleaned_sequence'].str.replace('M\[15.99\]', 'M', regex=True)
            fdr_df['cleaned_sequence'] = fdr_df['cleaned_sequence'].str.replace('C\[57.02\]', 'C', regex=True)

        fdr_df.to_csv(os.path.join(self.pred_result_dir, f'fdr_psms.csv'), index=0)

        # fdr_precursor
        fdr_precursor = fdr_df.sort_values(by='score', ascending=False, ignore_index=True)
        fdr_precursor = fdr_precursor.drop_duplicates(subset=['precursor_id'], keep='first')
        fdr_precursor.to_csv(os.path.join(self.pred_result_dir, f'fdr_precursors.csv'),
                             index=False)

        fdr_peptide = fdr_precursor.drop_duplicates(subset=['cleaned_sequence'], keep='first')
        fdr_peptide = fdr_peptide[['cleaned_sequence', 'precursor_mz', 'precursor_charge', 'modified_sequence',
                                   'label', 'score', 'q_value', 'scan_number', 'cleaned_sequence']]
        fdr_peptide.to_csv(os.path.join(self.pred_result_dir, f'fdr_peptide.csv'), index=False)
        self.send_msg('Saved fdr result to {}'.format(self.pred_result_dir))

        return stat, fdr_df, fdr_precursor, fdr_peptide
