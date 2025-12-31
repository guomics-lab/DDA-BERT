import json
import os
import time

import pandas as pd

from ..common import constant
from ..common.runtime_data_info import runtime_data
from ..obj.Entity import MzMLFileInfo, ProcessMsg
from ..process_handler.common_process_handler import CommonProcessHandler


class ResultBuildProcessHandler(CommonProcessHandler):

    def __init__(self, f_info: MzMLFileInfo, base_output_path, logger, env, input_param, start_time=0):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.RESULT_BUILD)

    def deal_process(self):
        self.send_msg('Processing result build', status=constant.ProgressStepStatusEnum.RUNNING)
        pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.RESULT_BUILD,
                        constant.ProgressStepStatusEnum.RUNNING)
        print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
        try:
            if not runtime_data.current_is_success:
                self.send_msg(status=constant.ProgressStepStatusEnum.ERROR)
                return

            ss = time.time()
            self.send_msg('Processing build result file')

            eval_fdr_peptide_file = os.path.join(self.pred_result_dir, 'fdr_peptide.csv')
            eval_fdr_psms_file = os.path.join(self.pred_result_dir, 'fdr_psms.csv')

            target_fdr_peptide_file = os.path.join(self.base_output_path, 'fdr_peptide.csv')

            # stander column name
            fdr_peptide_df = pd.read_csv(eval_fdr_peptide_file)

            fdr_peptide_df.rename(columns={'cleaned_sequence': 'PeptideSequence',
                                           'precursor_mz': 'PrecursorMz',
                                           'precursor_charge': 'PrecursorCharge',
                                           'modified_sequence': 'ModifiedSequence',
                                           'label': 'Label',
                                           'score': 'Score',
                                           'q_value': 'Qvalue',
                                           'scan_number': 'ScanNumber'}, inplace=True)

            fdr_peptide_df.to_csv(target_fdr_peptide_file, index=False)

            peptide_protein_path = os.path.join(self.protein_infer_dir,
                                                '{}_protein.csv'.format(self.f_info.base_file_name))
            peptide_protein_df = pd.read_csv(peptide_protein_path)

            # out protein
            protein_result_file = os.path.join(self.base_output_path,
                                               '{}_protein_group.tsv'.format(self.f_info.base_file_name))
            peptide_protein_df['FileName'] = self.f_info.file_name
            peptide_protein_df.rename(columns={'precursor_id': 'PrecursorID',
                                               'precursor_charge': 'PrecursorCharge',
                                               'modified_sequence': 'ModifiedSequence',
                                               'sequence': 'PeptideSequence',
                                               'protein': 'ProteinGroup',
                                               'q_value': 'ProteinQvalue'}, inplace=True)
            protein_save_col = ['FileName', 'PrecursorID', 'PrecursorCharge', 'ModifiedSequence', \
                                'PeptideSequence', 'ProteinGroup', 'ProteinQvalue']

            precursor_id_to_protein_group = \
            peptide_protein_df[['PrecursorID', 'ProteinGroup']].set_index('PrecursorID')[
                'ProteinGroup'].to_dict()

            peptide_protein_df.drop_duplicates(subset=['ProteinGroup'], inplace=True)
            peptide_protein_df[protein_save_col].to_csv(protein_result_file, sep='\t', index=False)

            precursor_id_to_sequence = peptide_protein_df[['PrecursorID', 'PeptideSequence']].set_index('PrecursorID')[
                'PeptideSequence'].to_dict()


            result_file = os.path.join(self.base_output_path, '{}_psm.tsv'.format(self.f_info.base_file_name))
            fdr_df = pd.read_csv(eval_fdr_psms_file)
            fdr_df['PrecursorID'] = fdr_df['precursor_charge'].astype(str) + '_' + fdr_df['modified_sequence'].astype(str)
            fdr_df['PeptideSequence'] = fdr_df['PrecursorID'].apply(lambda x: precursor_id_to_sequence.get(x))
            fdr_df['ProteinGroup'] = fdr_df['PrecursorID'].apply(lambda x: precursor_id_to_protein_group.get(x))
            fdr_df['FileName'] = self.f_info.file_name
            #
            fdr_df.rename(columns={'scan_number': 'ScanNumber',
                                   'precursor_mz': 'PrecursorMz',
                                   'modified_sequence': 'ModifiedSequence',
                                   'precursor_charge': 'PrecursorCharge',
                                   'q_value': 'Qvalue'}, inplace=True)
            #
            psm_save_col = ['FileName', 'ScanNumber', 'PrecursorID', 'PrecursorMz', 'PeptideSequence',
                            'ModifiedSequence', \
                            'PrecursorCharge', 'ProteinGroup', 'Qvalue']
            fdr_df = fdr_df[psm_save_col]
            fdr_df.sort_values(by=['ScanNumber', 'PrecursorID'], inplace=True)
            fdr_df.reset_index(drop=True, inplace=True)
            fdr_df['PSMID'] = fdr_df.index + 1
            fdr_df.to_csv(result_file, sep='\t', index=False)

            ee = time.time()
            self.send_msg('Finished build result file, spend time: {}'.format(self.get_use_time(ss, ee)),
                          status=constant.ProgressStepStatusEnum.SUCCESS)
        except Exception as e:
            runtime_data.current_is_success = False
            self.logger.exception('Build result file exception')
            self.send_msg('Build result file exception: {}'.format(e), status=constant.ProgressStepStatusEnum.ERROR)
        finally:
            if runtime_data.current_is_success:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.RESULT_BUILD,
                                constant.ProgressStepStatusEnum.SUCCESS)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
            else:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.RESULT_BUILD,
                                constant.ProgressStepStatusEnum.ERROR)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
