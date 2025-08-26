import json
import os
import time

from src.common import constant
from src.common.runtime_data_info import runtime_data
from src.obj.Entity import MzMLFileInfo, ProcessMsg
from src.process_handler.common_process_handler import CommonProcessHandler
from src.services import alphapept_fasta_digest_pep
from src.services import alphapept_protien_infer_v2 as protein_info_service


class ProteinInferProcessHandler(CommonProcessHandler):

    def __init__(self, f_info: MzMLFileInfo, base_output_path, logger, env, input_param, start_time=0):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.PROTEIN_INFER)

    def deal_process(self):
        self.send_msg('Processing protein infer', status=constant.ProgressStepStatusEnum.RUNNING)
        pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.PROTEIN_INFER,
                        constant.ProgressStepStatusEnum.RUNNING)
        print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
        try:
            if not runtime_data.current_is_success:
                self.send_msg(status=constant.ProgressStepStatusEnum.ERROR)
                return
            ss = time.time()
            if not self.input_param.open_protein_infer:
                ee = time.time()
                self.send_msg('Finished protein infer, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.SUCCESS)
                return
            mzml_info = self.f_info
            self.send_msg('Processing protein infer file {}'.format(mzml_info.base_file_name))
            self.deal_one(mzml_info)
            ee = time.time()
            self.send_msg('Finished protein infer, spend time: {}'.format(self.get_use_time(ss, ee)),
                          status=constant.ProgressStepStatusEnum.SUCCESS)
        except Exception as e:
            runtime_data.current_is_success = False
            self.logger.exception('Protein infer exception')
            self.send_msg('Protein infer exception: {}'.format(e), status=constant.ProgressStepStatusEnum.ERROR)
        finally:
            if runtime_data.current_is_success:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.PROTEIN_INFER,
                                constant.ProgressStepStatusEnum.SUCCESS)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
            else:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.PROTEIN_INFER,
                                constant.ProgressStepStatusEnum.ERROR)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))

    def deal_one(self, mzml_info):

        alphepept_pept_dict_path = alphapept_fasta_digest_pep.deal(
            os.path.join(self.alphapept_result_output_path, f'{self.f_info.base_file_name}_alphapept_config.yaml'),
            os.path.dirname(self.input_param.fasta_path), self.logger)

        protein_info_service.merge(mzml_info.base_file_name,
                                   os.path.join(self.pred_result_dir,
                                                'protein_infer_input.csv'),
                                   os.path.join(self.sage_result_output_path,
                                                f'{self.f_info.base_file_name}.sage.tsv'),
                                   os.path.join(self.fragpipe_fp_first_clean_dir, f'{self.f_info.base_file_name}.csv'),
                                   alphepept_pept_dict_path,
                                   self.protein_infer_dir)
