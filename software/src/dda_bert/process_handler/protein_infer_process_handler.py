import json
import os
import time

from ..common import constant
from ..common.runtime_data_info import runtime_data
from ..obj.Entity import MzMLFileInfo, ProcessMsg
from ..process_handler.common_process_handler import CommonProcessHandler
from ..services import alphapept_fasta_digest_pep
from ..services import alphapept_protien_infer_v3 as protein_info_service


class ProteinInferProcessHandler(CommonProcessHandler):

    def __init__(self, f_info: MzMLFileInfo, base_output_path, logger, env, input_param, start_time=0, open_sage=False,
                 open_fp=False, open_ap=False):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.PROTEIN_INFER)

        self.open_sage = open_sage
        self.open_fp = open_fp
        self.open_ap = open_ap

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

        alphepept_pept_dict_path = None
        if self.open_ap:
            alphepept_pept_dict_path = alphapept_fasta_digest_pep.deal(
                os.path.join(self.alphapept_result_output_path, f'{self.f_info.base_file_name}_alphapept_config.yaml'),
                self.logger)

        protein_info_service.merge(mzml_info.base_file_name,
                                   os.path.join(self.pred_result_dir,
                                                'protein_infer_input.csv'),
                                   os.path.join(self.sage_result_output_path,
                                                f'{self.f_info.base_file_name}.sage.tsv'),
                                   os.path.join(self.fragpipe_fp_first_clean_dir, f'{self.f_info.base_file_name}.csv'),
                                   alphepept_pept_dict_path,
                                   self.protein_infer_dir, self.open_sage, self.open_fp, self.open_ap)
