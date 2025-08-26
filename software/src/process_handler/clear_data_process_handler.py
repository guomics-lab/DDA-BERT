import shutil
import time

from src.common import constant
from src.common.runtime_data_info import runtime_data
from src.obj.Entity import MzMLFileInfo
from src.process_handler.common_process_handler import CommonProcessHandler


class ClearDataProcessHandler(CommonProcessHandler):

    def __init__(self, f_info: MzMLFileInfo, base_output_path, logger, env, input_param, start_time=0):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.CLEAR_DATA)

    def deal_process(self):
        self.send_msg('Processing clear data', status=constant.ProgressStepStatusEnum.RUNNING)
        try:
            if not runtime_data.current_is_success:
                self.send_msg(status=constant.ProgressStepStatusEnum.ERROR)
                return
            ss = time.time()
            if not self.input_param.open_clear_data:
                ee = time.time()
                self.send_msg('Finished clear data, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.SUCCESS)
                return

            shutil.rmtree(self.psm_output_path, ignore_errors=True)
            shutil.rmtree(self.sage_output_path, ignore_errors=True)
            shutil.rmtree(self.alphapept_result_output_path, ignore_errors=True)
            shutil.rmtree(self.combine_alphapept_sage_output_dir, ignore_errors=True)
            shutil.rmtree(self.fragpipe_output_path, ignore_errors=True)

            shutil.rmtree(self.pred_result_dir, ignore_errors=True)
            shutil.rmtree(self.protein_infer_dir, ignore_errors=True)

            ee = time.time()
            self.send_msg('Finished clear data, spend time: {}'.format(self.get_use_time(ss, ee)),
                          status=constant.ProgressStepStatusEnum.SUCCESS)
        except Exception as e:
            runtime_data.current_is_success = False
            self.logger.exception('Clear data exception')
            self.send_msg('Clear data exception: {}'.format(e), status=constant.ProgressStepStatusEnum.ERROR)
