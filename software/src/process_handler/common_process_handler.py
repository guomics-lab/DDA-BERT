import os
import time

from src import logger_utils
from src.obj.Entity import InputParam
from src.obj.Entity import MzMLFileInfo
from src.utils import file_utils


class CommonProcessHandler(object):

    def __init__(self, f_info: MzMLFileInfo, base_output_path, logger, env, input_param: InputParam, start_time=0,
                 step=None):
        self.env = env
        self.base_output_path = base_output_path
        self.f_info = f_info

        self.is_d = file_utils.is_d(self.f_info.file_path)

        self.step = step

        self.input_param = input_param

        os.makedirs(self.base_output_path, exist_ok=True)

        self.psm_output_path = os.path.join(self.base_output_path, 'psm')
        os.makedirs(self.psm_output_path, exist_ok=True)

        self.fragpipe_output_path = os.path.join(self.base_output_path, 'fragpipe')
        os.makedirs(self.fragpipe_output_path, exist_ok=True)

        self.fragpipe_fp_first_clean_dir = os.path.join(self.fragpipe_output_path, 'first_clean')
        os.makedirs(self.fragpipe_fp_first_clean_dir, exist_ok=True)

        self.sage_output_path = os.path.join(self.base_output_path, 'sage')
        os.makedirs(self.sage_output_path, exist_ok=True)

        self.sage_result_output_path = os.path.join(self.sage_output_path, 'result')
        os.makedirs(self.sage_result_output_path, exist_ok=True)

        self.sage_first_clean_dir = os.path.join(self.sage_output_path, 'sage_first_clean')
        os.makedirs(self.sage_first_clean_dir, exist_ok=True)

        self.pred_rt_v2_result_dir = os.path.join(self.sage_output_path, 'pred_rt_v2')
        os.makedirs(self.pred_rt_v2_result_dir, exist_ok=True)

        self.sage_gen_pred_ipc_result_dir = os.path.join(self.sage_output_path, 'gen_pred_ipc')
        os.makedirs(self.sage_gen_pred_ipc_result_dir, exist_ok=True)

        self.alphapept_result_output_path = os.path.join(self.base_output_path, 'alphapept')
        os.makedirs(self.alphapept_result_output_path, exist_ok=True)

        self.alphapept_data_clean_output_dir = os.path.join(self.alphapept_result_output_path, 'data_clean')
        os.makedirs(self.alphapept_data_clean_output_dir, exist_ok=True)

        self.alphapept_gen_pred_output_dir = os.path.join(self.alphapept_result_output_path, 'gen_pred')
        os.makedirs(self.alphapept_gen_pred_output_dir, exist_ok=True)

        self.combine_alphapept_sage_output_dir = os.path.join(self.base_output_path, 'combine_alphapept_sage')
        os.makedirs(self.combine_alphapept_sage_output_dir, exist_ok=True)

        self.pred_result_dir = os.path.join(self.base_output_path, 'pred')
        os.makedirs(self.pred_result_dir, exist_ok=True)

        self.protein_infer_dir = os.path.join(self.base_output_path, 'protein_infer')
        os.makedirs(self.protein_infer_dir, exist_ok=True)

        self.logger = self.get_logger(logger)
        self.start_time = start_time

        # self.stop_flag = False

    # def stop_process(self):
    #     self.stop_flag = True

    def deal_process(self):
        pass

    def send_msg(self, msg=None, status=None, with_time=True):
        if with_time and msg:
            msg = self.get_now_use_time() + ' ' + msg
        self.logger.info(msg)

    def get_now_use_time(self):
        now_time = time.time()
        minutes, seconds = divmod(now_time - self.start_time, 60)
        minutes = int(minutes)
        seconds = int(seconds)
        return '[{}:{}]'.format(minutes, str(seconds).zfill(2))

    def get_use_time(self, start_time, end_time):
        minutes, seconds = divmod(end_time - start_time, 60)
        minutes = int(minutes)
        seconds = int(seconds)
        return '[{}:{}]'.format(minutes, str(seconds).zfill(2))

    '''
    get logger
    '''

    def get_logger(self, logger):
        if not logger:
            log_file_path, logger = logger_utils.get_current_logger()
        return logger
