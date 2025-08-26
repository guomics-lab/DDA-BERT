import os.path
import time

from src import logger_utils
from src.common.runtime_data_info import runtime_data
from src.obj.Entity import InputParam
from src.obj.Entity import MzMLFileInfo
from src.process_handler.alphapept_process_handler import AlphapeptProcessHandler
from src.process_handler.clear_data_process_handler import ClearDataProcessHandler
from src.process_handler.fdr_calc_no_finetune_process_handler import FdrCalcNoFinetuneProcessHandler
from src.process_handler.fragpipe_process_handler import FragpipeProcessHandler
from src.process_handler.getdata_process_handler import GetDataProcessHandler
from src.process_handler.predict_process_handler import PredictProcessHandler
from src.process_handler.protein_infer_process_handler import ProteinInferProcessHandler
from src.process_handler.result_build_process_handler import ResultBuildProcessHandler
from src.process_handler.sage_process_handler import SageProcessHandler
from src.utils import file_utils


class MainProcessHandler(object):

    def __init__(self, input_param: InputParam):
        self.input_param = input_param
        self.current_logger = input_param.logger
        if not self.current_logger:
            log_file_path, self.current_logger = logger_utils.get_current_logger(self.input_param.output_dir_path)

        mzml_file_list = []
        for mzml_path in input_param.mzml_file_list:
            f_info = MzMLFileInfo()
            file_name = os.path.split(mzml_path)[-1]
            f_info.file_name = file_name
            f_info.base_file_name = file_name.removesuffix('.mzML').removesuffix('.d')
            f_info.file_path = mzml_path
            mzml_file_list.append(f_info)

        self.mzml_file_list = mzml_file_list

    def deal_one(self, current_mzml_index, each_mzml_info):
        runtime_data.current_mzml_index = current_mzml_index
        runtime_data.current_is_success = True

        is_d = file_utils.is_d(each_mzml_info.file_path)

        if self.input_param.env == 'win':
            if is_d:
                sage_exe_path = os.path.join(os.getcwd(), 'sage', 'windows_sage.exe')
                mzml_raw_spec_abs_path = os.path.join(os.getcwd(), 'sage', 'win_d_rawspectrum.exe')
            else:
                sage_exe_path = os.path.join(os.getcwd(), 'sage', 'windows_sage.exe')
                mzml_raw_spec_abs_path = os.path.join(os.getcwd(), 'sage', 'win_mzml_rawspecturm.exe')
        elif self.input_param.env == 'linux':
            if is_d:
                sage_exe_path = os.path.join(os.getcwd(), 'sage', 'linux_sage')
                mzml_raw_spec_abs_path = os.path.join(os.getcwd(), 'sage', 'linux_d_rawspectrum')
            else:
                sage_exe_path = os.path.join(os.getcwd(), 'sage', 'linux_sage')
                mzml_raw_spec_abs_path = os.path.join(os.getcwd(), 'sage', 'linux_mzml_rawspectrum')

        else:
            raise Exception('Env error')
        fragpipe_exe_path = os.path.join(os.getcwd(), 'FragPipe22_0', 'bin', 'fragpipe')

        base_alphapept_config_path = 'config/AlphaPept.yaml'

        base_output_path = os.path.join(self.input_param.output_dir_path, each_mzml_info.base_file_name)
        start_time = time.time()

        self.sage_process = SageProcessHandler(sage_exe_path, mzml_raw_spec_abs_path, each_mzml_info,
                                               self.input_param.fasta_path,
                                               base_output_path,
                                               self.current_logger, self.input_param.env, self.input_param,
                                               start_time=start_time)

        self.alphapept_process_handler = AlphapeptProcessHandler(base_alphapept_config_path, each_mzml_info,
                                                                 self.input_param.fasta_path, base_output_path,
                                                                 self.current_logger, self.input_param.env,
                                                                 self.input_param,
                                                                 start_time=start_time)

        self.fragpipe_process_handler = FragpipeProcessHandler(fragpipe_exe_path, each_mzml_info,
                                                               self.input_param.fasta_path, base_output_path,
                                                               self.current_logger, self.input_param.env,
                                                               self.input_param,
                                                               start_time=start_time)

        self.get_data_process = GetDataProcessHandler(each_mzml_info, base_output_path, self.current_logger,
                                                      self.input_param.env,
                                                      self.input_param,
                                                      start_time=start_time)

        self.pred_process = PredictProcessHandler(each_mzml_info, base_output_path, self.input_param.model_path,
                                                  self.current_logger,
                                                  self.input_param.env,
                                                  self.input_param,
                                                  start_time=start_time)

        self.fdr_calc_no_finetune_process = FdrCalcNoFinetuneProcessHandler(each_mzml_info, base_output_path,
                                                                            self.current_logger,
                                                                            self.input_param.env, self.input_param,
                                                                            start_time=start_time)
        #
        self.protein_infer_process = ProteinInferProcessHandler(each_mzml_info, base_output_path, self.current_logger,
                                                                self.input_param.env, self.input_param,
                                                                start_time=start_time)
        #
        self.result_build_process = ResultBuildProcessHandler(each_mzml_info, base_output_path, self.current_logger,
                                                              self.input_param.env, self.input_param,
                                                              start_time=start_time)

        self.fragpipe_process_handler.deal_process()
        self.sage_process.deal_process()
        self.alphapept_process_handler.deal_process()

        self.get_data_process.deal_process()
        self.pred_process.deal_process()

        self.fdr_calc_no_finetune_process.deal_process()
        self.protein_infer_process.deal_process()
        #
        self.result_build_process.deal_process()

        # # clear data
        self.clear_data_process = ClearDataProcessHandler(each_mzml_info, base_output_path, self.current_logger,
                                                          self.input_param.env, self.input_param,
                                                          start_time=start_time)
        self.clear_data_process.deal_process()

    def deal_processing(self):
        ts = time.time()
        for current_mzml_index, each_mzml_info in enumerate(self.mzml_file_list):
            try:
                self.deal_one(current_mzml_index, each_mzml_info)
            except Exception:
                self.current_logger.exception('Process file {} exception'.format(each_mzml_info.file_name))
                runtime_data.current_is_success = False
            # important, do not delete
            print('ProgressNumMsg:{}'.format(current_mzml_index + 1))
        te = time.time()
        self.current_logger.info('Done, all time {}s'.format(te - ts))
        return 1

    def processing(self):
        try:
            self.deal_processing()
        except Exception:
            self.current_logger.exception('Process exception')
