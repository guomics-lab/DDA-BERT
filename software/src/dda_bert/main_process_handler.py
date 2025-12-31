import os.path
import shutil
import time

from . import logger_utils
from .common.runtime_data_info import runtime_data
from .obj.Entity import InputParam
from .obj.Entity import MzMLFileInfo
from .process_handler.alphapept_process_handler import AlphapeptProcessHandler
from .process_handler.clear_data_process_handler import ClearDataProcessHandler
from .process_handler.fdr_calc_no_finetune_process_handler import FdrCalcNoFinetuneProcessHandler
from .process_handler.fragpipe_process_handler import FragpipeProcessHandler
from .process_handler.getdata_process_handler import GetDataProcessHandler
from .process_handler.predict_process_handler import PredictProcessHandler
from .process_handler.protein_infer_process_handler import ProteinInferProcessHandler
from .process_handler.result_build_process_handler import ResultBuildProcessHandler
from .process_handler.sage_process_handler import SageProcessHandler


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
        if input_param.only_run_predict:
            run_identify_type_arr = str(input_param.search_engines).split(',')
            run_identify_type_arr = list(set(run_identify_type_arr))
            if len(run_identify_type_arr) == 2 and 'ap' in run_identify_type_arr:
                raise Exception('Identify type error')
            if len(run_identify_type_arr) == 1 and run_identify_type_arr[0] == 'ap':
                raise Exception('Identify type error, not support only ap')

            self.open_sage = 'sage' in run_identify_type_arr
            self.open_fp = 'fp' in run_identify_type_arr
            self.open_ap = 'ap' in run_identify_type_arr
            # 加个强制判断，如果是raw文件可以3个都开启，如果是非raw，则只能开启sage，并且ap不能单独开启
            if input_param.mass_file_format == 'd' or input_param.mass_file_format == 'wiff':
                self.open_fp = False
                self.open_ap = False
        else:
            if input_param.mass_file_format == 'raw':
                self.open_sage = True
                self.open_fp = True
                self.open_ap = True
            else:
                self.open_sage = True
                self.open_fp = False
                self.open_ap = False

    def _load_alphapept_config(self):
        pass

    def _deal_one(self, current_mzml_index, each_mzml_info):
        runtime_data.current_mzml_index = current_mzml_index
        runtime_data.current_is_success = True

        # is_d = file_utils.is_d(each_mzml_info.file_path)
        is_d = self.input_param.mass_file_format == 'd'

        if is_d:
            sage_exe_path = os.path.join(os.getcwd(), 'sage', 'linux_sage')
            mzml_raw_spec_abs_path = os.path.join(os.getcwd(), 'sage', 'linux_d_rawspectrum')
        else:
            sage_exe_path = os.path.join(os.getcwd(), 'sage', 'linux_sage')
            mzml_raw_spec_abs_path = os.path.join(os.getcwd(), 'sage', 'linux_mzml_rawspectrum')

        fragpipe_exe_path = os.path.join(os.getcwd(), 'FragPipe22_0', 'bin', 'fragpipe')

        base_output_path = os.path.join(self.input_param.output_dir_path, each_mzml_info.base_file_name)
        start_time = time.time()

        self.sage_process = SageProcessHandler(sage_exe_path, mzml_raw_spec_abs_path, each_mzml_info,
                                               self.input_param.fasta_path,
                                               self.input_param.sage_config_file,
                                               base_output_path,
                                               self.current_logger, self.input_param.env, self.input_param,
                                               start_time=start_time)

        self.alphapept_process_handler = AlphapeptProcessHandler(each_mzml_info,
                                                                 self.input_param.fasta_path,
                                                                 self.input_param.ap_config_file,
                                                                 base_output_path,
                                                                 self.current_logger, self.input_param.env,
                                                                 self.input_param,
                                                                 start_time=start_time)

        self.fragpipe_process_handler = FragpipeProcessHandler(fragpipe_exe_path, each_mzml_info,
                                                               self.input_param.fasta_path,
                                                               self.input_param.fp_workflow_file,
                                                               base_output_path,
                                                               self.current_logger, self.input_param.env,
                                                               self.input_param,
                                                               start_time=start_time)

        self.get_data_process = GetDataProcessHandler(each_mzml_info, mzml_raw_spec_abs_path, base_output_path,
                                                      self.current_logger,
                                                      self.input_param.env,
                                                      self.input_param,
                                                      start_time=start_time, open_sage=self.open_sage,
                                                      open_fp=self.open_fp, open_ap=self.open_ap)

        self.pred_process = PredictProcessHandler(each_mzml_info, base_output_path, self.input_param.model_path,
                                                  self.current_logger,
                                                  self.input_param.env,
                                                  self.input_param,
                                                  start_time=start_time, open_sage=self.open_sage,
                                                  open_fp=self.open_fp, open_ap=self.open_ap)

        self.fdr_calc_no_finetune_process = FdrCalcNoFinetuneProcessHandler(each_mzml_info, base_output_path,
                                                                            self.current_logger,
                                                                            self.input_param.env, self.input_param,
                                                                            start_time=start_time,
                                                                            open_sage=self.open_sage,
                                                                            open_fp=self.open_fp, open_ap=self.open_ap)
        #
        self.protein_infer_process = ProteinInferProcessHandler(each_mzml_info, base_output_path, self.current_logger,
                                                                self.input_param.env, self.input_param,
                                                                start_time=start_time, open_sage=self.open_sage,
                                                                open_fp=self.open_fp, open_ap=self.open_ap)
        #
        self.result_build_process = ResultBuildProcessHandler(each_mzml_info, base_output_path, self.current_logger,
                                                              self.input_param.env, self.input_param,
                                                              start_time=start_time)
        if not self.input_param.only_run_predict:
            self.sage_process.deal_process()
            self.alphapept_process_handler.deal_process()
            self.fragpipe_process_handler.deal_process()

        else:
            # copy input file to output folder
            if self.input_param.mzml_ipc_file_dir:
                source_psm_ipc_file_path = os.path.join(self.input_param.mzml_ipc_file_dir,
                                                        f'{each_mzml_info.base_file_name}.ipc')
                target_psm_ipc_file_path = os.path.join(self.sage_process.psm_output_path, f'{each_mzml_info.base_file_name}.ipc')
                shutil.copy(source_psm_ipc_file_path, target_psm_ipc_file_path)

            if self.open_sage:
                sage_tsv = os.path.join(self.input_param.sage_file_dir, f'{each_mzml_info.base_file_name}.sage.tsv')
                target_sage_tsv = os.path.join(self.sage_process.sage_result_output_path, 'results.sage.tsv')

                shutil.copy(sage_tsv, target_sage_tsv)
            if self.open_fp:
                source_fp_csv = os.path.join(self.input_param.fp_file_dir, f'{each_mzml_info.base_file_name}.pin')

                target_pin_file_path = os.path.join(self.pred_process.fragpipe_output_path, 'exp',
                                                    f'{each_mzml_info.base_file_name}.pin')
                target_pin_file_dir = os.path.dirname(target_pin_file_path)
                os.makedirs(target_pin_file_dir, exist_ok=True)

                shutil.copy(source_fp_csv, target_pin_file_path)
            if self.open_ap:
                source_ap_pred_ipc = os.path.join(self.input_param.ap_file_dir,
                                                  f'{each_mzml_info.base_file_name}.ms_data.hdf')
                target_hdf_file_path = os.path.join(self.pred_process.alphapept_result_output_path,
                                                    f'{each_mzml_info.base_file_name}.ms_data.hdf')

                source_config_file = os.path.join(self.input_param.ap_file_dir,
                                                  f'{each_mzml_info.base_file_name}_results.yaml')
                target_config_file = os.path.join(self.pred_process.alphapept_result_output_path,
                                                  f'{each_mzml_info.base_file_name}_alphapept_config.yaml')

                shutil.copy(source_ap_pred_ipc, target_hdf_file_path)
                shutil.copy(source_config_file, target_config_file)

        self.get_data_process.deal_process()
        self.pred_process.deal_process()

        self.fdr_calc_no_finetune_process.deal_process()
        self.protein_infer_process.deal_process()
        #
        self.result_build_process.deal_process()
        #
        # # clear data
        self.clear_data_process = ClearDataProcessHandler(each_mzml_info, base_output_path, self.current_logger,
                                                          self.input_param.env, self.input_param,
                                                          start_time=start_time)
        self.clear_data_process.deal_process()

    def _deal_processing(self):
        ts = time.time()
        for current_mzml_index, each_mzml_info in enumerate(self.mzml_file_list):
            try:
                self._deal_one(current_mzml_index, each_mzml_info)
            except Exception:
                self.current_logger.exception('Process file {} exception'.format(each_mzml_info.file_name))
                runtime_data.current_is_success = False
            # important, do not delete
            print('ProgressNumMsg:{}'.format(current_mzml_index + 1))
        te = time.time()
        self.current_logger.info('Done, all time {}s'.format(te - ts))
        return 1

    def _processing(self):
        try:
            self._deal_processing()
        except Exception:
            self.current_logger.exception('Process exception')
