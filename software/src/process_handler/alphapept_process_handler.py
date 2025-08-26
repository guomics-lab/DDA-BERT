import json
import os.path
import shutil
import time

import alphapept.interface
import alphapept.interface
import yaml
from alphapept.settings import load_settings

from src.common import constant
from src.common.runtime_data_info import runtime_data
from src.obj.Entity import MzMLFileInfo, ProcessMsg
from src.process_handler.common_process_handler import CommonProcessHandler


class AlphapeptProcessHandler(CommonProcessHandler):

    def __init__(self, base_alphapept_config_path, f_info: MzMLFileInfo, fasta_path, base_output_path, logger,
                 env, input_param,
                 start_time=0):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.ALPHAPEPT)
        self.base_alphapept_config_path = base_alphapept_config_path
        self.fasta_path = fasta_path

    def deal_process(self):
        pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.ALPHAPEPT,
                        constant.ProgressStepStatusEnum.RUNNING)
        print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
        self.send_msg('Processing alphapept', status=constant.ProgressStepStatusEnum.RUNNING)
        try:
            ss = time.time()
            if not self.input_param.open_alphapept:
                ee = time.time()
                self.send_msg(msg='Finished alphapept, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.SUCCESS)
                return

            if self.is_d:
                #
                ee = time.time()
                self.send_msg(msg='D file skip alphapept, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.SUCCESS)
                return

            ms_data_path = os.path.join(self.alphapept_result_output_path, f'{self.f_info.base_file_name}.ms_data.hdf')

            run_alphapept_config_path = self.build_config()
            settings = load_settings(run_alphapept_config_path)

            alphapept.interface.run_complete_workflow(settings)

            ee = time.time()
            if not os.path.exists(ms_data_path):
                self.send_msg(msg='Finished alphapept, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.ERROR)
                runtime_data.current_is_success = False
                return

            self.send_msg(msg='Finished alphapept, spend time: {}'.format(self.get_use_time(ss, ee)),
                          status=constant.ProgressStepStatusEnum.SUCCESS)

        except Exception:
            runtime_data.current_is_success = False
            self.logger.exception('Processing alphapept exception')
            self.send_msg(msg='Alphapept exception', status=constant.ProgressStepStatusEnum.ERROR)
        finally:
            if runtime_data.current_is_success:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.ALPHAPEPT,
                                constant.ProgressStepStatusEnum.SUCCESS)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
            else:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.ALPHAPEPT,
                                constant.ProgressStepStatusEnum.ERROR)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))

    '''
    
    '''

    def build_config(self):
        current_mzml_base_name = self.f_info.base_file_name
        run_alphapept_config_path = os.path.join(self.alphapept_result_output_path,
                                                 '{}_alphapept_config.yaml'.format(current_mzml_base_name))
        self.send_msg('Processing build alphapept config, path is {}'.format(run_alphapept_config_path))

        #
        dest_temp_file_path = os.path.join(self.alphapept_result_output_path,
                                           str(self.f_info.file_name).replace('.mzML', '.raw'))
        source_raw_file_path = str(self.f_info.file_path).replace('.mzML', '.raw')
        shutil.copy(source_raw_file_path, dest_temp_file_path)

        results_path = os.path.join(self.alphapept_result_output_path, f'{self.f_info.base_file_name}_results.hdf')

        #
        run_config = self.get_config_dict(dest_temp_file_path, results_path)
        yaml.dump(run_config, open(run_alphapept_config_path, 'w'))
        self.send_msg('Finished build alphapept config')
        return run_alphapept_config_path

    '''
    
    '''

    def get_config_dict(self, dest_temp_file_path, results_path):

        base_config = load_settings(self.base_alphapept_config_path)
        base_config['experiment']['fasta_paths'] = [self.fasta_path]
        base_config['experiment']['file_paths'] = [dest_temp_file_path]
        base_config['experiment']['results_path'] = results_path
        base_config['experiment']['shortnames'] = f'{self.f_info.base_file_name}'

        mods_fixed = []
        if self.input_param.static_mods_c:
            mods_fixed.append('cC')

        variable_mods = []
        if self.input_param.variable_mods_m:
            variable_mods.append('oxM')
        if self.input_param.variable_mods_n:
            variable_mods.append('deamN')
        if self.input_param.variable_mods_q:
            variable_mods.append('deamQ')

        mods_variable_terminal_prot = []
        if self.input_param.variable_mods_a:
            mods_variable_terminal_prot.append('a<^')

        base_config['fasta']['n_missed_cleavages'] = self.input_param.enzyme_missed_cleavages
        base_config['fasta']['pep_length_min'] = self.input_param.digest_min_len
        base_config['fasta']['pep_length_max'] = self.input_param.digest_max_len
        base_config['fasta']['mods_fixed'] = mods_fixed
        base_config['fasta']['mods_variable'] = variable_mods
        base_config['fasta']['mods_variable_terminal_prot'] = mods_variable_terminal_prot

        base_config['search']['top_n'] = self.input_param.report_psms
        base_config['search']['prec_tol'] = self.input_param.max_precursor_tol_ppm
        base_config['search']['frag_tol'] = self.input_param.max_fragment_tol_ppm

        base_config['raw']['n_most_abundant'] = self.input_param.max_peaks

        base_config['features']['iso_charge_min'] = self.input_param.min_precursor_charge
        base_config['features']['iso_charge_max'] = self.input_param.max_precursor_charge

        return base_config
