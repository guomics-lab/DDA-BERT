import json
import os.path
import shutil
import time

import alphapept.interface
import alphapept.interface
import yaml
from alphapept.settings import load_settings

from ..common import constant
from ..common.runtime_data_info import runtime_data
from ..common_config import get_alpha_config_path
from ..obj.Entity import MzMLFileInfo, ProcessMsg
from ..process_handler.common_process_handler import CommonProcessHandler


class AlphapeptProcessHandler(CommonProcessHandler):

    def __init__(self, f_info: MzMLFileInfo, fasta_path, ap_config_file, base_output_path, logger,
                 env, input_param,
                 start_time=0):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.ALPHAPEPT)
        # self.base_alphapept_config_path = base_alphapept_config_path
        self.fasta_path = fasta_path
        self.ap_config_file = ap_config_file

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

            if not self.is_raw:
                self.send_msg(msg='Not raw file, skip alphapept')
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

        if not self.ap_config_file:
            self.send_msg('Custom ap config not set, use default config')
            base_config = load_settings(get_alpha_config_path())
        else:
            base_config = load_settings(self.ap_config_file)

        base_config['experiment']['fasta_paths'] = [self.fasta_path]
        base_config['experiment']['file_paths'] = [dest_temp_file_path]
        base_config['experiment']['results_path'] = results_path
        base_config['experiment']['shortnames'] = f'{self.f_info.base_file_name}'

        if not self.ap_config_file:
            mods_fixed = []
            mods_fixed.append('cC')

            variable_mods = []
            variable_mods.append('oxM')
            variable_mods.append('deamN')
            variable_mods.append('deamQ')

            mods_variable_terminal_prot = []
            mods_variable_terminal_prot.append('a<^')

            base_config['fasta']['n_missed_cleavages'] = 2
            base_config['fasta']['pep_length_min'] = 7
            base_config['fasta']['pep_length_max'] = 50
            base_config['fasta']['mods_fixed'] = mods_fixed
            base_config['fasta']['mods_variable'] = variable_mods
            base_config['fasta']['mods_variable_terminal_prot'] = mods_variable_terminal_prot

            base_config['search']['top_n'] = 10
            base_config['search']['prec_tol'] = 20
            base_config['search']['frag_tol'] = 20

            base_config['raw']['n_most_abundant'] = 150

            base_config['features']['iso_charge_min'] = 2
            base_config['features']['iso_charge_max'] = 5

        return base_config
