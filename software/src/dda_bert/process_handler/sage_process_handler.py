import json
import os.path
import subprocess
import time

from ..common import constant
from ..common.runtime_data_info import runtime_data
from ..obj.Entity import MzMLFileInfo, ProcessMsg
from ..process_handler.common_process_handler import CommonProcessHandler


class SageProcessHandler(CommonProcessHandler):

    def __init__(self, exe_abs_path, mzml_raw_spec_abs_path, f_info: MzMLFileInfo, fasta_path, sage_config_file,
                 base_output_path, logger,
                 env, input_param,
                 start_time=0):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.SAGA)
        self.exe_abs_path = exe_abs_path
        self.mzml_raw_spec_abs_path = mzml_raw_spec_abs_path
        self.fasta_path = fasta_path
        self.sage_config_file = sage_config_file

    def deal_process(self):
        pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.SAGA,
                        constant.ProgressStepStatusEnum.RUNNING)
        print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
        self.send_msg('Processing sage', status=constant.ProgressStepStatusEnum.RUNNING)
        try:
            ss = time.time()
            if not self.input_param.open_sage:
                ee = time.time()
                self.send_msg(msg='Finished sage, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.SUCCESS)
                return
            sage_config_path = self.build_config(self.sage_config_file)
            self.run_sage(sage_config_path)
            # if self.is_d:
            #     self.run_mzml_raw_spect(self.f_info.file_path)
            sage_result_file_path = os.path.join(self.sage_result_output_path, 'results.sage.tsv')
            ee = time.time()
            if not os.path.exists(sage_result_file_path):
                self.send_msg(msg='Finished sage, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.ERROR)
                runtime_data.current_is_success = False
                return
            self.send_msg(msg='Finished sage, spend time: {}'.format(self.get_use_time(ss, ee)),
                          status=constant.ProgressStepStatusEnum.SUCCESS)

        except Exception:
            runtime_data.current_is_success = False
            self.logger.exception('Processing sage exception')
            self.send_msg(msg='Sage exception', status=constant.ProgressStepStatusEnum.ERROR)
        finally:
            if runtime_data.current_is_success:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.SAGA,
                                constant.ProgressStepStatusEnum.SUCCESS)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
            else:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.SAGA,
                                constant.ProgressStepStatusEnum.ERROR)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))

    def build_config(self, custom_config_path=None):
        current_mzml_base_name = self.f_info.base_file_name
        sage_config_path = os.path.join(self.sage_output_path, '{}_sage_config.json'.format(current_mzml_base_name))
        self.send_msg('Processing build sage config, path is {}'.format(sage_config_path))
        with open(sage_config_path, 'w+') as f:
            json.dump(self.get_config_dict(custom_config_path), f)
        self.send_msg('Finished build sage config')
        return sage_config_path

    def get_config_dict(self, custom_config_path=None):
        mzml_paths = [self.f_info.file_path]
        mzml_output_path = self.sage_result_output_path

        if custom_config_path:
            config_json_dict = json.load(open(custom_config_path))
        else:
            self.send_msg('Custom sage config not set, use default config')
            config_json_dict = self.get_default_dict()

        config_json_dict['database']['fasta'] = self.fasta_path

        config_json_dict['mzml_paths'] = mzml_paths
        config_json_dict['output_directory'] = mzml_output_path

        return config_json_dict

    def get_default_dict(self):

        static_mods_dict = {}
        static_mods_dict['C'] = 57.0216

        variable_mods_dict = {}
        variable_mods_dict["["] = [42.0]
        variable_mods_dict["M"] = [15.9949]
        variable_mods_dict["N"] = [0.98]
        variable_mods_dict["Q"] = [0.98]

        config_json_dict = {
            "database": {
                "bucket_size": 8192,
                "fragment_min_mz": 150.0,
                "fragment_max_mz": 1500.0,
                "enzyme": {
                    "missed_cleavages": 2,
                    "min_len": 7,
                    "max_len": 50,
                    "cleave_at": "KR",
                    "restrict": "P"
                },
                "static_mods": static_mods_dict,
                "variable_mods": variable_mods_dict,
                "max_variable_mods": 3,
                "decoy_tag": "rev_",
                "generate_decoys": True
            },
            "deisotope": True,
            "chimera": False,
            "min_peaks": 15,
            "max_peaks": 150,
            "min_matched_peaks": 4,
            "max_fragment_charge": 2,
            "report_psms": 10,
            "precursor_tol": {
                "ppm": [
                    -20,
                    20
                ]
            },
            "fragment_tol": {
                "ppm": [
                    -20,
                    20
                ]
            },
            "precursor_charge": [2, 5],
            "isotope_errors": [
                0,
                2
            ]
        }
        return config_json_dict

    def run_sage(self, sage_config_path):
        run_env = os.environ.copy()
        if self.env == 'win':
            cmd = '{} {}'.format(self.exe_abs_path, sage_config_path)
            self.send_msg('Processing run sage, command is {}'.format(cmd))
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, creationflags=134217728)
        elif self.env == 'linux':
            run_env['RUST_MIN_STACK'] = '8388608'
            cmd = [self.exe_abs_path, sage_config_path]
            self.send_msg('Processing run sage, command is {}'.format(cmd))
            p = subprocess.Popen(cmd, env=run_env, stdout=subprocess.PIPE)

        stdout = p.stdout
        while True:
            output = stdout.readline()
            if output == b'' or (output == '' and p.poll() is not None):
                break
            self.logger.info(output)
            if output:
                info_msg = output.decode('utf-8')
                info_msg = info_msg.rstrip()
                if len(info_msg) == 0:
                    continue
        mzml_output_path = os.path.join(self.sage_result_output_path,
                                        self.f_info.base_file_name)
        self.send_msg('Finished run sage, result path is {}'.format(mzml_output_path))

    def run_mzml_raw_spect(self, mzml_path):
        raw_spect_output_path = os.path.join(self.sage_result_output_path,
                                             '{}_raw_spec.tsv'.format(self.f_info.base_file_name))
        if self.env == 'win':
            cmd = '{} "{}" "{}" '.format(self.mzml_raw_spec_abs_path, mzml_path, raw_spect_output_path)
            # self.logger.info('Processing run raw_spect, command is {}'.format(cmd))
            self.send_msg('Processing run raw_spect, command is {}'.format(cmd))
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, creationflags=134217728)
        elif self.env == 'linux':
            cmd = [self.mzml_raw_spec_abs_path, mzml_path, raw_spect_output_path]
            # self.logger.info('Processing run raw_spect, command is {}'.format(cmd))
            self.send_msg('Processing run raw_spect, command is {}'.format(cmd))
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        stdout = p.stdout
        while True:
            output = stdout.readline()
            if output == b'' or (output == '' and p.poll() is not None):
                break
            # self.logger.info(output)
            if output:
                info_msg = output.decode('utf-8')
                info_msg = info_msg.rstrip()
                if len(info_msg) == 0:
                    continue

        self.send_msg('Finished Raw spect, result path is {}'.format(raw_spect_output_path))
