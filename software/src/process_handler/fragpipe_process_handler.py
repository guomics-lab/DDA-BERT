import json
import os.path
import subprocess
import time

from src.common import constant
from src.common.runtime_data_info import runtime_data
from src.obj.Entity import MzMLFileInfo, ProcessMsg
from src.process_handler.common_process_handler import CommonProcessHandler


class FragpipeProcessHandler(CommonProcessHandler):

    def __init__(self, exe_abs_path, f_info: MzMLFileInfo, fasta_path, base_output_path, logger,
                 env, input_param,
                 start_time=0):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.SAGA)
        self.exe_abs_path = exe_abs_path
        self.fasta_path = fasta_path
        self.base_work_flow_file_path = 'config/FP_MSBooster.workflow'

    def deal_process(self):
        pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.SAGA,
                        constant.ProgressStepStatusEnum.RUNNING)
        print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
        self.send_msg('Processing fragpipe', status=constant.ProgressStepStatusEnum.RUNNING)
        try:
            ss = time.time()
            if not self.input_param.open_sage:
                ee = time.time()
                self.send_msg(msg='Finished fragpipe, spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.SUCCESS)
                return

            # build decoy fasta
            decoy_fast_path = self.build_fasta_with_decoy()
            manifest_path = self.build_manifest()
            workflow_path = self.build_workflow(decoy_fast_path)
            self.run_fragpipe(manifest_path, workflow_path)
            ee = time.time()

            self.send_msg(msg='Finished fragpipe, spend time: {}'.format(self.get_use_time(ss, ee)),
                          status=constant.ProgressStepStatusEnum.SUCCESS)

        except Exception:
            runtime_data.current_is_success = False
            self.logger.exception('Processing fragpipe exception')
            self.send_msg(msg='fragpipe exception', status=constant.ProgressStepStatusEnum.ERROR)
        finally:
            if runtime_data.current_is_success:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.SAGA,
                                constant.ProgressStepStatusEnum.SUCCESS)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
            else:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.SAGA,
                                constant.ProgressStepStatusEnum.ERROR)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))

    def build_fasta_with_decoy(self):

        philosopher_exe_path = os.path.join(os.path.join(os.getcwd(), 'FragPipe22_0', 'philosopher-v5.1.1'))
        cmd2 = [philosopher_exe_path, 'workspace', '--init']
        self.run_cmd(cmd2, cwd=self.fragpipe_output_path)
        cmd3 = [philosopher_exe_path, 'database', '--custom', self.fasta_path, '--contam', '--contamprefix', 'rev_']
        self.run_cmd(cmd3, cwd=self.fragpipe_output_path)
        cmd4 = [philosopher_exe_path, 'workspace', '--clean']
        self.run_cmd(cmd4, cwd=self.fragpipe_output_path)

        all_file_list = os.listdir(self.fragpipe_output_path)
        decoy_fasta_file_name_list = [file_name for file_name in all_file_list if
                                      file_name.endswith(f'-decoys-contam-{os.path.basename(self.fasta_path)}.fas')]
        decoy_fasta_file_name = decoy_fasta_file_name_list[-1]
        decoy_fast_path = os.path.join(self.fragpipe_output_path, decoy_fasta_file_name)
        self.logger.info(f'Fasta with decoy: {decoy_fast_path}')
        return decoy_fast_path

    def run_cmd(self, cmd, cwd=None):
        env = os.environ.copy()
        java_bin_path = os.path.join(os.getcwd(), 'jdk-11.0.26', 'bin')
        env['PATH'] = java_bin_path + os.pathsep + env['PATH']
        java_home = os.path.join(os.getcwd(), 'jdk-11.0.26')
        env['JAVA_HOME'] = java_home

        p = subprocess.Popen(cmd, env=env, cwd=cwd, stdout=subprocess.PIPE)
        stdout = p.stdout
        while True:
            output = stdout.readline()
            if output == b'' or (output == '' and p.poll() is not None):
                break
            self.logger.info(output)
            if output:
                #
                info_msg = output.decode('utf-8')
                info_msg = info_msg.rstrip()
                if len(info_msg) == 0:
                    continue

    '''
    
    '''

    def build_manifest(self):
        manifest_path = os.path.join(self.fragpipe_output_path, 'fragpipe-files.fp-manifest')
        self.send_msg('Processing build fragpipe manifest, path is {}'.format(manifest_path))
        #
        with open(manifest_path, 'w+') as f:
            f.write(f'{self.f_info.file_path}\texp\t\tDDA')
        self.send_msg('Finished build fragpipe manifest')
        return manifest_path

    '''
    
    '''

    def build_workflow(self, decoy_fast_path):
        workflow_path = os.path.join(self.fragpipe_output_path, 'fragpipe.workflow')
        self.send_msg('Processing build fragpipe workflow, path is {}'.format(workflow_path))
        with open(self.base_work_flow_file_path, 'r', encoding='utf-8') as fin:
            content = fin.read()

        with open(workflow_path, 'w+', encoding='utf-8') as fout:
            fout.write(f'database.db-path={decoy_fast_path}\n')
            fout.write(content)
        self.send_msg('Finished build fragpipe workflow')
        return workflow_path

    '''
    
    '''

    def run_fragpipe(self, manifest_path, workflow_path):
        if self.env == 'win':
            # p = subprocess.Popen(cmd, stdout=subprocess.PIPE, creationflags=134217728)
            pass
        elif self.env == 'linux':
            frag_base_dir = os.path.dirname(os.path.dirname(self.exe_abs_path))
            cmd = [self.exe_abs_path, '--headless', '--workflow', workflow_path, '--manifest',
                   manifest_path, '--workdir', self.fragpipe_output_path,
                   '--config-tools-folder', frag_base_dir, '--threads', str(self.input_param.thread_num)]

            self.send_msg('Processing run fragpipe, command is {}'.format(cmd))
            self.run_cmd(cmd)
        self.send_msg('Finished run fragpipe, result path is {}'.format(self.fragpipe_output_path))
