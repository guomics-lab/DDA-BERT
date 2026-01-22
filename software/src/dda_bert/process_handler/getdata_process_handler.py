import bisect
import json
import os
import shutil
import time
from pathlib import Path
import subprocess

import pickle

import numpy as np
import pandas as pd
import polars as pl
from pyteomics import mzml

from ..common import constant
from ..common.runtime_data_info import runtime_data
from ..obj.Entity import MzMLFileInfo, ProcessMsg
from ..process_handler.common_process_handler import CommonProcessHandler
from ..services import alphapept_gen_pred_ipc
from ..services import alphapept_hdf_to_csv_param
from ..services import alphapept_search_clean
from ..services import fp_gen_pred_ipc_replace_residues
from ..services import fp_search_clean
from ..services import gen_pred_ipc
from ..services import pred_rt_v4_3
from ..services import sage_search_clean
from alphapept.settings import load_settings


def clean_psm_func(peptide, residues_dict):
    for key, value in residues_dict.items():
        if value not in peptide:
            peptide = peptide.replace(key, value)
    return peptide


class GetDataProcessHandler(CommonProcessHandler):

    def __init__(self, f_info: MzMLFileInfo, mzml_raw_spec_abs_path, base_output_path, logger, env, input_param, start_time=0, open_sage=False,
                 open_fp=False, open_ap=False):
        CommonProcessHandler.__init__(self, f_info, base_output_path, logger, env, input_param, start_time,
                                      step=constant.ProgressStepEnum.GET_DATA)
        self.mzml_raw_spec_abs_path = mzml_raw_spec_abs_path

        self.open_sage = open_sage
        self.open_fp = open_fp
        self.open_ap = open_ap

    def deal_process(self):
        self.send_msg('Processing split sage result', status=constant.ProgressStepStatusEnum.RUNNING)
        pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.GET_DATA,
                        constant.ProgressStepStatusEnum.RUNNING)
        print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
        try:
            if not runtime_data.current_is_success:
                self.logger.error('Current is error')
                self.send_msg(status=constant.ProgressStepStatusEnum.ERROR)
                return
            ss = time.time()
            if not self.input_param.open_get_data:
                ee = time.time()
                self.send_msg('Finished build ipc data， spend time: {}'.format(self.get_use_time(ss, ee)),
                              status=constant.ProgressStepStatusEnum.SUCCESS)
                return

            '''********************************'''
            base_file_name = self.f_info.base_file_name
            psm_ipc_file_path = os.path.join(self.psm_output_path, f'{self.f_info.base_file_name}.ipc')
            scan_id_dict = None
            if self.is_d:
                self.run_d_file_spect(self.f_info.file_path)
                self.load_dfile_mzml_data(psm_ipc_file_path)
            elif self.is_wiff:
                _, scan_id_dict = self.load_wiff_mzml_data(self.f_info.file_path, psm_ipc_file_path)
                with open(os.path.join(self.psm_output_path, 'scan_id_dict.pkl'), mode='wb') as f:
                    pickle.dump(scan_id_dict, f)
            else:
                self.load_mzml_data(self.f_info.file_path, psm_ipc_file_path)

            self.deal_sage_data(scan_id_dict)

            if self.open_fp:
                self.deal_fragpipe_data()
            if self.open_ap:
                self.deal_alphapept_data()

            ee = time.time()
            self.send_msg('Finished build ipc data， spend time: {}'.format(self.get_use_time(ss, ee)),
                          status=constant.ProgressStepStatusEnum.SUCCESS)

        except Exception as e:
            runtime_data.current_is_success = False
            self.logger.exception('Build ipc data exception')
            self.send_msg('Build ipc data exception: {}'.format(e),
                          status=constant.ProgressStepStatusEnum.ERROR)
        finally:
            if runtime_data.current_is_success:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.GET_DATA,
                                constant.ProgressStepStatusEnum.SUCCESS)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))
            else:
                pm = ProcessMsg(runtime_data.current_mzml_index, constant.ProgressStepEnum.GET_DATA,
                                constant.ProgressStepStatusEnum.ERROR)
                print('ProcessMsg:{}'.format(json.dumps(pm.__dict__)))

    def deal_sage_data(self, scan_id_dict):
        #
        org_result_file_path = os.path.join(self.sage_result_output_path, 'results.sage.tsv')
        self.send_msg('Finished split sage result')
        base_file_name = self.f_info.base_file_name
        sage_result_file_path = os.path.join(self.sage_result_output_path, f'{base_file_name}.sage.tsv')
        shutil.copy(org_result_file_path, sage_result_file_path)

        self.logger.info(f'Process sage_first_clean, result dir: {self.sage_first_clean_dir}')
        sage_search_clean.sage_first_clean_one_file(base_file_name, sage_result_file_path, self.sage_first_clean_dir, self.is_wiff, scan_id_dict)
        self.logger.info(f'Finish sage_first_clean, result dir: {self.sage_first_clean_dir}')

        psm_ipc_file_path = os.path.join(self.psm_output_path, f'{self.f_info.base_file_name}.ipc')
        self.logger.info(f'Process gen_pred_ipc {base_file_name}, result dir: {self.sage_gen_pred_ipc_result_dir}')
        gen_pred_ipc.proc_one_file(base_file_name, psm_ipc_file_path, self.sage_first_clean_dir,
                                   self.sage_gen_pred_ipc_result_dir)
        self.logger.info(f'Finish gen_pred_ipc {base_file_name}, result dir: {self.sage_gen_pred_ipc_result_dir}')

    def deal_alphapept_data(self):
        base_file_name = self.f_info.base_file_name
        try:
            run_alphapept_config_path = os.path.join(self.alphapept_result_output_path,
                                                     '{}_alphapept_config.yaml'.format(base_file_name))
            settings = load_settings(run_alphapept_config_path)
            min_frag_hits = settings['search']['min_frag_hits']
            self.logger.info(f'load ap set file, get min_frag_hits config is: {min_frag_hits}')
        except Exception:
            self.logger.exception('load ap set file error, get default min_frag_hits config.')
            min_frag_hits = 5

        #
        hdf_file_path = os.path.join(self.alphapept_result_output_path, f'{base_file_name}.ms_data.hdf')
        alphapept_csv_path = os.path.join(self.alphapept_result_output_path, f'{base_file_name}.csv')
        self.logger.info(
            f'Process alphapept_hdf_to_csv {base_file_name}, hdf_file: {hdf_file_path}, result dir: {self.alphapept_result_output_path}')
        alphapept_hdf_to_csv_param.hdf_to_csv(hdf_file_path, alphapept_csv_path, min_frag_hits)
        self.logger.info(
            f'Finish alphapept_hdf_to_csv {base_file_name}, result dir: {self.alphapept_result_output_path}')

        sage_result_file_path = os.path.join(self.sage_result_output_path, f'{base_file_name}.sage.tsv')
        self.logger.info(
            f'Process pred_rt_v4_3 {base_file_name}, base_sage_result: {self.sage_result_output_path}, result dir: {self.alphapept_result_output_path}')
        pred_rt_v4_3.deal_process_with_exception(sage_result_file_path, alphapept_csv_path,
                                                 self.alphapept_result_output_path)
        self.logger.info(f'Finish pred_rt_v4_3 {base_file_name}, result dir: {self.alphapept_result_output_path}')

        with_fitting_csv = os.path.join(self.alphapept_result_output_path, f'{base_file_name}.with_fitting.csv')
        self.logger.info(
            f'Process alphapept_search_clean {base_file_name}, with_fitting_csv: {with_fitting_csv}, result dir: {self.alphapept_data_clean_output_dir}')
        alphapept_search_clean.alphapept_first_clean_one(base_file_name, with_fitting_csv,
                                                         self.alphapept_data_clean_output_dir)
        self.logger.info(
            f'Finish alphapept_search_clean {base_file_name}, result dir: {self.alphapept_data_clean_output_dir}')

        psm_file_path = os.path.join(self.psm_output_path, f'{self.f_info.base_file_name}.ipc')
        self.logger.info(
            f'Process alphapept_gen_pred_ipc {base_file_name}, psm_file_path: {psm_file_path}, result dir: {self.alphapept_gen_pred_output_dir}')
        alphapept_gen_pred_ipc.proc_one(psm_file_path, base_file_name, self.alphapept_data_clean_output_dir,
                                        self.alphapept_gen_pred_output_dir)
        self.logger.info(
            f'Finish alphapept_gen_pred_ipc {base_file_name}, result dir: {self.alphapept_gen_pred_output_dir}')

    def deal_fragpipe_data(self):
        base_file_name = self.f_info.base_file_name
        pin_file_path = os.path.join(self.fragpipe_output_path, 'exp', f'{base_file_name}.pin')
        clean_file_path = os.path.join(self.fragpipe_fp_first_clean_dir, f'{base_file_name}.csv')
        self.logger.info(
            f'Process fp_frist_clean: {pin_file_path}, clean_file_path: {clean_file_path}')
        fp_search_clean.fp_frist_clean(pin_file_path, clean_file_path)
        self.logger.info(f'Finish fp_frist_clean')

        sage_result_file_path = os.path.join(self.sage_result_output_path, f'{base_file_name}.sage.tsv')
        self.logger.info(
            f'Process fragpipe pred_rt_v4_3 {base_file_name}, result dir: {self.fragpipe_fp_first_clean_dir}')
        pred_rt_v4_3.deal_process_with_exception(sage_result_file_path, clean_file_path,
                                                 self.fragpipe_fp_first_clean_dir)
        self.logger.info(f'Finish fragpipe pred_rt_v4_3')

        psm_file_path = os.path.join(self.psm_output_path, f'{self.f_info.base_file_name}.ipc')
        fp_clean_fitting_file_path = os.path.join(self.fragpipe_fp_first_clean_dir,
                                                  f'{base_file_name}.with_fitting.csv')
        fp_proc_file_path = os.path.join(self.fragpipe_fp_first_clean_dir, f'{base_file_name}.ipc')
        self.logger.info(f'Process fragpipe proc')
        fp_gen_pred_ipc_replace_residues.proc(psm_file_path, fp_clean_fitting_file_path, fp_proc_file_path)
        self.logger.info(f'Finish fragpipe proc, result file: {fp_proc_file_path}')


    def run_d_file_spect(self, mzml_path):
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


    def load_dfile_mzml_data(self, target):
        #
        raw_spect_output_path = os.path.join(self.sage_result_output_path,
                                             '{}_raw_spec.tsv'.format(self.f_info.base_file_name))
        spec_df = pd.read_csv(raw_spect_output_path, sep='\t')

        self.logger.info('Finished load mzml data')
        self.send_msg('Finished load mzml data, path is {}'.format(raw_spect_output_path))

        spec_df['mz_array'] = spec_df['mz_array'].apply(lambda x: [float(dd) for dd in x.split(',')])
        spec_df['intensity_array'] = spec_df['intensity_array'].apply(lambda x: [float(dd) for dd in x.split(',')])

        cols = ['scan', 'ms1_scan', 'precursor_mz', 'mz_array', 'intensity_array']

        # save result
        dl = pl.from_pandas(spec_df[cols])
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        dl.write_ipc(target)
        return spec_df[cols]

    def load_wiff_mzml_data(self, source: str, target: str):
        """Load data from an mzml file as a dictionary.

        Args:
            source (Path | str): The name or path of .mzML
            target (int): the filename of output.

        """

        try:
            reader = mzml.read(source, use_index=True)
            spec_indices = np.array(range(1, len(reader) + 1))
        except OSError:
            self.logger.info('Could not open the file. Please, specify the correct path to the file.')
            return

        scan_list = []
        ms1_scan_list = []
        ms_order_list = []
        mz_list = []
        intentisy_list = []
        precursor_mz_list = []

        scan_id_dict = {}

        pre_ms1_scan = 0

        for idx, _ in enumerate(spec_indices):
            try:
                spec = next(reader)
                scan = spec.get('index')
                scan_id = spec.get('id')
                scan_id_dict[str(scan_id)] = scan

                rt, mzs, intensities, ms_order, precursor_mz = self.extract_mzml_info(spec)

                # 把scan替换为index

                if ms_order == 1:
                    pre_ms1_scan = scan

                ms_order_list.append(ms_order)

                #
                scan_list.append(scan)
                #  scan
                ms1_scan_list.append(pre_ms1_scan)

                # precursor moz
                precursor_mz_list.append(precursor_mz)

                #
                sortindex = np.argsort(mzs)
                mzs = mzs[sortindex]
                intensities = intensities[sortindex]

                # Remove zero intensities
                to_keep = intensities > 0
                mzs = mzs[to_keep]
                intensities = intensities[to_keep]

                #
                mz_list.append(np.round(mzs, 6))
                intentisy_list.append(np.round(intensities, 6))

            except KeyboardInterrupt as e:
                raise e
            except SystemExit as e:
                raise e
            except Exception as e:
                self.logger.info(f"Bad scan={scan} in mzML file '{source}' {e}")

        #
        query_data = {}
        query_data["scan"] = scan_list
        query_data["ms1_scan"] = ms1_scan_list
        query_data["ms_level"] = ms_order_list
        query_data["mz_array"] = mz_list  #
        query_data["intensity_array"] = intentisy_list  #
        query_data["precursor_mz"] = precursor_mz_list

        #
        df = pd.DataFrame.from_dict(query_data)

        # split ms1 and ms2
        ms1_df = df[df['ms_level'] == 1]
        ms2_df = df[df['ms_level'] == 2]

        for col in ['scan']:
            ms2_df[col] = ms2_df[col].astype(int)

        # calc ms1 spectrum
        ms2_df['precursor_mz_array'] = ms2_df['ms1_scan'].map(ms1_df.set_index('scan')['mz_array'])
        ms2_df['precursor_intensity_array'] = ms2_df['ms1_scan'].map(ms1_df.set_index('scan')['intensity_array'])

        ms2_df['ms1_spectrum'] = ms2_df.apply(lambda row: parse_spectrum(row['precursor_mz'], row['precursor_mz_array'],
                                                                         row['precursor_intensity_array']), axis=1)
        ms2_df['ms1_mz_array'] = ms2_df['ms1_spectrum'].apply(lambda x: x[0])
        ms2_df['ms1_intensity_array'] = ms2_df['ms1_spectrum'].apply(lambda x: x[1])

        # combine ms1 and ms2 spectrum
        ms2_df['combine_mz_array'] = ms2_df.apply(lambda x: x['ms1_mz_array'].tolist() + x['mz_array'].tolist(), axis=1)
        ms2_df['combine_intensity_array'] = ms2_df.apply(
            lambda x: x['ms1_intensity_array'].tolist() + x['intensity_array'].tolist(), axis=1)

        cols = ['scan', 'ms1_scan', 'precursor_mz', 'mz_array', 'intensity_array', 'combine_mz_array',
                'combine_intensity_array']

        # save result
        dl = pl.from_pandas(ms2_df[cols])
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        dl.write_ipc(target)
        return ms2_df[cols], scan_id_dict


    def load_mzml_data(self, source: str, target: str):
        """Load data from an mzml file as a dictionary.

        Args:
            source (Path | str): The name or path of .mzML
            target (int): the filename of output.

        """

        try:
            reader = mzml.read(source, use_index=True)
            spec_indices = np.array(range(1, len(reader) + 1))
        except OSError:
            self.logger.info('Could not open the file. Please, specify the correct path to the file.')
            return

        scan_list = []
        # spectrum_scan_list = []
        ms1_scan_list = []
        ms_order_list = []
        mz_list = []
        intentisy_list = []
        precursor_mz_list = []

        pre_ms1_scan = 0

        for idx, scan in enumerate(spec_indices):
            try:
                spec = next(reader)

                rt, mzs, intensities, ms_order, precursor_mz = self.extract_mzml_info(spec)

                if ms_order == 1:
                    pre_ms1_scan = scan

                ms_order_list.append(ms_order)

                #
                scan_list.append(scan)
                #
                # spectrum_scan_list.append(spectrum_scan)
                #  scan
                ms1_scan_list.append(pre_ms1_scan)

                # precursor moz
                precursor_mz_list.append(precursor_mz)

                #
                sortindex = np.argsort(mzs)
                mzs = mzs[sortindex]
                intensities = intensities[sortindex]

                # Remove zero intensities
                to_keep = intensities > 0
                mzs = mzs[to_keep]
                intensities = intensities[to_keep]

                #
                mz_list.append(np.round(mzs, 6))
                intentisy_list.append(np.round(intensities, 6))

            except KeyboardInterrupt as e:
                raise e
            except SystemExit as e:
                raise e
            except Exception as e:
                self.logger.info(f"Bad scan={scan} in mzML file '{source}' {e}")

        #
        query_data = {}
        query_data["scan"] = scan_list
        # query_data["spectrum_scan"] = spectrum_scan_list
        query_data["ms1_scan"] = ms1_scan_list
        query_data["ms_level"] = ms_order_list
        query_data["mz_array"] = mz_list  #
        query_data["intensity_array"] = intentisy_list  #
        query_data["precursor_mz"] = precursor_mz_list

        #
        df = pd.DataFrame.from_dict(query_data)

        # split ms1 and ms2
        ms1_df = df[df['ms_level'] == 1]
        ms2_df = df[df['ms_level'] == 2]

        for col in ['scan']:
            ms2_df[col] = ms2_df[col].astype(int)
        # assert ms2_df['scan'].equals(ms2_df['spectrum_scan'])

        # calc ms1 spectrum
        ms2_df['precursor_mz_array'] = ms2_df['ms1_scan'].map(ms1_df.set_index('scan')['mz_array'])
        ms2_df['precursor_intensity_array'] = ms2_df['ms1_scan'].map(ms1_df.set_index('scan')['intensity_array'])

        ms2_df['ms1_spectrum'] = ms2_df.apply(lambda row: parse_spectrum(row['precursor_mz'], row['precursor_mz_array'],
                                                                         row['precursor_intensity_array']), axis=1)
        ms2_df['ms1_mz_array'] = ms2_df['ms1_spectrum'].apply(lambda x: x[0])
        ms2_df['ms1_intensity_array'] = ms2_df['ms1_spectrum'].apply(lambda x: x[1])

        # combine ms1 and ms2 spectrum
        ms2_df['combine_mz_array'] = ms2_df.apply(lambda x: x['ms1_mz_array'].tolist() + x['mz_array'].tolist(), axis=1)
        ms2_df['combine_intensity_array'] = ms2_df.apply(
            lambda x: x['ms1_intensity_array'].tolist() + x['intensity_array'].tolist(), axis=1)

        cols = ['scan', 'ms1_scan', 'precursor_mz', 'mz_array', 'intensity_array', 'combine_mz_array',
                'combine_intensity_array']

        # save result
        dl = pl.from_pandas(ms2_df[cols])
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        dl.write_ipc(target)
        return ms2_df[cols]

    def extract_mzml_info(self, input_dict: dict) -> tuple:
        """Extract basic MS coordinate arrays from a dictionary.

        Args:
            input_dict (dict): A dictionary obtained by iterating over a Pyteomics mzml.read function.

        Returns:
            tuple: The rt, masses, intensities, ms_order, prec_mass, mono_mz, charge arrays retrieved from the input_dict.
                If the `ms level` in the input dict does not equal 2, the charge, mono_mz and prec_mass will be equal to 0.

        """
        rt = float(input_dict.get('scanList').get('scan')[0].get('scan start time'))  # rt_list_ms1/2
        masses = input_dict.get('m/z array')
        intensities = input_dict.get('intensity array')
        ms_order = input_dict.get('ms level')
        # spectrum_scan = input_dict.get('spectrum title').split('.')[1]

        precursor_mz = 0
        if ms_order == 2:
            precursor_mz = \
                input_dict.get('precursorList').get('precursor')[0].get('selectedIonList').get('selectedIon')[
                    0].get(
                    'selected ion m/z')
        elif ms_order == 1:
            precursor_mz = input_dict.get('base peak m/z')

        return rt, masses, intensities, ms_order, precursor_mz


def parse_spectrum(precursor_mz, precursor_mz_array, precursor_intensity_array, mz_unit='Da', mz_tol=10,
                   padding_len=50):
    if precursor_mz == -1:
        return [0 for _ in range(len(precursor_mz_array))]

    if mz_unit == "Da":
        extract_width = [precursor_mz - mz_tol / 2, precursor_mz + mz_tol / 2]
    elif mz_unit == "ppm":
        mz_tol_da = precursor_mz * mz_tol * 0.000001
        extract_width = [precursor_mz - mz_tol_da / 2, precursor_mz + mz_tol_da / 2]

    extract_left_index = bisect.bisect_left(precursor_mz_array, extract_width[0])
    extract_right_index = bisect.bisect_right(precursor_mz_array, extract_width[1])

    extract_mz = np.array(precursor_mz_array[extract_left_index: extract_right_index])
    extract_intensity = np.array(precursor_intensity_array[extract_left_index: extract_right_index])

    if len(extract_intensity) >= padding_len:
        #
        topn_index = extract_intensity.argsort()[-padding_len:][::-1]
        extract_intensity = extract_intensity[topn_index]
        extract_mz = extract_mz[topn_index]
    else:
        #
        topn_index = extract_intensity.argsort()[::-1]
        extract_intensity = extract_intensity[topn_index]
        extract_mz = extract_mz[topn_index]

        pad_array = np.array([0 for _ in range(padding_len - len(extract_intensity))])
        extract_intensity = np.append(extract_intensity, pad_array)
        extract_mz = np.append(extract_mz, pad_array)

    return (extract_mz, extract_intensity)
