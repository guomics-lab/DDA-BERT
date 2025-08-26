import os

import numpy as np
import pandas as pd
import polars as pl


def mkdir_p(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)

    return True, 'OK'


def gen_alphapept_psm(psm, candidate, output_dir, file_name):
    psm = psm[['scan', 'mz_array', 'intensity_array', 'precursor_mz']]
    psm.columns = ['scan_number', 'mz_array', 'intensity_array', 'precursor_mz']

    #
    ## alphapept
    candidate = candidate[
        ['scan_number', 'charge', 'precursor_sequence', 'label', 'index', 'predicted_rt', 'delta_rt_model']].copy()
    candidate['scan_number'] = candidate['scan_number'].astype(np.int64)
    candidate.columns = ['scan_number', 'precursor_charge', 'modified_sequence', 'label', 'index', 'predicted_rt',
                         'delta_rt_model']
    candidate = pl.from_pandas(candidate)

    #
    cols = ['scan_number', 'mz_array', 'intensity_array', 'precursor_mz']
    psm = psm.select(pl.col(cols))

    #
    psm = psm.join(candidate, on='scan_number')
    psm = psm.filter(~pl.col("modified_sequence").is_null())
    print('psm: ', len(psm), 'final scan: ', len(set(psm['scan_number'])), 'sequence: ',
          len(set(psm['modified_sequence'])))

    final_out_file = os.path.join(output_dir, f'{file_name}.ipc')
    psm.write_ipc(final_out_file)


def proc_one(psm_file_path, base_file_name, alphapept_data_clean_output_dir, output_dir):
    # file_path = os.path.join(psm_dir, f'{base_file_name}.ipc')
    if os.path.exists(psm_file_path):
        psm = pl.read_ipc(psm_file_path)
        candidate = pd.read_csv(os.path.join(alphapept_data_clean_output_dir, f'{base_file_name}.csv'))
        gen_alphapept_psm(psm, candidate, output_dir, base_file_name)
    else:
        print(f'file: {psm_file_path} dont exist!!!')
