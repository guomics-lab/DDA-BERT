import os

import numpy as np
import pandas as pd
import polars as pl


def gen_sage_psm(psm, candidate, out_dir, file_name):
    psm = psm[['scan', 'mz_array', 'intensity_array', 'precursor_mz']]
    psm.columns = ['scan_number', 'mz_array', 'intensity_array', 'precursor_mz']

    ## sage
    candidate = candidate[
        ['scan_number', 'charge', 'precursor_sequence', 'label', 'index', 'delta_rt_model', 'predicted_rt']].copy()
    candidate['scan_number'] = candidate['scan_number'].astype(np.int64)
    candidate.columns = ['scan_number', 'precursor_charge', 'modified_sequence', 'label', 'index', 'delta_rt_model',
                         'predicted_rt']
    candidate = pl.from_pandas(candidate)

    cols = ['scan_number', 'mz_array', 'intensity_array', 'precursor_mz']
    psm = psm.select(pl.col(cols))

    psm = psm.join(candidate, on='scan_number')
    psm = psm.filter(~pl.col("modified_sequence").is_null())

    final_out_file = os.path.join(out_dir, f'{file_name}.ipc')
    psm.write_ipc(final_out_file)


def proc_one_file(file_name, psm_file_path, candidate_dir, out_dir):
    print('parse: ', file_name)
    candidate_path = os.path.join(candidate_dir, f'{file_name}.csv')
    if os.path.exists(psm_file_path):
        psm = pl.read_ipc(psm_file_path)
        candidate = pd.read_csv(candidate_path)
        gen_sage_psm(psm, candidate, out_dir, file_name)
    else:
        print(f'file: {psm_file_path} dont exist!!!')
