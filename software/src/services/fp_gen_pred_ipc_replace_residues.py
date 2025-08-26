import os
from functools import reduce

import numpy as np
import pandas as pd
import polars as pl


def mkdir_p(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)

    return True, 'OK'


def gen_fp_psm(psm, candidate, final_out_file):
    psm = psm[['scan', 'mz_array', 'intensity_array', 'precursor_mz']]
    psm.columns = ['scan_number', 'mz_array', 'intensity_array', 'precursor_mz']

    required_cols = ['scan_number', 'charge', 'precursor_sequence', 'label', 'index',
                     'predicted_rt_1', 'aligned_rt', 'delta_rt_1']

    available_cols = [col for col in required_cols if col in candidate.columns]

    candidate = candidate[available_cols].copy()
    rename_map = {
        'charge': 'precursor_charge',
        'precursor_sequence': 'modified_sequence',
        'predicted_rt_1': 'predicted_rt',
        'delta_rt_1': 'delta_rt_model'
    }

    candidate = candidate.rename(columns=rename_map)
    candidate['scan_number'] = candidate['scan_number'].astype(np.int64)
    candidate = pl.from_pandas(candidate)

    psm = psm.join(candidate, on='scan_number')
    psm = psm.filter(~pl.col("modified_sequence").is_null())

    residues_frag = {
        'oxM': 'M[15.99]',
        'cC': 'C[57.02]',
        '^a': 'n[42]',
        'deamN': 'N[.98]',
        'deamQ': 'Q[.98]'
    }

    replacement_expr = reduce(
        lambda acc, item: acc.str.replace_all(item[0], item[1], literal=False),
        residues_frag.items(),
        pl.col("modified_sequence")
    )

    psm = psm.with_columns(
        replacement_expr.alias("modified_sequence")
    )

    print('psm: ', len(psm), 'final scan: ', len(set(psm['scan_number'])), 'sequence: ',
          len(set(psm['modified_sequence'])))

    # final_out_file = f'{output_dir}{file_name}.ipc'
    psm.write_ipc(final_out_file)


def proc(psm_file_path, candidate_file_path, final_out_file):
    psm = pl.read_ipc(psm_file_path)
    candidate_data = pd.read_csv(candidate_file_path)
    gen_fp_psm(psm, candidate_data, final_out_file)
