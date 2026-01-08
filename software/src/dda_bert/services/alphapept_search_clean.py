# pip install peptdeep;import argparse
import os
import re

import numpy as np
import pandas as pd


def mkdir_p(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)

    return True, 'OK'


residues_alphapept = {
    '-': 'C[57.02]',
    '=': 'M[15.99]',
    '*': 'n[42]',
    '.': 'N[.98]',
    '+': 'Q[.98]'
}


def clean_psm_func(peptide: str, residues_dict: dict) -> str:
    peptide = re.sub(r'(?<!\d)\.(?!\d)', residues_dict['.'], peptide)

    for key, value in residues_dict.items():
        if key != '.':
            peptide = peptide.replace(key, value)
    return peptide


def clean_alphapept_frist(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path)

    df['precursor_sequence'] = df['sequenceNew'].apply(lambda x: clean_psm_func(x, residues_alphapept))

    df['label'] = df['label'].astype(int)
    df['charge'] = df['charge'].astype(int)

    if 'scan_no' in df.columns:
        df.rename(columns={'scan_no': 'scan_number'}, inplace=True)
        df['scan_number'] = df['scan_number'].astype(np.int32)

    df = df.rename(columns={
        'predicted_rt_1': 'predicted_rt',
        'delta_rt_1': 'delta_rt_model'
    })

    df['cleaned_sequence'] = df['precursor_sequence']
    df['cleaned_sequence'] = df['cleaned_sequence'].str.replace('n[42]', '', regex=False)
    df['cleaned_sequence'] = df['cleaned_sequence'].str.replace('N[.98]', 'N', regex=False)
    df['cleaned_sequence'] = df['cleaned_sequence'].str.replace('Q[.98]', 'Q', regex=False)
    df['cleaned_sequence'] = df['cleaned_sequence'].str.replace('M[15.99]', 'M', regex=False)
    df['cleaned_sequence'] = df['cleaned_sequence'].str.replace('C[57.02]', 'C', regex=False)

    df['sequence_len'] = df['cleaned_sequence'].apply(len)

    df['precursor_id'] = df['charge'].astype(str) + '_' + df['precursor_sequence']
    # df['psm_id'] = df['scan_number'].astype(str) + '_' + df['precursor_id']

    if 'scan_number' in df.columns:
        df['psm_id'] = df['scan_number'].astype(str) + '_' + df['precursor_id']

    df = df[(df['sequence_len'] <= 50) & (df['sequence_len'] >= 7)]
    df = df[(df['charge'] <= 5) & (df['charge'] >= 2)]

    df = df.reset_index(drop=True)
    df['index'] = range(len(df))

    return df


def alphapept_first_clean_one(base_file_name, file, alphapept_first_clean_dir):
    data = clean_alphapept_frist(file)
    data.to_csv(os.path.join(alphapept_first_clean_dir, f'{base_file_name}.csv'), index=False)
