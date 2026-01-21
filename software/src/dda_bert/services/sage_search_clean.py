import numpy as np
import os
import pandas as pd

residues_sage = {
    'C[+57.0216]': 'C[57.02]',
    'M[+15.9949]': 'M[15.99]',
    '[+42]-': 'n[42]',
    'N[+0.98]': 'N[.98]',
    'Q[+0.98]': 'Q[.98]'
}


def clean_psm_func(peptide, residues_dict):
    for key, value in residues_dict.items():
        if value not in peptide:
            peptide = peptide.replace(key, value)
    return peptide


def clean_sage_frist(file_path):
    df = pd.read_table(file_path)
    df = df.rename(columns={'label': 'label_source', 'psm_id': 'psm_id_source'})

    df['precursor_sequence'] = df['peptide'].apply(lambda x: clean_psm_func(x, residues_sage))
    df.loc[df['label_source'] == 1, 'label'] = 1
    df.loc[df['label_source'] == -1, 'label'] = 0
    df['label'] = df['label'].astype(int)

    df['charge'] = df['charge'].astype(int)

    try:
        df['scan_number'] = df['scannr'].apply(lambda x: x.split('=')[-1]).astype(int)
    except:
        ##For tims
        df['scan_number'] = df['scannr'].astype(np.int32)

    df['scan_number'] = df['scan_number'].astype(np.int32)

    #     # M C n N Q
    df['cleaned_sequence'] = df['precursor_sequence'].str.replace('n[42]', '')
    df['cleaned_sequence'] = df['cleaned_sequence'].str.replace('N[.98]', 'N').str.replace('Q[.98]', 'Q').str.replace(
        'M[15.99]', 'M').str.replace('C[57.02]', 'C')
    df['sequence_len'] = df['cleaned_sequence'].apply(len)

    df['precursor_id'] = df['charge'].astype(str) + '_' + df['precursor_sequence']
    df['psm_id'] = df['scan_number'].astype(str) + '_' + df['precursor_id']

    df = df[(df['sequence_len'] <= 50) & (df['sequence_len'] >= 7)]

    # charge in [2,5]
    df = df[(df['charge'] <= 5) & (df['charge'] >= 2)]

    df = df.reset_index()

    return df


def clean_sage_frist_wiff(file_path, scan_id_dict):
    df = pd.read_table(file_path)
    df = df.rename(columns={'label': 'label_source', 'psm_id': 'psm_id_source'})

    df['precursor_sequence'] = df['peptide'].apply(lambda x: clean_psm_func(x, residues_sage))
    df.loc[df['label_source'] == 1, 'label'] = 1
    df.loc[df['label_source'] == -1, 'label'] = 0
    df['label'] = df['label'].astype(int)

    df['charge'] = df['charge'].astype(int)

    df['scan_number'] = df['scannr'].apply(lambda x: scan_id_dict[x]).astype(int)

    df['scan_number'] = df['scan_number'].astype(np.int32)

    #     # M C n N Q
    df['cleaned_sequence'] = df['precursor_sequence'].str.replace('n[42]', '')
    df['cleaned_sequence'] = df['cleaned_sequence'].str.replace('N[.98]', 'N').str.replace('Q[.98]', 'Q').str.replace(
        'M[15.99]', 'M').str.replace('C[57.02]', 'C')
    df['sequence_len'] = df['cleaned_sequence'].apply(len)

    df['precursor_id'] = df['charge'].astype(str) + '_' + df['precursor_sequence']
    df['psm_id'] = df['scan_number'].astype(str) + '_' + df['precursor_id']

    df = df[(df['sequence_len'] <= 50) & (df['sequence_len'] >= 7)]

    # charge in [2,5]
    df = df[(df['charge'] <= 5) & (df['charge'] >= 2)]

    df = df.reset_index()

    return df


def sage_first_clean_one_file(base_file_name, sage_result_file_path, sage_first_clean_dir, is_wiff, scan_id_dict):
    if is_wiff:
        data = clean_sage_frist_wiff(sage_result_file_path, scan_id_dict)
    else:
        data = clean_sage_frist(sage_result_file_path)
    data.to_csv(os.path.join(sage_first_clean_dir, f'{base_file_name}.csv'), index=False)
