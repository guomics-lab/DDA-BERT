import numpy as np
import pandas as pd

residues_frag = {
    'M[147]': 'oxM',
    'M[15.9949]': 'oxM',
    'C[57.0215]': 'cC',
    'N[0.9800]': 'deamN',
    'Q[0.9800]': 'deamQ',
    'n[42.0106]': 'a',
    'n[420106]': 'a'
}


def clean_psm_func(peptide, residues_dict):
    if 'n[42.0106]' in residues_dict:
        peptide = peptide.replace('n[42.0106]', residues_dict['n[42.0106]'])
    for key, value in residues_dict.items():
        if value not in peptide:
            peptide = peptide.replace(key, value)
    return peptide


def clean_fp_frist(file_path):
    df = pd.read_table(file_path)

    df['Peptide'] = df['Peptide'].apply(lambda x: x[2:-2]).astype(str)

    df['precursor_sequence'] = df['Peptide'].apply(lambda x: clean_psm_func(x, residues_frag))

    df['precursor_sequence'] = df['precursor_sequence'].str.replace(r'\d+$', '', regex=True)

    df.loc[df['Label'] == 1, 'Label'] = 1
    df.loc[df['Label'] == -1, 'Label'] = 0
    df.rename(columns={'Label': 'label'}, inplace=True)
    df['label'] = df['label'].astype(int)

    df['charge'] = df['SpecId'].astype(str).apply(lambda x: x.split('.')[-1].split('_')[0]).astype(int)

    df['scan_number'] = df['ScanNr'].astype(np.int32)

    if 'pred_RT_real_units' in df.columns and 'retentiontime' in df.columns:
        df = df.rename(columns={
            'pred_RT_real_units': 'predicted_rt',
            'retentiontime': 'rt'
        })

    df.rename(columns={'Label': 'label'}, inplace=True)
    df['label'] = df['label'].astype(int)

    df['cleaned_sequence'] = df['precursor_sequence']

    df['sequence'] = df.apply(
        lambda row: row['precursor_sequence'] + '_decoy' if row['label'] == 0 else row['precursor_sequence'],
        axis=1
    )
    df['rt'] = df['retentiontime']

    df['cleaned_sequence'] = df['cleaned_sequence'].str.replace('oxM', 'M').str.replace('cC', 'C').str.replace('deamN',
                                                                                                               'N').str.replace(
        'deamQ', 'Q').str.replace('a', '')

    df['sequence_len'] = df['cleaned_sequence'].apply(len)

    df['precursor_id'] = df['charge'].astype(str) + '_' + df['precursor_sequence']
    df['psm_id'] = df['scan_number'].astype(str) + '_' + df['precursor_id']

    df = df[(df['sequence_len'] <= 50) & (df['sequence_len'] >= 7)]
    df = df[(df['charge'] <= 5) & (df['charge'] >= 2)]

    df = df.reset_index()
    df = df[[col for col in df.columns if not col.startswith(('charge_', 'length_', 'group_'))]]

    return df


def fp_frist_clean(pin_file, save_file_path):
    data = clean_fp_frist(pin_file)

    data.to_csv(save_file_path, index=False)
