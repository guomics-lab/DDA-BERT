import os.path

import matplotlib
import matplotlib.pyplot as plt
from peptdeep.model.rt import *
from peptdeep.pretrained_models import ModelManager

matplotlib.use('pdf')


def deal_process_with_exception(sage_result_file_path, alphapept_csv_path, output_dir='output'):
    try:
        print(sage_result_file_path, alphapept_csv_path)
        deal_process(sage_result_file_path, alphapept_csv_path, output_dir)
    except Exception as e:
        print(e)


def draw_pic(x_list, y_list, out_file_dir, csv_name):
    plt.figure(figsize=(6, 6))
    plt.scatter(x_list, y_list)
    # plt.plot(line_x_list, line_y_list, c="red")
    plt.xlabel("pred irt")
    plt.ylabel("aligned_rt")
    plt.title("DDA-BERT RT normalization, {}".format(os.path.split(out_file_dir)[-1]))
    plt.savefig(os.path.join(out_file_dir, f"{csv_name}_fitting_tutorials.pdf"))
    plt.close()


def find_special_characters_positions(sequence, characters):
    positions = {char: [] for char in characters}

    for index, char in enumerate(sequence):
        if char in characters:
            positions[char].append(index)

    return positions


def convert_seq(seq):
    seq_new = (str(seq).replace('C[+57.0216]', '-').replace('M[+15.9949]', '=').
               replace('Q[+0.98]', '+').
               replace('[+42]-', '*').
               replace('N[+0.98]', '.'))
    characters = ['-', '=', '+', '*', '.']
    mods = []
    mod_sites = []
    positions = find_special_characters_positions(seq_new, characters)
    pos_list1 = positions.get('-')
    pos_list2 = positions.get('=')
    pos_list3 = positions.get('+')
    pos_list4 = positions.get('*')
    pos_list5 = positions.get('.')

    mods.extend(['Carbamidomethyl@C'] * len(pos_list1))
    mods.extend(['Oxidation@M'] * len(pos_list2))
    mods.extend(['Deamidated@Q'] * len(pos_list3))
    mods.extend(['Acetyl@Any_N-term'] * len(pos_list4))
    mods.extend(['Deamidated@N'] * len(pos_list5))

    mod_sites.extend([str(nn) for nn in pos_list1])
    mod_sites.extend([str(nn) for nn in pos_list2])

    seq_fin = seq_new.replace('-', 'C').replace('=', 'M').replace('+', 'Q').replace('*', '').replace('.', 'N')

    mods_str = ';'.join(mods)
    mod_sites_str = ';'.join(mod_sites)
    return seq_fin, mods_str, mod_sites_str, len(seq_fin)


def convert_seq_len(seq):
    seq_new = (str(seq).replace('C[+57.0216]', '-').replace('M[+15.9949]', '=').
               replace('Q[+0.98]', '+').
               replace('[+42]-', '*').
               replace('N[+0.98]', '.'))
    seq_fin = seq_new.replace('-', 'C').replace('=', 'M').replace('+', 'Q').replace('*', '').replace('.', 'N')
    return len(seq_fin)


def clear_pp(seq):
    return str(seq).removeprefix('b\'').removesuffix('\'').replace('_decoy', '')


def calc_label(seq):
    if '_decoy' in seq:
        return 0
    else:
        return 1


def only_convert_seq_frag(seq):
    seq_new = (str(seq).replace('_decoy', '').replace('cC', '-').replace('oxM', '=').
               replace('deamQ', '+').
               replace('a<^', '*').
               replace('deamN', '.'))
    return seq_new


def convert_seq_frag(seq):
    seq_new = only_convert_seq_frag(seq)
    characters = ['-', '=', '+', '*', '.']
    mods = []
    mod_sites = []
    positions = find_special_characters_positions(seq_new, characters)
    pos_list1 = positions.get('-')
    pos_list2 = positions.get('=')
    pos_list3 = positions.get('+')
    pos_list4 = positions.get('*')
    pos_list5 = positions.get('.')

    mods.extend(['Carbamidomethyl@C'] * len(pos_list1))
    mods.extend(['Oxidation@M'] * len(pos_list2))
    mods.extend(['Deamidated@Q'] * len(pos_list3))
    mods.extend(['Acetyl@Any_N-term'] * len(pos_list4))
    mods.extend(['Deamidated@N'] * len(pos_list5))

    mod_sites.extend([str(nn) for nn in pos_list1])
    mod_sites.extend([str(nn) for nn in pos_list2])

    seq_fin = seq_new.replace('-', 'C').replace('=', 'M').replace('+', 'Q').replace('*', '').replace('.', 'N')
    mods_str = ';'.join(mods)
    mod_sites_str = ';'.join(mod_sites)
    return seq_fin, mods_str, mod_sites_str, len(seq_fin)


def get_pred_model(df, output_dir, csv_name):
    df['precursor_id'] = df['peptide'].astype(str) + df['charge'].astype(str)
    df['sequence_len'] = df['peptide'].apply(lambda x: convert_seq_len(x))

    counts = df.groupby('precursor_id').size()
    filtered_precursor_ids = counts[counts >= 3].index
    filtered_df = df[~df['precursor_id'].isin(filtered_precursor_ids)]
    filtered_df['f1'] = filtered_df['matched_peaks'] / (filtered_df['sequence_len'] + 6)
    filtered_df = filtered_df[filtered_df['f1'] > 0.7]

    filtered_df.sort_values(by='peptide_q', ascending=True, inplace=True)
    filtered_df = filtered_df[filtered_df['predicted_rt'] != 1]
    filtered_df = filtered_df.head(10000)

    seq_list = filtered_df['peptide'].tolist()
    aligned_rt_list = filtered_df['aligned_rt'].tolist()

    ddd = [convert_seq(ss) for ss in seq_list]
    train_df_data = pd.DataFrame({
        'sequence': [nn[0] for nn in ddd],
        'mods': [nn[1] for nn in ddd],
        'mod_sites': [nn[2] for nn in ddd],
        'nAA': [nn[3] for nn in ddd],
        'rt_norm': aligned_rt_list
    })

    # | hide
    models = ModelManager(device='cpu')
    models.load_installed_models()
    models.rt_model.train(train_df_data)

    df_w_rt_prediction = models.rt_model.predict(train_df_data)
    df_w_irt_prediction_added = models.rt_model.add_irt_column_to_precursor_df(train_df_data)

    irt_pred_list = df_w_irt_prediction_added['irt_pred'].tolist()

    draw_pic(irt_pred_list, df_w_irt_prediction_added['rt_norm'].tolist(), output_dir, csv_name)

    return models


def deal_process(sage_csv_path, alphapept_csv_path, output_dir='output'):
    os.makedirs(output_dir, exist_ok=True)
    sage_csv_name = os.path.basename(sage_csv_path)
    base_csv_name = sage_csv_name.removesuffix('.sage.tsv')

    save_path = os.path.join(output_dir, base_csv_name + '.with_fitting.csv')
    if os.path.exists(save_path):
        return
    sage_df = pd.read_csv(sage_csv_path, sep='\t')
    fragpipe_df = pd.read_csv(alphapept_csv_path)
    fragpipe_df['label'] = fragpipe_df['sequence'].apply(calc_label)
    fragpipe_df['sequence'] = fragpipe_df['sequence'].apply(clear_pp)

    pred_model = get_pred_model(sage_df, output_dir, base_csv_name)

    fragpipe_df['sequenceNew'] = fragpipe_df['sequence'].apply(only_convert_seq_frag)

    fragpipe_df.to_csv(save_path + 'no_filter_.csv', index=False)

    fragpipe_df = fragpipe_df[~fragpipe_df['sequenceNew'].str.contains('a')]

    seq_list = fragpipe_df['sequence'].tolist()

    fragpipe_df['aligned_rt'] = fragpipe_df['rt'] / fragpipe_df['rt'].max()

    ddd = [convert_seq_frag(ss) for ss in seq_list]
    pred_df_data = pd.DataFrame({
        'sequence': [nn[0] for nn in ddd],
        'mods': [nn[1] for nn in ddd],
        'mod_sites': [nn[2] for nn in ddd],
        'nAA': [nn[3] for nn in ddd],
        'rt_norm': fragpipe_df['aligned_rt'].tolist()
    })

    df_w_rt_prediction = pred_model.rt_model.predict(pred_df_data)
    df_w_irt_prediction_added = pred_model.rt_model.add_irt_column_to_precursor_df(pred_df_data)

    fragpipe_df['predicted_rt_1'] = df_w_rt_prediction['rt_pred'].tolist()
    fragpipe_df['delta_rt_1'] = abs(fragpipe_df['aligned_rt'] - fragpipe_df['predicted_rt_1'])
    fragpipe_df['predicted_irt'] = df_w_irt_prediction_added['irt_pred'].tolist()

    fragpipe_df.to_csv(save_path, index=False)
