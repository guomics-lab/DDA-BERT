import os.path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

matplotlib.use('pdf')
from sklearn.preprocessing import KBinsDiscretizer
from ropwr import RobustPWRegression


def deal_process_with_exception(csv_path_list, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for csv_path in csv_path_list:
        try:
            csv_name = os.path.basename(csv_path)
            save_path = os.path.join(output_dir, csv_name + '.with_fitting.csv')
            deal_process(csv_path, save_path)
        except Exception as e:
            print(e)


def draw_pic(pw, splits, x_list, y_list, out_file_dir, csv_name):
    line_X = np.arange(min(x_list) - 0.1, max(x_list) + 0.1, 0.01)
    line_y = pw.predict(line_X)
    plt.figure(figsize=(6, 6))
    plt.scatter(x_list, y_list)
    for s in splits:
        plt.axvline(s, color="grey", linestyle="--")

    plt.plot(line_X, line_y, c="red")
    plt.xlabel("predicted_rt")
    plt.ylabel("aligned_rt")
    plt.title("DDA-BERT RT normalization, {}".format(os.path.split(out_file_dir)[-1]))
    plt.savefig(os.path.join(out_file_dir, f"{csv_name}_fitting_tutorials.pdf"))
    plt.close()


def deal_process(csv_path, save_path):
    # if os.path.exists(save_path):
    #     return

    df = pd.read_csv(csv_path)
    counts = df.groupby('precursor_id').size()
    filtered_precursor_ids = counts[counts >= 3].index
    filtered_df = df[~df['precursor_id'].isin(filtered_precursor_ids)]
    filtered_df.sort_values(by='peptide_q', ascending=True, inplace=True)
    filtered_df = filtered_df[filtered_df['predicted_rt'] != 1]
    # filtered_df = filtered_df[filtered_df['predicted_rt'] != 0]
    filtered_df = filtered_df.head(10000)

    filtered_df = filtered_df[filtered_df['peptide_q'] <= 0.1]

    aligned_predicted_rt_list = filtered_df[['aligned_rt', 'predicted_rt']]

    n_neighbors = 5
    aligned_predicted_rt_matrix_np = np.array(aligned_predicted_rt_list)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree').fit(aligned_predicted_rt_matrix_np)
    distances, indices = nbrs.kneighbors(aligned_predicted_rt_matrix_np)
    # print(distances)
    mean_distances = np.mean(distances, axis=1)
    distance_matrix = np.column_stack([indices[:, 0], mean_distances])

    distance_matrix = distance_matrix[np.argsort(distance_matrix[:, 1])]

    peak_num = 2200
    if len(distance_matrix) < peak_num:
        peak_num = int(len(distance_matrix) * 0.85)
    min_distance_matrix = distance_matrix[:peak_num]

    peak_data = aligned_predicted_rt_matrix_np[min_distance_matrix[:, 0].astype(int)]

    x_list = peak_data[:, 1]

    y_list = peak_data[:, 0]

    bins = np.linspace(min(x_list), max(x_list), 11)

    cut_df = pd.cut(x_list, bins=bins)

    cut_df_count_list = list(cut_df.value_counts())
    bins_list = [[bins[nn], bins[nn + 1], cut_df_count_list[nn]] for nn in range(len(bins) - 1)]
    bins_list.sort(key=lambda x: x[2])

    last_first_bin = bins_list[0]
    last_second_bin = bins_list[1]
    last_third_num = bins_list[2][2]

    last_first_add_list = []
    last_second_add_list = []

    no_use_distance_matrix = aligned_predicted_rt_matrix_np[distance_matrix[peak_num:, 0].astype(int)]

    for y_point, x_point in no_use_distance_matrix.tolist():
        if x_point > last_first_bin[0] and x_point <= last_first_bin[1] and len(last_first_add_list) < last_third_num - \
                last_first_bin[2]:
            last_first_add_list.append([y_point, x_point])
        elif x_point > last_second_bin[0] and x_point <= last_second_bin[1] and len(
                last_second_add_list) < last_third_num - last_second_bin[2]:
            last_second_add_list.append([y_point, x_point])

    vstack_list = []
    vstack_list.append(peak_data)
    if len(last_first_add_list) > 0:
        vstack_list.append(np.array(last_first_add_list))
    if len(last_second_add_list) > 0:
        vstack_list.append(np.array(last_second_add_list))
    vstack_list.append(np.array([[0, 0]] * 10))
    vstack_list.append(np.array([[1, 1]] * 3))
    peak_data = np.vstack(vstack_list)

    pw_x_list = peak_data[:, 1]
    pw_y_list = peak_data[:, 0]

    pw, splits = get_tutorials_param(pw_x_list, pw_y_list)

    predicted_rt_1_list = pw.predict(df['predicted_rt'])
    df['predicted_rt_1'] = predicted_rt_1_list
    df['delta_rt_1'] = abs(df['aligned_rt'] - df['predicted_rt_1'])
    df.to_csv(save_path, index=False)


def get_tutorials_param(x_list, y_list):
    strategy = 'uniform'
    n_bins = 10
    pw = RobustPWRegression(objective="huber", degree=1, continuous_deriv=False,
                            monotonic_trend="ascending", reg_l1=0, reg_l2=0, h_epsilon=1)

    x_list = np.array(x_list)
    y_list = np.array(y_list)

    est = KBinsDiscretizer(n_bins=n_bins, strategy=strategy)
    est.fit(x_list.reshape(-1, 1), y_list)
    splits = est.bin_edges_[0][1:-1]
    pw.fit(x_list, y_list, splits=splits)
    return pw, splits
