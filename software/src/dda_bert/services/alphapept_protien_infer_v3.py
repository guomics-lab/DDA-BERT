import warnings
from collections import defaultdict

import networkx as nx
from joblib import load
from numba import njit

warnings.filterwarnings("ignore")

import os
import numpy as np
import re
import pandas as pd


def assign_proteins(data: pd.DataFrame, pept_dict: dict) -> (pd.DataFrame, dict):
    """
    Assign psms to proteins.
    This function appends the dataframe with a column 'n_possible_proteins' which indicates how many proteins a psm could be matched to.
    It returns the appended dataframe and a dictionary `found_proteins` where each protein is mapped to the psms indices.

    Args:
        data (pd.DataFrame): psms table of scored and filtered search results from alphapept.
        pept_dict (dict): dictionary that matches peptide sequences to proteins

    Returns:
        pd.DataFrame: psms table of search results from alphapept appended with the number of matched proteins.
        dict: dictionary mapping psms indices to proteins.

    """

    data = data.reset_index(drop=True)

    data['n_possible_proteins'] = data['sequence'].apply(lambda x: len(pept_dict[x]))
    unique_peptides = (data['n_possible_proteins'] == 1).sum()
    shared_peptides = (data['n_possible_proteins'] > 1).sum()

    # logging.info(f'A total of {unique_peptides:,} unique and {shared_peptides:,} shared peptides.')

    sub = data[data['n_possible_proteins'] == 1]
    psms_to_protein = sub['sequence'].apply(lambda x: pept_dict[x])

    found_proteins = {}
    for idx, _ in enumerate(psms_to_protein):
        idx_ = psms_to_protein.index[idx]
        p_str = '@' + str(_[0])
        if p_str in found_proteins:
            found_proteins[p_str] = found_proteins[p_str] + [str(idx_)]
        else:
            found_proteins[p_str] = [str(idx_)]

    return data, found_proteins


def get_shared_proteins(data: pd.DataFrame, found_proteins: dict, pept_dict: dict) -> dict:
    """
    Assign peptides to razor proteins.

    Args:
        data (pd.DataFrame): psms table of scored and filtered search results from alphapept, appended with `n_possible_proteins`.
        found_proteins (dict): dictionary mapping psms indices to proteins
        pept_dict (dict): dictionary mapping peptide indices to the originating proteins as a list

    Returns:
        dict: dictionary mapping peptides to razor proteins

    """

    G = nx.Graph()

    sub = data[data['n_possible_proteins'] > 1]

    for i in range(len(sub)):
        seq, score = sub.iloc[i][['sequence', 'score']]
        idx = sub.index[i]
        possible_proteins = pept_dict[seq]

        for p in possible_proteins:
            G.add_edge(str(idx), '@' + str(p), score=score)

    connected_groups = np.array([list(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)],
                                dtype=object)
    n_groups = len(connected_groups)

    # logging.info('A total of {} ambigious proteins'.format(len(connected_groups)))

    # Solving with razor:
    found_proteins_razor = {}
    for a in connected_groups[::-1]:
        H = G.subgraph(a).copy()
        shared_proteins = list(np.array(a)[np.array(list(i[0] == '@' for i in a))])

        while len(shared_proteins) > 0:
            neighbors_list = []

            for node in shared_proteins:
                shared_peptides = list(H.neighbors(node))

                if node in G:
                    if node in found_proteins.keys():
                        shared_peptides += found_proteins[node]

                n_neigbhors = len(shared_peptides)

                neighbors_list.append((n_neigbhors, node, shared_peptides))

            # Check if we have a protein_group (e.g. they share the same everythin)
            neighbors_list.sort()

            # Check for protein group
            node_ = [neighbors_list[-1][1]]
            idx = 1
            while idx < len(neighbors_list):  # Check for protein groups
                if neighbors_list[-idx][0] == neighbors_list[-idx - 1][0]:  # lenght check
                    if set(neighbors_list[-idx][2]) == set(neighbors_list[-idx - 1][2]):  # identical peptides
                        node_.append(neighbors_list[-idx - 1][1])
                        idx += 1
                    else:
                        break
                else:
                    break

            # Remove the last entry:
            shared_peptides = neighbors_list[-1][2]
            for node in node_:
                shared_proteins.remove(node)

            for _ in shared_peptides:
                if _ in H:
                    H.remove_node(_)

            if len(shared_peptides) > 0:
                if len(node_) > 1:
                    node_ = tuple(node_)
                else:
                    node_ = node_[0]

                found_proteins_razor[node_] = shared_peptides

    return found_proteins_razor


def get_protein_groups(data: pd.DataFrame, pept_dict: dict, decoy=False, decoy_symbol='REV_', **kwargs) -> pd.DataFrame:
    """
    Function to perform protein grouping by razor approach.
    This function calls `assign_proteins` and `get_shared_proteins`.
    Each protein is indicated with a p -> protein index

    Args:
        data (pd.DataFrame): psms table of scored and filtered search results from alphapept.
        pept_dict (dict): A dictionary mapping peptide indices to the originating proteins as a list.
        decoy (bool, optional): Defaults to False.
    Returns:
        pd.DataFrame: alphapept results table now including protein level information.
    """
    data, found_proteins = assign_proteins(data, pept_dict)
    found_proteins_razor = get_shared_proteins(data, found_proteins, pept_dict)

    report = data.copy()

    assignment = np.zeros(len(report), dtype=object)
    assignment[:] = ''
    assignment_pg = assignment.copy()

    assignment_idx = assignment.copy()
    assignment_idx[:] = ''

    razor = assignment.copy()
    razor[:] = False

    if decoy:
        add = decoy_symbol
    else:
        add = ''

    for protein_str in found_proteins.keys():
        protein_name = add + protein_str
        indexes = [int(_) for _ in found_proteins[protein_str]]
        assignment[indexes] = protein_name
        assignment_pg[indexes] = protein_name
        assignment_idx[indexes] = protein_name

    for protein_str in found_proteins_razor.keys():
        indexes = [int(_) for _ in found_proteins_razor[protein_str]]

        if isinstance(protein_str, tuple):
            protein_name = ','.join([add + _ for _ in protein_str])
            protein = ','.join([str(_) for _ in protein_str])

        else:
            protein = protein_str
            protein_name = add + protein_str

        assignment[indexes] = protein_name
        assignment_pg[indexes] = protein_name
        assignment_idx[indexes] = str(protein)
        razor[indexes] = True

    report['protein'] = assignment
    report['protein_group'] = assignment_pg
    report['razor'] = razor
    report['protein_idx'] = assignment_idx

    for col in ['protein', 'protein_group', 'protein_idx']:
        report[col] = report[col].str.replace("@", "")

    return report


def perform_protein_grouping(data: pd.DataFrame, pept_dict: dict, **kwargs) -> pd.DataFrame:
    """
    Wrapper function to perform protein grouping by razor approach

    Args:
        data (pd.DataFrame): psms table of scored and filtered search results from alphapept.
        pept_dict (dict): A dictionary mapping peptide indices to the originating proteins as a list.
        fasta_dict (dict): A dictionary with fasta sequences.

    Returns:
        pd.DataFrame: alphapept results table now including protein level information.
    """
    data_sub = data[['sequence', 'score', 'decoy']]
    data_sub_unique = data_sub.groupby(['sequence', 'decoy'], as_index=False).agg({"score": "max"})

    targets = data_sub_unique[data_sub_unique.decoy == False]
    targets = targets.reset_index(drop=True)
    protein_targets = get_protein_groups(targets, pept_dict, **kwargs)

    protein_targets['decoy_protein'] = False

    decoys = data_sub_unique[data_sub_unique.decoy == True]
    decoys = decoys.reset_index(drop=True)
    protein_decoys = get_protein_groups(decoys, pept_dict, decoy=True, **kwargs)

    protein_decoys['decoy_protein'] = True

    protein_groups = pd.concat([protein_targets, protein_decoys])
    protein_groups_app = protein_groups[
        ['sequence', 'decoy', 'protein', 'protein_group', 'razor', 'protein_idx', 'decoy_protein',
         'n_possible_proteins']]
    protein_report = pd.merge(data,
                              protein_groups_app,
                              how='inner',
                              on=['sequence', 'decoy'],
                              validate="many_to_one")

    return protein_report


@njit
def get_q_values(fdr_values: np.ndarray) -> np.ndarray:
    """
    Calculate q-values from fdr_values.

    Args:
        fdr_values (np.ndarray): np.ndarray of fdr values.

    Returns:
        np.ndarray: np.ndarray of q-values.
    """
    q_values = np.zeros_like(fdr_values)
    min_q_value = np.max(fdr_values)
    for i in range(len(fdr_values) - 1, -1, -1):
        fdr = fdr_values[i]
        if fdr < min_q_value:
            min_q_value = fdr
        q_values[i] = min_q_value

    return q_values


# Note that the test function for cut_fdr is further down in the notebook to also test protein-level FDR.
def cut_fdr(df: pd.DataFrame, fdr_level: float = 0.01, cut: bool = True) -> (float, pd.DataFrame):
    """
    Cuts a dataframe with a given fdr level

    Args:
        df (pd.DataFrame): psms table of search results from alphapept.
        fdr_level (float, optional): fdr level that should be used for filtering. The value should lie between 0 and 1. Defaults to 0.01.
        cut (bool, optional): flag to cut above fdr threshold. Defaults to 'True'.

    Returns:
        float: numerical value of the applied score cutoff
        pd.DataFrame: df with psms within fdr

    """

    df["target"] = ~df["decoy"]

    df = df.sort_values(by=["score", "decoy"], ascending=False)
    df = df.reset_index()

    df["target_cum"] = np.cumsum(df["target"])
    df["decoys_cum"] = np.cumsum(df["decoy"])

    df["fdr"] = df["decoys_cum"] / df["target_cum"]
    df["q_value"] = get_q_values(df["fdr"].values)

    last_q_value = df["q_value"].iloc[-1]
    first_q_value = df["q_value"].iloc[0]

    if last_q_value <= fdr_level:
        # logging.info('Last q_value {:.3f} of dataset is smaller than fdr_level {:.3f}'.format(last_q_value, fdr_level))
        cutoff_index = len(df) - 1

    elif first_q_value >= fdr_level:
        # logging.info('First q_value {:.3f} of dataset is larger than fdr_level {:.3f}'.format(last_q_value, fdr_level))
        cutoff_index = 0

    else:
        cutoff_index = df[df["q_value"].gt(fdr_level)].index[0] - 1

    cutoff_value = df.loc[cutoff_index]["score"]

    if cut:
        cutoff = df[df["score"] >= cutoff_value]
    else:
        cutoff = df

    targets = df.loc[cutoff_index, "target_cum"]
    decoy = df.loc[cutoff_index, "decoys_cum"]

    fdr = df.loc[cutoff_index, "fdr"]

    # logging.info(
    #     f"{targets:,} target ({decoy:,} decoy) of {len(df):,} PSMs. FDR {fdr:.6f} for a cutoff of {cutoff_value:.2f} (set FDR was {fdr_level}).")

    cutoff = cutoff.reset_index(drop=True)
    return cutoff_value, cutoff


def cut_global_fdr(data: pd.DataFrame, analyte_level: str = 'sequence', fdr_level: float = 0.01, strategy='max',
                   **kwargs) -> pd.DataFrame:
    """
    Function to estimate and filter by global peptide or protein fdr

    Args:
        data (pd.DataFrame): psms table of search results from alphapept.
        analyte_level (str, optional): string specifying the analyte level to apply the fdr threshold. Options include: 'precursor', 'sequence', 'protein_group' and 'protein'. Defaults to 'sequence'.
        fdr_level (float, optional): fdr level that should be used for filtering. The value should lie between 0 and 1. Defaults to 0.01.
        strategy(str): min, max, median, mean
    Returns:
        pd.DataFrame: df with filtered results

    """
    # logging.info('Global FDR on {}'.format(analyte_level))
    data_sub = data[[analyte_level, 'score', 'decoy']]
    data_sub_unique = data_sub.groupby([analyte_level, 'decoy'], as_index=False).agg({"score": strategy})

    analyte_levels = ['precursor', 'sequence', 'protein_group', 'protein']

    if analyte_level in analyte_levels:
        agg_score = data_sub_unique.groupby([analyte_level, 'decoy']).agg({"score": strategy}).reset_index()
    else:
        raise Exception('analyte_level should be either sequence or protein. The selected analyte_level was: {}'.format(
            analyte_level))

    agg_cval, agg_cutoff = cut_fdr(agg_score, fdr_level=fdr_level)

    # logging.info(f'Global FDR cutoff at {agg_cval:.3f}.')

    agg_report = data.reset_index(drop=True).merge(agg_cutoff,
                                                   how='inner',
                                                   on=[analyte_level, 'decoy'],
                                                   suffixes=('', '_' + analyte_level),
                                                   validate="many_to_one").set_index(
        'index')  # retain the original index
    agg_report = agg_report[agg_report['decoy'] == False]
    return agg_report


def remove_charge_modi(text):
    result = re.sub(r"\([^()]*\)|\[.*?\]|[^A-Z]", '', str(text))
    return result


def map_alphepept_pept_dict(pep, alphepept_pept_dict):
    value = alphepept_pept_dict.get(pep, '')

    if isinstance(value, str) and len(value) >= 2 and value.startswith('[') and value.endswith(']'):
        return value[1:-1]
    else:
        return value


def clean_protein_list(s):
    if not isinstance(s, str):
        return ''

    s_cleaned = re.sub(r"[\[\]']", "", s)
    proteins = []
    seen = set()
    for p in s_cleaned.split(';'):
        p = p.strip()
        if p and p not in seen:
            proteins.append(p)
            seen.add(p)
    return ';'.join(proteins)


def merge(base_file_name, protein_infer_input_csv, sage_path, fp_clean_file, alphepept_pept_dict_path, out_path,
          open_sage=True, open_fp=True, open_ap=True):
    DDABert_precursor = None
    if open_sage and open_fp and open_ap:
        DDABert_precursor = merge_all(protein_infer_input_csv, sage_path, fp_clean_file, alphepept_pept_dict_path)
    elif open_sage and open_fp:
        DDABert_precursor = merge_all(protein_infer_input_csv, sage_path, fp_clean_file, None)
    elif open_sage:
        DDABert_precursor = merge_all(protein_infer_input_csv, sage_path, None, None)

    process(base_file_name, DDABert_precursor, out_path)


def merge_all(protein_infer_input_csv, sage_path, fp_clean_file, alphepept_pept_dict_path):
    inf_input = pd.read_csv(protein_infer_input_csv)

    inf_input['sequence'] = inf_input['modified_sequence'].astype(str).apply(remove_charge_modi)
    inf_input['decoy'] = np.where(inf_input['label'] == 1, False, True)
    inf_input['psm_id'] = inf_input['scan_number'].astype(str) + '_' + inf_input[
        'precursor_charge'].astype(str) + '_' + inf_input['sequence'].astype(str)

    dd_list = []

    ap_DDABert_precursor = inf_input[inf_input['engine'] == 'ap']
    fp_DDABert_precursor = inf_input[inf_input['engine'] == 'fp']
    sage_DDABert_precursor = inf_input[inf_input['engine'] == 'sage']

    # from sage
    sage_precursor = pd.read_table(sage_path)
    sage_precursor['sequence'] = sage_precursor['peptide'].astype(str).apply(remove_charge_modi)
    sage_precursor["scan_number"] = sage_precursor["scannr"].astype(str).apply(lambda x: x.split('scan=')[-1])
    sage_precursor['psm_id'] = sage_precursor['scan_number'].astype(str) + '_' + sage_precursor['charge'].astype(
        str) + '_' + sage_precursor['sequence'].astype(str)
    sage_precursor = sage_precursor[['psm_id', 'proteins']].drop_duplicates(subset=['psm_id'])

    # 关联sage
    sage_DDABert_precursor['proteins'] = sage_DDABert_precursor['psm_id'].map(
        sage_precursor.set_index('psm_id')['proteins'])
    sage_DDABert_precursor['proteins'] = sage_DDABert_precursor['proteins'].astype(str).apply(
        lambda x: x.replace(',', ';')).apply(clean_protein_list)
    dd_list.append(sage_DDABert_precursor)

    if fp_clean_file:
        fp_precursor = pd.read_csv(fp_clean_file)
        fp_precursor['psm_id'] = fp_precursor['ScanNr'].astype(str) + '_' + fp_precursor['charge'].astype(str) + '_' + \
                                 fp_precursor['cleaned_sequence'].astype(str)
        fp_precursor = fp_precursor[['psm_id', 'Proteins']].drop_duplicates(subset=['psm_id'])

        fp_DDABert_precursor['proteins'] = fp_DDABert_precursor['psm_id'].map(
            fp_precursor.set_index('psm_id')['Proteins'])
        fp_DDABert_precursor['proteins'] = fp_DDABert_precursor['proteins'].astype(str).apply(
            lambda x: x.replace(',', ';')).apply(clean_protein_list)
        dd_list.append(fp_DDABert_precursor)

    if alphepept_pept_dict_path:

        alphepept_pept_dict = load(alphepept_pept_dict_path)
        # 关联ap
        ap_DDABert_precursor['proteins'] = ap_DDABert_precursor['sequence'].apply(
            lambda x: map_alphepept_pept_dict(x, alphepept_pept_dict))
        ap_DDABert_precursor['proteins'] = ap_DDABert_precursor['proteins'].astype(str).apply(
            lambda x: x.replace(',', ';')).apply(clean_protein_list)
        dd_list.append(ap_DDABert_precursor)

    # concat
    DDABert_precursor = pd.concat(dd_list)
    return DDABert_precursor


def process(base_file_name, DDABert_precursor, out_path):
    DDABert_precursor = DDABert_precursor[
        (DDABert_precursor['proteins'] != '') & (DDABert_precursor['proteins'].notnull())]

    # protein explode
    DDABert_precursor['proteins'] = DDABert_precursor['proteins'].astype(str).apply(
        lambda x: x.replace(',', ';').split(';'))
    DDABert_precursor_explode = DDABert_precursor.explode('proteins')
    DDABert_precursor_explode['proteins'] = DDABert_precursor_explode['proteins'].astype(str).replace('rev_', 'REV_')
    DDABert_precursor_explode = DDABert_precursor_explode[DDABert_precursor_explode['proteins'] != '']

    # peptide_2_protein
    peptide2protein = DDABert_precursor_explode[['sequence', 'proteins']].drop_duplicates(
        subset=['sequence', 'proteins'])
    peptide2protein = peptide2protein.groupby('sequence').apply(lambda x: x['proteins'].tolist()).reset_index()
    peptide2protein.columns = ['sequence', 'protein']
    peptide_dict = dict(
        [(sequence, protein) for sequence, protein in zip(peptide2protein["sequence"], peptide2protein["protein"])])

    default_pept_dict = defaultdict(lambda: '')
    default_pept_dict.update(peptide_dict)

    if 'fdr' in DDABert_precursor.columns:
        DDABert_precursor.drop(columns=['fdr'], inplace=True)

    DDABert_precursor = DDABert_precursor[DDABert_precursor['sequence'].isin(set(default_pept_dict.keys()))]

    protein_report = perform_protein_grouping(DDABert_precursor, default_pept_dict)
    protein_res = cut_global_fdr(protein_report, analyte_level="protein")

    save_cols = ['psm_id', 'sequence', 'protein', 'razor', 'decoy_protein', 'n_possible_proteins', 'score_protein',
                 'target', 'target_cum',
                 'decoys_cum', 'fdr', 'q_value_protein']
    protein_res = protein_res[save_cols]
    protein_res.drop_duplicates(inplace=True)

    protein_res = protein_res[
        (protein_res['protein'] != '') & (protein_res['protein'] != "['']") & (protein_res['protein'].notnull())]
    protein_res['protein'] = protein_res['protein'].astype(str)

    invalid = protein_res['protein'].isin(['nan', 'None', ''])

    starts_with_rev = protein_res['protein'].str.startswith(('rev_', 'REV_'))

    protein_res = protein_res[~(invalid | starts_with_rev)]

    protein_res.to_csv(os.path.join(out_path, f'{base_file_name}_protein.csv'), index=False)
