# -*- coding: utf-8 -*-
"""
跨引擎( fp / ap / sage )按 scan 竞争 -> PSM 级打分
修改：Peptide 层基于 PSM@1%FDR 去修饰后得到（不再独立计算 FDR）
"""

import os
import re
import warnings
from typing import Optional

import numpy as np
import pandas as pd
import polars as pl

warnings.filterwarnings("ignore")


def clean_mods(seq: str) -> str:
    if not isinstance(seq, str):
        return ""
    s = seq
    s = re.sub(r"[A-Z]\[\s*[-+0-9.]+\s*\]", lambda m: m.group(0)[0], s)
    s = re.sub(r"n\[\s*[-+0-9.]+\s*\]", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\[.*?\]", "", s)
    return s


def safe_read_ipc(path: str, logger) -> Optional[pl.DataFrame]:
    try:
        dl = pl.read_ipc(path)
        if "index" not in dl.columns:
            dl = dl.with_columns(pl.arange(0, pl.count()).alias("index"))
        return dl
    except Exception as e:
        logger.error(f"✗ IPC read failed: {path}, error={e}")
        return None


def safe_read_pred(pred_path: str, logger) -> Optional[pd.DataFrame]:
    if not pred_path or not os.path.exists(pred_path):
        return None
    try:
        df = pd.read_csv(pred_path, header=None)
        if df.shape[1] >= 2:
            df = df.iloc[:, :2]
            df.columns = ["index", "score"]
        else:
            df = pd.read_csv(pred_path, header="infer")
            if {"index", "score"}.issubset(df.columns):
                df = df[["index", "score"]]
            elif df.shape[1] >= 2:
                df = df.iloc[:, :2]
                df.columns = ["index", "score"]
            else:
                return None
        df = df.drop_duplicates(subset=["index"], keep="first")
        df["index"] = df["index"].astype(int)
        return df
    except Exception as e:
        logger.error(f"✗ Pred read failed: {pred_path}, error={e}")
        return None


def get_psm_q_values(df):
    df = df.sort_values(by='score', ascending=False, ignore_index=True)
    df['decoy'] = np.where(df['label'] == 1, 0, 1)
    target_num = (df.decoy == 0).cumsum()
    decoy_num = (df.decoy == 1).cumsum()
    target_num[target_num == 0] = 1
    decoy_num[decoy_num == 0] = 1
    df['q_value'] = decoy_num / target_num
    df['q_value'] = df['q_value'][::-1].cummin()
    return df


def process_dl(dl, pred_df, need_cols):
    pdf = dl.select(need_cols).to_pandas()
    pdf["index"] = pdf["index"].astype(int)

    matched = pdf["index"].isin(pred_df["index"])
    match_rate = matched.sum() / len(pdf) if len(pdf) > 0 else 0

    pdf["score"] = pdf["index"].map(pred_df.set_index("index")["score"])
    pdf = pdf[~pdf["score"].isna()].copy()
    return pdf


def load_all_from_dirs(
        sage_ipc_path,
        fp_ipc_path,
        ap_ipc_path,
        sage_pred_result_path,
        fp_pred_result_path,
        ap_pred_result_path,
        open_sage, open_fp, open_ap, logger
):
    dfs = []
    need_cols = ["scan_number", "precursor_charge", "modified_sequence", "label", "index", 'precursor_mz']

    if open_sage:
        sage_dl = safe_read_ipc(sage_ipc_path, logger)
        sage_pred_df = safe_read_pred(sage_pred_result_path, logger)
        sage_pdf = process_dl(sage_dl, sage_pred_df, need_cols)
        dfs.append(sage_pdf)

    if open_fp:
        fp_dl = safe_read_ipc(fp_ipc_path, logger)
        fp_pred_df = safe_read_pred(fp_pred_result_path, logger)
        fp_pdf = process_dl(fp_dl, fp_pred_df, need_cols)
        dfs.append(fp_pdf)

    if open_ap:
        ap_dl = safe_read_ipc(ap_ipc_path, logger)
        ap_pred_df = safe_read_pred(ap_pred_result_path, logger)
        ap_pdf = process_dl(ap_dl, ap_pred_df, need_cols)
        dfs.append(ap_pdf)

    return dfs


def compete_across_engines_keep_best(df_all: pd.DataFrame) -> pd.DataFrame:
    df_sorted = df_all.sort_values(["scan_number", "score"], ascending=[True, False], ignore_index=True)
    return df_sorted.drop_duplicates(subset=["scan_number"], keep="first")


def evaluate_one_base(base_file_name, dfs, out_dir, logger):
    all_df = pd.concat(dfs, ignore_index=True)
    if all_df.empty:
        return

    best_df = compete_across_engines_keep_best(all_df)
    best_df['psm_id'] = best_df['scan_number'].astype(str) + '_' + best_df['precursor_charge'].astype(str) + '_' + \
                        best_df['modified_sequence']
    if best_df.empty:
        return

    psm_sorted = best_df.sort_values("score", ascending=False, ignore_index=True)
    psm_q = get_psm_q_values(psm_sorted)
    psm_1pct = psm_q[psm_q["q_value"] <= 0.01].copy()

    psm_1pct['precursor_id'] = psm_1pct['precursor_charge'].astype(str) + '_' + psm_1pct['modified_sequence']
    psm_1pct = psm_1pct.sort_values(by='score', ascending=False, ignore_index=True)
    precursor = psm_1pct.drop_duplicates(subset=['precursor_id'], keep='first')

    before_rows = len(psm_1pct)
    before_unique_precursors = psm_1pct["precursor_id"].nunique(dropna=False)
    after_rows = len(precursor)
    after_unique_precursors = precursor["precursor_id"].nunique(dropna=False)

    logger.info(f"[precursor dedup] rows: {before_rows} -> {after_rows} (removed {before_rows - after_rows})")
    logger.info(f"[precursor dedup] unique precursor_id: {before_unique_precursors} -> {after_unique_precursors}")

    precursor["cleaned_sequence"] = precursor["modified_sequence"].map(clean_mods)

    pep_1pct = (
        precursor.sort_values("score", ascending=False, ignore_index=True)
        .drop_duplicates(subset=["cleaned_sequence"], keep="first")
        .reset_index(drop=True)
    )
    # ==================================================

    psm_1pct.to_csv(os.path.join(out_dir, "psms.csv"), index=False)
    pep_1pct.to_csv(os.path.join(out_dir, "fdr_peptide.csv"), index=False)

    logger.info(
        f"[{base_file_name}]: 1%FDR_PSM={len(psm_1pct)}, 1%FDR_Peptide={len(pep_1pct)}"
    )

    protein_infer_input = psm_1pct[
        ['psm_id', 'precursor_mz', 'precursor_charge', 'modified_sequence',
         'label', 'score', 'q_value', 'scan_number', 'precursor_id']]
    protein_infer_input.to_csv(os.path.join(out_dir, 'protein_infer_input.csv'), index=False)

    target_fdr_psm = psm_1pct[psm_1pct['label'] == 1]
    save_target_fdr_psm = target_fdr_psm.drop(['mz_array', 'intensity_array'], errors='ignore')
    save_target_fdr_psm.to_csv(os.path.join(out_dir, 'fdr_psms.csv'), index=False)


def deal_process(base_file_name, sage_ipc, ap_ipc, fp_ipc, sage_pred, ap_pred, fp_pred, output_dir, open_sage, open_fp,
                 open_ap, logger):
    dfs = load_all_from_dirs(sage_ipc, fp_ipc, ap_ipc, sage_pred, fp_pred, ap_pred, open_sage, open_fp, open_ap, logger)
    evaluate_one_base(base_file_name, dfs, output_dir, logger)
