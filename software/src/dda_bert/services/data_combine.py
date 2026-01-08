import os

import numpy as np
import polars as pl


def get_fdr_result(df):
    df = df.to_pandas()
    df['decoy'] = np.where(df['label'] == 1, 0, 1)
    target_num = (df.decoy == 0).cumsum()
    decoy_num = (df.decoy == 1).cumsum()
    target_num[target_num == 0] = 1
    decoy_num[decoy_num == 0] = 1
    df['q_value'] = decoy_num / target_num
    df['q_value'] = df['q_value'][::-1].cummin()
    return pl.from_pandas(df)



    # ---------- helper that post-processes one DF ----------
def _post(df: pl.DataFrame, source_name: str) -> pl.DataFrame:
        df = (
            df.with_columns(
                pl.concat_str(
                    [pl.col("precursor_charge").cast(pl.Utf8),
                     pl.col("modified_sequence")],
                    separator="_"
                ).alias("precursor_id"),
                pl.concat_str(
                    [pl.col("scan_number").cast(pl.Utf8),
                     pl.col("precursor_charge").cast(pl.Utf8),
                     pl.col("modified_sequence")],
                    separator="_"
                ).alias("psm_id"),
                pl.lit(source_name).alias("source")
            )
            .sort("score", descending=True, maintain_order=True)
            .unique("scan_number", keep="first")
            .sort("score", descending=True, maintain_order=True)
        )
        df = get_fdr_result(df)
        df = df.filter(pl.col("q_value") < 0.1)
        # print(source_name, len(df))
        return df


def process_file(base_file_name, sage_ipc, ap_ipc, fp_ipc, sage_pred, ap_pred, fp_pred, output_dir, open_sage, open_fp, open_ap) -> pl.DataFrame:


    ipc_list = []

    if open_sage:
        sage_ipc = pl.read_ipc(sage_ipc)
        sage_score = pl.read_csv(sage_pred,
                                 has_header=False,
                                 new_columns=["index", "score", "label_y", "weight"]).unique('index')
        assert len(sage_ipc) == len(sage_score)
        sage_ipc = sage_ipc.join(sage_score.select(['index', 'score']), on="index", how="left")
        sage_ipc = _post(sage_ipc, "sage")
        ipc_list.append(sage_ipc)


    if open_fp:
        fp_ipc = pl.read_ipc(fp_ipc)
        fp_score = pl.read_csv(fp_pred,
                               has_header=False,
                               new_columns=["index", "score", "label_y", "weight"]).unique('index')

        assert len(fp_ipc) == len(fp_score)

        # ---------- join ----------

        fp_ipc = fp_ipc.join(fp_score.select(['index', 'score']), on="index", how="left").drop('aligned_rt')

        fp_ipc = _post(fp_ipc, "fp")
        ipc_list.append(fp_ipc)

        # common_cols = set.intersection(*[set(x.columns) for x in (sage_ipc, ap_ipc, fp_ipc)])
        # sage_ipc = sage_ipc.select(sorted(common_cols))
        # ap_ipc = ap_ipc.select(sorted(common_cols))
        # fp_ipc = fp_ipc.select(sorted(common_cols))
        # df = pl.concat([fp_ipc, sage_ipc, ap_ipc])

    if open_ap:
        ap_ipc = pl.read_ipc(ap_ipc)
        ap_score = pl.read_csv(ap_pred,
                               has_header=False,
                               new_columns=["index", "score", "label_y", "weight"]).unique('index')
        # ---------- sanity checks ----------
        assert len(ap_ipc) == len(ap_score)
        # ---------- join ----------
        ap_ipc = ap_ipc.join(ap_score.select(['index', 'score']), on="index", how="left")
        ap_ipc = _post(ap_ipc, "ap")
        ipc_list.append(ap_ipc)

        # common_cols = set.intersection(*[set(x.columns) for x in (sage_ipc, ap_ipc, fp_ipc)])
        # sage_ipc = sage_ipc.select(sorted(common_cols))
        # ap_ipc = ap_ipc.select(sorted(common_cols))
        # fp_ipc = fp_ipc.select(sorted(common_cols))
        # df = pl.concat([fp_ipc, sage_ipc, ap_ipc])
    #
    # else:
    #     common_cols = set.intersection(*[set(x.columns) for x in [sage_ipc]])
    #     sage_ipc = sage_ipc.select(sorted(common_cols))
    #     df = pl.concat([sage_ipc])

    common_cols = set.intersection(*[set(x.columns) for x in ipc_list])
    selected_ipc_list = []
    for ipc in ipc_list:
        selected_ipc_list.append(ipc.select(sorted(common_cols)))
    df = pl.concat(selected_ipc_list)

    # deduplicate at PSM level
    df = df.unique("scan_number", keep='first').sort("score", descending=True, maintain_order=True)

    df = get_fdr_result(df)

    df = df.filter(pl.col("q_value") < 0.01)
    target_fdr_psm = df.filter(pl.col('label') == 1)

    psm_count = len(target_fdr_psm)
    print("pre_filter PSM", psm_count)

    save_target_fdr_psm = target_fdr_psm.drop(['mz_array', 'intensity_array'], strict=False)
    # target_fdr_psm.write_ipc(os.path.join(output_dir, 'fdr_psms.ipc'))
    save_target_fdr_psm.write_csv(os.path.join(output_dir, 'fdr_psms.csv'))

    protein_infer_input = df[
        ['psm_id', 'precursor_mz', 'precursor_charge', 'modified_sequence',
         'label', 'score', 'q_value', 'scan_number', 'precursor_id']]
    protein_infer_input.write_csv(os.path.join(output_dir, 'protein_infer_input.csv'))

    # split targets / decoys, dedupe targets on precursor
    targets = df.filter(pl.col("label") == 1).unique("precursor_id", keep='first')
    decoys = df.filter(pl.col("label") == 0)
    df = pl.concat([targets, decoys]).sort("score", descending=True, maintain_order=True)

    df = df.with_columns(
        pl.col("modified_sequence")
        .str.replace_all("n[42]", "", literal=True)
        .str.replace_all("N[.98]", "N", literal=True)
        .str.replace_all("Q[.98]", "Q", literal=True)
        .str.replace_all("M[15.99]", "M", literal=True)
        .str.replace_all("C[57.02]", "C", literal=True)
        .alias("cleaned_sequence")
    )

    peptide_count = df.filter(pl.col('label') == 1)['cleaned_sequence'].n_unique()

    print("peptide", peptide_count)

    df = df.drop(['mz_array', 'intensity_array'], strict=False)
    df.write_csv(os.path.join(output_dir, 'psms.csv'))
    peptide_df = df.filter(pl.col('label') == 1).unique('cleaned_sequence')
    peptide_df.write_csv(os.path.join(output_dir, 'fdr_peptide.csv'))

    return psm_count, peptide_count


def deal_process(base_file_name, sage_ipc, ap_ipc, fp_ipc, sage_pred, ap_pred, fp_pred, output_dir, open_sage, open_fp, open_ap):
    psm_count, peptide_count = process_file(base_file_name, sage_ipc, ap_ipc, fp_ipc, sage_pred, ap_pred, fp_pred,
                                            output_dir, open_sage, open_fp, open_ap)
    df = pl.DataFrame({"base_file_name": [base_file_name], "psm": [psm_count], "peptide": [peptide_count]})
    df.write_csv(os.path.join(output_dir, "1fdr_stats.csv"))
