import warnings
from functools import reduce

import polars as pl
from joblib import load

warnings.filterwarnings("ignore")

import os


def merge(base_file_name, protein_infer_input_csv, sage_path, fp_clean_file, alphepept_pept_dict_path, out_path):

    alphepept_pept_dict = pl.DataFrame(list(load(alphepept_pept_dict_path).items()),
                                       schema=["cleaned_sequence", "pro"]).with_columns(
        pl.col("pro").list.join(";").alias("Proteins")
    )

    inf_input = pl.read_csv(protein_infer_input_csv)
    inf_input = inf_input.with_columns(
        pl.col("modified_sequence").str.replace_all(r"\([^()]*\)|\[.*?\]|[^A-Z]", '').alias('cleaned_sequence'),
        pl.col("precursor_charge")
        .alias("charge"),
    )

    fp_clean = pl.read_csv(fp_clean_file).with_columns(
        pl.col("precursor_sequence").alias("modified_sequence"))

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
        pl.col("psm_id")
    )

    fp_clean = fp_clean.with_columns(
        replacement_expr.alias("psm_id")
    )

    sage_clean = pl.read_csv(sage_path,
                             separator='\t')
    sage_clean = sage_clean.with_columns(
        pl.col("peptide").str.replace_all(r"\([^()]*\)|\[.*?\]|[^A-Z]", '').alias('cleaned_sequence'),
        pl.col("peptide")
        .str.replace_all('C[+57.0216]', 'C[57.02]', literal=True)
        .str.replace_all('M[+15.9949]', 'M[15.99]', literal=True)
        .str.replace_all('Q[+0.98]', 'Q[.98]', literal=True)
        .str.replace_all('[+42]-', 'n[42]', literal=True)
        .str.replace_all('N[+0.98]', 'N[.98]', literal=True)
        .alias("modified_sequence"),
        pl.col("scannr")
        .str.split("=")
        .list.get(-1)
        .alias("scan_number"),
        pl.col("proteins").alias("Proteins")
    ).with_columns(
        pl.concat_str(
            [pl.col("scan_number").cast(pl.Utf8),
             pl.col("charge").cast(pl.Utf8),
             pl.col("modified_sequence")],
            separator="_"
        ).alias("psm_id")
    )

    merged = inf_input.join(fp_clean.select(['psm_id', 'Proteins']), on='psm_id', how='left')
    # cols = ['cleaned_sequence', 'Proteins', 'label', 'score', 'precursor_id', 'precursor_mz', 'precursor_charge', 'modified_sequence', 'q_value']
    cols = ['cleaned_sequence', 'Proteins', 'label', 'score', 'precursor_id', 'precursor_charge', 'modified_sequence',
            'q_value']
    unmatched = merged.filter(pl.col("Proteins").is_null())
    matched = merged.filter(pl.col("Proteins").is_not_null())

    unmatched = unmatched.join(sage_clean.select(['psm_id', 'Proteins']), on='psm_id', how='left').with_columns(
        pl.when(pl.col("Proteins_right").is_not_null())
        .then(pl.concat_str([pl.col("Proteins"), pl.col("Proteins_right")], separator=';'))
        .otherwise(pl.col("Proteins"))
        .alias("Proteins")
    ).drop(["Proteins_right"])
    matched = pl.concat([matched.select(cols), unmatched.filter(pl.col('Proteins').is_not_null()).select(cols)])
    unmatched = unmatched.filter(pl.col("Proteins").is_null())

    unmatched = unmatched.join(
        matched.select("cleaned_sequence"),
        on="cleaned_sequence",
        how="anti"
    )
    unmatched = unmatched.join(alphepept_pept_dict.select(['Proteins', 'cleaned_sequence']), on='cleaned_sequence',
                               how='left').with_columns(
        pl.when(pl.col("Proteins_right").is_not_null())
        .then(pl.concat_str([pl.col("Proteins"), pl.col("Proteins_right")], separator=';'))
        .otherwise(pl.col("Proteins"))
        .alias("Proteins")
    ).drop(["Proteins_right"])

    matched = pl.concat([matched.select(cols), unmatched.filter(pl.col('Proteins').is_not_null()).select(cols)])
    matched = matched.rename({'cleaned_sequence': 'sequence', 'Proteins': 'protein'})
    # print(unmatched['Proteins'].null_count() / len(merged))
    matched.write_csv(os.path.join(out_path, f'{base_file_name}_protein.csv'))
    return matched
