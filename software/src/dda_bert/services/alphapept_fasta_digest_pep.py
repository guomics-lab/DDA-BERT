import functools
import os
import re
from multiprocessing import Pool

import tqdm
from alphapept import constants
from alphapept.fasta import generate_peptides, add_to_pept_dict, generate_fasta_list, block_idx, merge_pept_dicts
from alphapept.performance import set_worker_count
from alphapept.settings import load_settings
from joblib import dump

mass_dict = constants.mass_dict


def keep_uppercase_letters(s):
    return re.sub(r'[^A-Z]', '', s)


def clean_peptide(peptide: str) -> str:
    # tag decoy
    is_decoy = peptide.endswith("_decoy")

    # replace decoy
    peptide = peptide.replace("_decoy", "")

    peptide = keep_uppercase_letters(peptide)
    return peptide


def tqdm_wrapper(pbar, update: float) -> None:
    """Update a qdm progress bar.

    Args:
        pbar (type): a tqd,.tqdm objet.
        update (float): The new value for the progressbar.

    """
    current_value = pbar.n
    delta = update - current_value
    pbar.update(delta)


# This function is a wrapper function and to be tested by the integration test
def digest_fasta_block(to_process: tuple) -> (list, dict):
    """
    Digest and create spectra for a whole fasta_block for multiprocessing. See generate_database_parallel.
    """
    fasta_index, fasta_block, settings = to_process

    pept_dict = {}
    for element in fasta_block:
        sequence = element["sequence"]
        try:
            mod_peptides = generate_peptides(sequence, **settings['fasta'])

            mod_peptides = list(set([clean_peptide(pep) for pep in mod_peptides]))

            pept_dict, _ = add_to_pept_dict(pept_dict, mod_peptides, element["name"])
        except:
            continue

    return pept_dict


# This function is a wrapper function and to be tested by the integration test
def generate_database_parallel(settings: dict, logger, callback=None):
    """
    Function to generate a database from a fasta file in parallel.
    Args:
        settings: alphapept settings.
    Returns:
        list: theoretical spectra. See generate_spectra()
        dict: peptide dict. See add_to_pept_dict()
        dict: fasta_dict. See generate_fasta_list()
    """

    n_processes = set_worker_count(
        worker_count=settings['general']['n_processes'],
        set_global=False
    )

    fasta_list, fasta_dict = generate_fasta_list(fasta_paths=settings['experiment']['fasta_paths'], **settings['fasta'])
    logger.info(f'FASTA contains {len(fasta_list):,} entries.')

    blocks = block_idx(len(fasta_list), settings['fasta']['fasta_block'])

    to_process = [(idx_start, fasta_list[idx_start:idx_end], settings) for idx_start, idx_end in blocks]

    pept_dicts = []
    with Pool(n_processes) as p:
        max_ = len(to_process)
        for i, ret in enumerate(p.imap_unordered(digest_fasta_block, to_process)):
            if callback:
                callback((i + 1) / max_)
            pept_dicts.append(ret)

    pept_dict = merge_pept_dicts(pept_dicts)
    return pept_dict


def deal(alphapept_config_path, logger):
    settings = load_settings(alphapept_config_path)
    fasta_path = settings['experiment']['fasta_paths'][0]
    fasta_name = os.path.basename(fasta_path).split('.')[0]

    output_dir = os.path.dirname(fasta_path)

    alphepept_pept_dict_path = os.path.join(output_dir, f'{fasta_name}_pept_dict.pkl')
    if os.path.exists(alphepept_pept_dict_path):
        return alphepept_pept_dict_path
    pept_dict = generate_database_parallel(settings, logger,
                                           callback=functools.partial(tqdm_wrapper, tqdm.tqdm(total=1)))

    dump(pept_dict, alphepept_pept_dict_path)

    return alphepept_pept_dict_path
