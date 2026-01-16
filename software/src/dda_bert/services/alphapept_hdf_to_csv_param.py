import h5py
import pandas as pd


def hdf_to_csv(file_path, alphapept_csv_path, min_frag_hits):
    file = h5py.File(file_path, 'r')
    first_search_col = ['charge', 'mz', 'rt', 'scan_no', 'sequence', 'n_fragments_matched']
    data = {}

    for key in file['first_search'].keys():
        data[key] = file['first_search'][key][:]

    print(f"Columns in 'first_search' for file {file_path}: {list(data.keys())}")
    first_search = pd.DataFrame(data)

    first_search = first_search[first_search['n_fragments_matched'] >= min_frag_hits]
    available_cols = [col for col in first_search_col if col in first_search.columns]
    first_search = first_search[available_cols]

    first_search.to_csv(alphapept_csv_path, index=False)
