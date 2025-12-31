import os


def read_txt(txt_path):
    mzml_list = []
    with open(txt_path, encoding='utf-8') as f:
        lines = f.readlines()
        for line in lines:
            mzml_list.append(line.strip())
    return mzml_list


def is_d(file_path):
    file_path = str(file_path)
    is_d = False
    if os.path.isdir(file_path) and str(file_path).endswith('.d'):
        is_d = True
    return is_d
