import yaml

from ..obj.Entity import InputParam


def read_config_file(config_file_path):
    with open(config_file_path) as f:
        config_data = yaml.load(f.read(), Loader=yaml.FullLoader)
        input_param = InputParam()
        input_param.env = 'win'

        mzml_list = config_data.get('selectedFiles')
        if len(mzml_list) == 0:
            raise Exception('mzML path must get')

        input_param.mzml_file_list = mzml_list

        fasta = config_data.get('input').get('fastaPath')
        input_param.fasta_path = fasta
        if not fasta:
            raise Exception('fasta path must get')

        output_path = config_data.get('input').get('outputDirPath')
        if not output_path:
            raise Exception('output path must get')
        input_param.output_dir_path = output_path

        input_param.fragment_min_mz = config_data.get('datasetParam').get('fragment_min_mz')
        input_param.fragment_max_mz = config_data.get('datasetParam').get('fragment_max_mz')
        input_param.enzyme_missed_cleavages = config_data.get('datasetParam').get('enzyme_missed_cleavages')
        input_param.digest_min_len = config_data.get('datasetParam').get('digest_min_len')
        input_param.digest_max_len = config_data.get('datasetParam').get('digest_max_len')
        input_param.enzyme_cleave_at = config_data.get('datasetParam').get('enzyme_cleave_at')
        input_param.enzyme_restrict = config_data.get('datasetParam').get('enzyme_restrict')
        input_param.max_variable_mods = config_data.get('datasetParam').get('max_variable_mods')
        input_param.decoy_tag = config_data.get('datasetParam').get('decoy_tag')
        if input_param.decoy_tag and len(input_param.decoy_tag) > 0:
            input_param.generate_decoys = True
        else:
            input_param.generate_decoys = False
        input_param.min_matched_peaks = config_data.get('datasetParam').get('min_matched_peaks')
        input_param.max_fragment_charge = config_data.get('datasetParam').get('max_fragment_charge')

        input_param.min_precursor_charge = config_data.get('datasetParam').get('min_precursor_charge')
        input_param.max_precursor_charge = config_data.get('datasetParam').get('max_precursor_charge')

        input_param.min_precursor_tol_ppm = config_data.get('datasetParam').get('min_precursor_tol_ppm')
        input_param.max_precursor_tol_ppm = config_data.get('datasetParam').get('max_precursor_tol_ppm')

        input_param.min_fragment_tol_ppm = config_data.get('datasetParam').get('min_fragment_tol_ppm')
        input_param.max_fragment_tol_ppm = config_data.get('datasetParam').get('max_fragment_tol_ppm')

        input_param.static_mods_c = config_data.get('staticMods').get('staticModsC')
        input_param.variable_mods_a = config_data.get('staticMods').get('staticModsA')
        input_param.variable_mods_m = config_data.get('staticMods').get('staticModsM')
        input_param.variable_mods_n = config_data.get('staticMods').get('staticModsN')
        input_param.variable_mods_q = config_data.get('staticMods').get('staticModsQ')

        return input_param
