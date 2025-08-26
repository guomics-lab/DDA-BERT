import multiprocessing
from optparse import OptionParser

from src import common_config
from src.common import constant
from src.main_process_handler import MainProcessHandler
from src.obj.Entity import InputParam
from src.utils import file_utils, config_util
from src.utils.gpu_utils import get_usage_device, get_top_free_max_device

config_data = common_config.read_yml()

if __name__ == '__main__':
    multiprocessing.freeze_support()
    multiprocessing.set_start_method('spawn', force=True)

    parser = OptionParser(version=constant.VERSION)

    parser.add_option('-m', "--mzml_paths", type="string", default=None,
                      help="mzML file absolute path list, split by ',' ")

    parser.add_option('-t', "--txt_path", type="string", default=None,
                      help="mzML file absolute path list txt file path. Each mzml in one line")

    parser.add_option('-f', "--fasta", type="string", default=None,
                      help="Path to FASTA database.")

    parser.add_option('-c', "--config", type="string", default=None,
                      help="Path to config file.")

    parser.add_option('-o', "--output_path", type="string", default=None,
                      help="Path where search and quant results will be written. ")

    parser.add_option("--thread_num", type="int", default='20',
                      help="Run thread number.")

    parser.add_option("--fragment_min_mz", type="float",
                      default='{}'.format(config_data['database_param']['fragment_min_mz']),
                      help="Sage database fragment min mz.")
    parser.add_option("--fragment_max_mz", type="float",
                      default='{}'.format(config_data['database_param']['fragment_max_mz']),
                      help="Sage database fragment max mz.")
    parser.add_option("--enzyme_missed_cleavages", type="int",
                      default='{}'.format(config_data['database_param']['enzyme_missed_cleavages']),
                      help="Sage enzyme missed cleavages. ")
    parser.add_option("--report_psms", type="int",
                      default='{}'.format(config_data['database_param']['report_psms']),
                      help="Sage report psms. ")
    parser.add_option("--digest_min_len", type="int",
                      default='{}'.format(config_data['database_param']['digest_min_len']),
                      help="Sage enzyme min len. ")
    parser.add_option("--digest_max_len", type="int",
                      default='{}'.format(config_data['database_param']['digest_max_len']),
                      help="Sage digest max len. ")
    parser.add_option("--enzyme_cleave_at", type="string",
                      default='{}'.format(config_data['database_param']['enzyme_cleave_at']),
                      help="Sage enzyme cleave at. ")
    parser.add_option("--enzyme_restrict", type="string",
                      default='{}'.format(config_data['database_param']['enzyme_restrict']),
                      help="Sage enzyme restrict. ")
    parser.add_option("--max_variable_mods", type="int",
                      default='{}'.format(config_data['database_param']['max_variable_mods']),
                      help="Sage max variable mods. ")
    parser.add_option("--decoy_tag", type="string", default='{}'.format(config_data['database_param']['decoy_tag']),
                      help="Sage decoy tag, if not blank generate decoys is false, else is true ")
    parser.add_option("--min_matched_peaks", type="int",
                      default='{}'.format(config_data['database_param']['min_matched_peaks']),
                      help="Sage min matched peaks. ")
    parser.add_option("--max_fragment_charge", type="int",
                      default='{}'.format(config_data['database_param']['max_fragment_charge']),
                      help="Sage max matched peaks. ")

    parser.add_option("--min_precursor_charge", type="int",
                      default='{}'.format(config_data['database_param']['min_precursor_charge']),
                      help="Min precursor charge. ")
    parser.add_option("--max_precursor_charge", type="int",
                      default='{}'.format(config_data['database_param']['max_precursor_charge']),
                      help="Max precursor charge. ")

    parser.add_option("--min_precursor_tol_ppm", type="int",
                      default='{}'.format(config_data['database_param']['min_precursor_tol_ppm']),
                      help="Min precursor tol ppm. ")
    parser.add_option("--max_precursor_tol_ppm", type="int",
                      default='{}'.format(config_data['database_param']['max_precursor_tol_ppm']),
                      help="Max precursor tol ppm. ")
    parser.add_option("--min_fragment_tol_ppm", type="int",
                      default='{}'.format(config_data['database_param']['min_fragment_tol_ppm']),
                      help="Min fragment tol ppm. ")
    parser.add_option("--max_fragment_tol_ppm", type="int",
                      default='{}'.format(config_data['database_param']['max_fragment_tol_ppm']),
                      help="Max fragment tol ppm. ")

    parser.add_option("--min_peaks", type="int",
                      default='{}'.format(config_data['database_param']['min_peaks']),
                      help="Min peak. ")
    parser.add_option("--max_peaks", type="int",
                      default='{}'.format(config_data['database_param']['max_peaks']),
                      help="Max peak. ")

    parser.add_option("--model_path", type="string", default='{}'.format(config_data['model_file']),
                      help="Predict model path. ")

    parser.add_option("--cancel_mods_c", action="store_true", default=False,
                      help="Carbamidomethyl (C). ")
    parser.add_option("--cancel_mods_a", action="store_true", default=False,
                      help="Acetyl (N-term). ")
    parser.add_option("--cancel_mods_m", action="store_true", default=False,
                      help="Oxidation (M). ")
    parser.add_option("--cancel_mods_n", action="store_true", default=False,
                      help="Deamidation (N). ")
    parser.add_option("--cancel_mods_q", action="store_true", default=False,
                      help="Deamidation (Q). ")

    parser.add_option("--finetune_base_model", type="string", default='resource/model/mp_rank_00_model_states.pt',
                      help="Finetune base model. ")

    parser.add_option("--gpu_devices", type="string", default='auto',
                      help="GPU devices index list. Split by ',' ")

    (options, args) = parser.parse_args()

    config_file = options.config
    if config_file:
        # if set config file, use config file
        input_param = config_util.read_config_file(options.config)
        input_param.env = 'linux'
        input_param.model_path = options.model_path
        input_param.thread_num = options.thread_num

        mph = MainProcessHandler(input_param)
        mph.processing()
    else:
        input_param = InputParam()
        input_param.env = 'linux'
        input_param.thread_num = options.thread_num
        mzml_paths = options.mzml_paths
        txt_path = options.txt_path
        if mzml_paths:
            mzml_list = mzml_paths.split(',')
            input_param.mzml_file_list = mzml_list
        elif txt_path:
            mzml_list = file_utils.read_txt(txt_path)
            input_param.mzml_file_list = mzml_list
        else:
            raise Exception('--mzml_paths and --txt_path must get one')

        fasta = options.fasta
        if not fasta:
            raise Exception('fasta path must get')
        input_param.fasta_path = fasta

        output_path = options.output_path
        if not output_path:
            raise Exception('output path must get')
        input_param.output_dir_path = output_path

        input_param.model_path = options.model_path

        input_param.fragment_min_mz = options.fragment_min_mz
        input_param.fragment_max_mz = options.fragment_max_mz
        input_param.enzyme_missed_cleavages = options.enzyme_missed_cleavages
        input_param.report_psms = options.report_psms
        input_param.digest_min_len = options.digest_min_len
        input_param.digest_max_len = options.digest_max_len
        input_param.enzyme_cleave_at = options.enzyme_cleave_at
        input_param.enzyme_restrict = options.enzyme_restrict
        input_param.max_variable_mods = options.max_variable_mods
        input_param.decoy_tag = options.decoy_tag
        if input_param.decoy_tag and len(input_param.decoy_tag) > 0:
            input_param.generate_decoys = True
        else:
            input_param.generate_decoys = False
        input_param.min_matched_peaks = options.min_matched_peaks
        input_param.max_fragment_charge = options.max_fragment_charge

        input_param.min_precursor_charge = options.min_precursor_charge
        input_param.max_precursor_charge = options.max_precursor_charge

        input_param.min_precursor_tol_ppm = options.min_precursor_tol_ppm
        input_param.max_precursor_tol_ppm = options.max_precursor_tol_ppm

        input_param.min_fragment_tol_ppm = options.min_fragment_tol_ppm
        input_param.max_fragment_tol_ppm = options.max_fragment_tol_ppm

        input_param.min_peaks = options.min_peaks
        input_param.max_peaks = options.max_peaks

        input_param.static_mods_c = not options.cancel_mods_c
        input_param.variable_mods_a = not options.cancel_mods_a
        input_param.variable_mods_m = not options.cancel_mods_m
        input_param.variable_mods_n = not options.cancel_mods_n
        input_param.variable_mods_q = not options.cancel_mods_q

        input_param.finetune_base_model = options.finetune_base_model

        if options.gpu_devices == 'auto':
            try:
                max_member_use_rate = config_data['gpu']['max_member_use_rate']
            except Exception as e:
                max_member_use_rate = 0.5
            usage_device_list, min_free_memory = get_usage_device(max_member_use_rate)
            input_param.gpu_devices = usage_device_list
            print(usage_device_list)
            input_param.device = get_top_free_max_device(usage_device_list)
            print(f'Use device list: {input_param.gpu_devices}, max free member device is: {input_param.device}')
        else:
            usage_device_list = options.gpu_devices.split(',')
            usage_device_list = [int(nn) for nn in usage_device_list]
            input_param.gpu_devices = usage_device_list
            print(usage_device_list)
            input_param.device = get_top_free_max_device(usage_device_list)
            print(f'Use device list: {input_param.gpu_devices}, max free member device is: {input_param.device}')

        mph = MainProcessHandler(input_param)
        mph.processing()
