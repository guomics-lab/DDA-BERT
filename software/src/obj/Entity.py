class MzMLFileInfo(object):

    def __init__(self):
        self.file_path = None
        self.file_name = None
        self.base_file_name = None


class InputParam(object):

    def __init__(self):
        self.mzml_file_list = []
        self.fasta_path = None
        self.output_dir_path = None
        self.model_path = None
        self.env = None

        self.batch_size = None

        self.logger = None
        self.logger_path = None

        self.fragment_min_mz = None
        self.fragment_max_mz = None
        self.enzyme_missed_cleavages = None
        self.report_psms = None
        self.digest_min_len = None
        self.digest_max_len = None
        self.enzyme_cleave_at = None
        self.enzyme_restrict = None
        self.max_variable_mods = None
        self.generate_decoys = None
        self.decoy_tag = None
        self.min_peaks = None
        self.max_peaks = None
        self.min_matched_peaks = None
        self.max_fragment_charge = None

        self.min_precursor_charge = None
        self.max_precursor_charge = None

        self.min_precursor_tol_ppm = None
        self.max_precursor_tol_ppm = None

        self.min_fragment_tol_ppm = None
        self.max_fragment_tol_ppm = None

        self.static_mods_c = True
        self.variable_mods_a = True
        self.variable_mods_m = True
        self.variable_mods_n = False
        self.variable_mods_q = False

        self.open_finetune_process = False

        self.thread_num = 20

        self.device = None
        self.gpu_devices = []

        self.open_sage = True
        self.open_alphapept = True
        self.open_get_data = True
        self.open_pred = True
        self.open_fdr_calc = True
        self.open_finetune = True
        self.open_eval = True
        self.open_protein_infer = True
        self.open_clear_data = True

        self.finetune_base_model = None


class SageConfigParam(object):

    def __init__(self):
        self.fragment_min_mz = None
        self.fragment_max_mz = None
        self.enzyme_missed_cleavages = None
        self.report_psms = None
        self.digest_min_len = None
        self.digest_max_len = None
        self.enzyme_cleave_at = None
        self.enzyme_restrict = None
        self.max_variable_mods = None
        self.generate_decoys = None
        self.decoy_tag = None
        self.min_matched_peaks = None
        self.max_fragment_charge = None

        self.min_precursor_charge = None
        self.max_precursor_charge = None

        self.min_precursor_tol_ppm = None
        self.max_precursor_tol_ppm = None

        self.min_fragment_tol_ppm = None
        self.max_fragment_tol_ppm = None

        self.static_mods_c = None
        self.variable_mods_a = None
        self.variable_mods_m = None
        self.variable_mods_n = None
        self.variable_mods_q = None


class InfoMsg(object):

    def __init__(self, step, status, msg=None, mzml_index=None):
        '''
        :param step:
            1 - start
            20 - sage
            30 - get data
            40 - predict
            50 - protein infer
            100 - stop
        :param status:
            1 - start

            900 - stopping
            988 - success exit
            999 - exception exit

        :param msg:
        '''
        self.step = step
        self.status = status
        self.msg = msg
        self.mzml_index = mzml_index

    @staticmethod
    def json_to_object(dct):
        return InfoMsg(mzml_index=dct['mzml_index'], step=dct['step'], status=dct['status'], msg=dct['msg'])


class ProcessMsg(object):
    '''
    status
        0-wait
        1-running
        5-success
        9-error
    '''

    def __init__(self, file_index, step, status, msg=''):
        self.file_index = file_index
        self.step = step
        self.status = status
        self.msg = msg
