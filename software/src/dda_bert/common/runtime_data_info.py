class RuntimeData(object):

    def __init__(self):
        self.mzml_list = []
        self.mzml_deal_count = 0

        self.mpt = None
        self.msg_sub_thread = None

        self.start_timestamp = None
        self.current_mzml_index = None
        self.current_is_success = True
        self.current_pred_num = None
        self.current_pred_all_num = None

        self.running_flag = True


runtime_data = RuntimeData()
