from enum import IntEnum

VERSION = 'DDA-BERT V2.0'

env_linux = 'linux'

env_win = 'win'


class ProgressStepEnum(IntEnum):
    START = 0
    PSM_PEAK = 5
    SAGA = 10
    ALPHAPEPT = 15
    GET_DATA = 20
    PREDICT = 30
    FDR_CALC = 40
    FINETUNE = 50
    EVAL = 60
    PROTEIN_INFER = 70
    RESULT_BUILD = 80
    CLEAR_DATA = 90
    END = 100


class ProgressStepStatusEnum(IntEnum):
    WAIT = 1
    RUNNING = 2
    SUCCESS = 3
    ERROR = 4
    IDENTIFY_NUM = 90
    END = 99
    FAIL_END = 98
    STOPPING = 900
    STOPPED = 901
    ALL_END = 999
