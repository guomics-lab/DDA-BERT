import importlib.resources
import json
import logging
import logging.config as log_config
import os
import os.path
import time

cwd = os.getcwd()
log_dir = os.path.join(cwd, 'logs')
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

try:
    with importlib.resources.open_text("dda_bert.config", 'log.config') as f:
        config = json.load(f)
        log_config.dictConfig(config)
except Exception:
    logging.basicConfig(level=logging.DEBUG)

logger = logging.getLogger()
logger.setLevel(logging.INFO)


def get_current_logger(log_dir='logs'):
    os.makedirs(log_dir, exist_ok=True)

    log_format = logging.Formatter('DDA-BERT: %(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
    current_logger = logging.getLogger('DDA-BERT')
    log_file_path = os.path.join(log_dir, 'DDA-BERT_{}.log'.format(time.time_ns()))
    fh = logging.FileHandler(log_file_path, mode='a', encoding='utf-8')
    fh.setFormatter(log_format)
    fh.setLevel(logging.INFO)
    current_logger.addHandler(fh)
    return log_file_path, current_logger
