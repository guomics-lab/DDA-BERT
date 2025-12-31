import importlib.resources
import os
from pathlib import Path

import yaml

ALPHA_PEPT_DEFAULT_CONFIG_PATH = 'config/AlphaPept.yaml'
DDA_BERT_DEFAULT_CONFIG_PATH = 'config/dda_bert_config.yml'
MODEL_DEFAULT_CONFIG_PATH = 'config/model.yaml'
WORK_FLOW_DEFAULT_CONFIG_PATH = 'config/FP_MSBooster.workflow'

os.makedirs('config', exist_ok=True)


def ensure_base_config():
    config_path = Path(DDA_BERT_DEFAULT_CONFIG_PATH)
    if not config_path.exists():
        # not exist, copy default
        with importlib.resources.open_text("dda_bert.config", 'dda_bert_config.yml') as f:
            default_content = f.read()
        config_path.write_text(default_content, encoding="utf-8")


def load_base_config():
    ensure_base_config()
    config = {}
    try:
        with importlib.resources.open_text("dda_bert.config", "dda_bert_config.yml") as f:
            config.update(yaml.safe_load(f) or {})
    except Exception:
        pass

    if os.path.exists(DDA_BERT_DEFAULT_CONFIG_PATH):
        with open(DDA_BERT_DEFAULT_CONFIG_PATH, encoding="utf-8") as f:
            user_config = yaml.safe_load(f) or {}
        config.update(user_config)

    return config


def ensure_alpha_config():
    config_path = Path(ALPHA_PEPT_DEFAULT_CONFIG_PATH)

    if not config_path.exists():
        # not exist, copy default
        with importlib.resources.open_text("dda_bert.config", 'AlphaPept.yaml') as f:
            default_content = f.read()
        config_path.write_text(default_content, encoding="utf-8")


def get_alpha_config_path():
    ensure_alpha_config()
    return ALPHA_PEPT_DEFAULT_CONFIG_PATH


def ensure_model_config():
    config_path = Path(MODEL_DEFAULT_CONFIG_PATH)

    if not config_path.exists():
        # not exist, copy default
        with importlib.resources.open_text("dda_bert.config", 'model.yaml') as f:
            default_content = f.read()
        config_path.write_text(default_content, encoding="utf-8")


def get_model_config_path():
    ensure_model_config()
    return MODEL_DEFAULT_CONFIG_PATH


def ensure_work_flow_config():
    config_path = Path(WORK_FLOW_DEFAULT_CONFIG_PATH)

    if not config_path.exists():
        # not exist, copy default
        with importlib.resources.open_text("dda_bert.config", 'model.yaml') as f:
            default_content = f.read()
        config_path.write_text(default_content, encoding="utf-8")


def get_work_fow_config_path():
    ensure_work_flow_config()
    return WORK_FLOW_DEFAULT_CONFIG_PATH
