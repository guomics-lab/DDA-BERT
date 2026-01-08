import importlib.resources
import os
from pathlib import Path

ALPHA_PEPT_DEFAULT_CONFIG_PATH = 'config/AlphaPept.yaml'
SAGE_DEFAULT_CONFIG_PATH = 'config/Sage.json'
MODEL_DEFAULT_CONFIG_PATH = 'config/Model.yaml'
WORK_FLOW_DEFAULT_CONFIG_PATH = 'config/FP_MSBooster.workflow'

os.makedirs('config', exist_ok=True)


def ensure_sage_config():
    config_path = Path(SAGE_DEFAULT_CONFIG_PATH)

    if not config_path.exists():
        # not exist, copy default
        with importlib.resources.open_text("dda_bert.config", 'Sage.json') as f:
            default_content = f.read()
        config_path.write_text(default_content, encoding="utf-8")


def get_sage_config_path():
    ensure_sage_config()
    return SAGE_DEFAULT_CONFIG_PATH


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
        with importlib.resources.open_text("dda_bert.config", 'Model.yaml') as f:
            default_content = f.read()
        config_path.write_text(default_content, encoding="utf-8")


def get_model_config_path():
    ensure_model_config()
    return MODEL_DEFAULT_CONFIG_PATH


def ensure_work_flow_config():
    config_path = Path(WORK_FLOW_DEFAULT_CONFIG_PATH)

    if not config_path.exists():
        # not exist, copy default
        with importlib.resources.open_text("dda_bert.config", 'FP_MSBooster.workflow') as f:
            default_content = f.read()
        config_path.write_text(default_content, encoding="utf-8")


def get_work_fow_config_path():
    ensure_work_flow_config()
    return WORK_FLOW_DEFAULT_CONFIG_PATH
