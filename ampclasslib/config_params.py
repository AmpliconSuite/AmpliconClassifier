#!/usr/bin/env python3

import json
import os

try:
    from pathlib import Path

    HAS_PATHLIB = True
except ImportError:
    HAS_PATHLIB = False

# Default configuration values
DEFAULT_CONFIG = {
    "tot_min_del": 5000,
    "minCycleSize": 5000,
    "compCycContCut": 50000,
    "anyCycContcut": 10000,
    "ampLenOverMinCN": 5000,
    "cycCut": 0.12,
    "compCut": 0.3,
    "min_upper_cn": 4.5,
    "sig_amp": 7.5,
    "decomposition_strictness": 0.1,
    "min_score_for_bfb": 0.25,
    "fb_dist_cut": 25000,
    "min_flow": 1.0
}


def get_default_config_path():
    """Get the path to the default config file in the ampclasslib directory."""
    # Get the directory containing this file
    if HAS_PATHLIB:
        ampclasslib_dir = Path(__file__).parent
        return str(ampclasslib_dir / "default_config.json")
    else:
        # Fallback for older Python versions
        ampclasslib_dir = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(ampclasslib_dir, "default_config.json")


def create_default_config_file():
    """Create the default config file if it doesn't exist."""
    config_path = get_default_config_path()
    if not os.path.exists(config_path):
        with open(config_path, 'w') as f:
            json.dump(DEFAULT_CONFIG, f, indent=4)
    return config_path


def validate_config(config):
    """
    Ensure config has all required keys and validate values.

    Args:
        config (dict): Configuration parameters

    Raises:
        ValueError: If any value is invalid
    """
    # Ensure all required keys exist, fill with defaults if missing
    for key, default_value in DEFAULT_CONFIG.items():
        if key not in config:
            config[key] = default_value

    # Validate numeric values are non-negative
    numeric_params = ['tot_min_del', 'minCycleSize', 'compCycContCut', 'anyCycContcut',
                      'ampLenOverMinCN', 'min_upper_cn', 'sig_amp', 'fb_dist_cut', 'min_flow']
    for param in numeric_params:
        if config[param] < 0:
            raise ValueError("{} must be non-negative".format(param))

    # Validate percentage/ratio values
    if not 0 <= config["decomposition_strictness"] <= 1:
        raise ValueError("decomposition_strictness must be between 0 and 1")

    if not 0 <= config["cycCut"] <= 1:
        raise ValueError("cycCut must be between 0 and 1")

    if not 0 <= config["compCut"] <= 1:
        raise ValueError("compCut must be between 0 and 1")

    return True


def load_config(config_file=None):
    """
    Load configuration from file and validate it.

    Args:
        config_file (str, optional): Path to custom config file. If None, uses default.

    Returns:
        dict: Validated configuration
    """
    if config_file is None:
        config_path = create_default_config_file()
    else:
        config_path = str(config_file) if HAS_PATHLIB and isinstance(config_file, Path) else config_file
        if not os.path.exists(config_path):
            raise FileNotFoundError("Config file not found: {}".format(config_path))

    with open(config_path, 'r') as f:
        config = json.load(f)

    validate_config(config)
    return config