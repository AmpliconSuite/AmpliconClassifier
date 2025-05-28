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
    "tot_min_del": 5000,      # non-trivial deletion size
    "minCycleSize": 5000,     # minimum size of AA cycle to consider as valid
    "compCycContCut": 50000,  # minimum total complex cycle size for classifying complex cycles
    "anyCycContcut": 10000,   # minimum total complex path/cycle size for complex amplicon
    "ampLenOverMinCN": 5000,  # amount of amplicon over minimum CN to trigger focal amp call
    "cycCut": 0.12,           # minimum proportion cycle weight in cyclic path for ecDNA
    "compCut": 0.3,           # minimum proportion cycle weight in complex cyclic path for complex ecDNA
    "min_amp_cn": 4.5,        # minimum CN for focal amp
    "sig_amp": 7,             # minimum CN for significantly amplified focal amp
    "high_amp": 12,           # minimum CN for high CN focal amp
    "max_segdup_size": 1000000,  # maximum allowed size for something to be a segmental dup
    "segdup_max_extra_fraction": 0.25,  # maximum additional CN ratio beyond baseline 2x for a segmental dup
    "decomposition_strictness": 0.1,    # for singleton paths/cycles, strictness scale of filtering path flow against CN
    "min_flow": 1.0,          # minimum flow to consider path for classification as valid focal amp

    # bfb-related items
    "min_fb_read_prop": 0.25,          # min proportion of SV reads in foldbacks to call BFB
    "fb_break_weight_prop": 0.3,       # min proportion of utilized path/cycle SVs supporting BFB foldbacks (weighted by flow)
    "fb_dist_cut": 25000,              # max distance between ends for inversion to be foldback
    "max_nonbfb_break_weight": 0.5,    # max proportion of utilized path/cycle SVs supporting non-foldback connections
    "min_bfb_cycle_weight_ratio": 0.6  # minimum flow*length weighted proportion of BFB-like paths/cycles for BFB
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
    """Create the default config file if it doesn't exist or if it differs from DEFAULT_CONFIG."""
    config_path = get_default_config_path()

    # Always write if file doesn't exist
    should_write = not os.path.exists(config_path)

    # If file exists, check if it matches DEFAULT_CONFIG
    if not should_write:
        try:
            with open(config_path, 'r') as f:
                existing_config = json.load(f)

            # Compare existing config with DEFAULT_CONFIG
            # Check if any values differ or if keys are missing/extra
            if (existing_config != DEFAULT_CONFIG or
                    set(existing_config.keys()) != set(DEFAULT_CONFIG.keys())):
                should_write = True

        except (json.JSONDecodeError, IOError):
            # If we can't read/parse the file, rewrite it
            should_write = True

    # Write the default config if needed
    if should_write:
        with open(config_path, 'w') as f:
            json.dump(DEFAULT_CONFIG, f, indent=4)

        print("Updated default config file")

    return config_path


def validate_config(config):
    """
    Ensure config has all required keys and validate values.

    Args:
        config (dict): Configuration parameters - will be modified in place

    Raises:
        ValueError: If any value is invalid
    """
    # Ensure all required keys exist, fill with defaults if missing
    for key, default_value in DEFAULT_CONFIG.items():
        if key not in config:
            config[key] = default_value

    # Validate numeric values are non-negative
    numeric_params = list(DEFAULT_CONFIG.keys())
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


def load_config(config_file=None, outloc=None):
    """
    Load configuration from file and validate it.

    Args:
        config_file (str, optional): Path to custom config file. If None, uses default.
        outloc (str, optional): Output location (currently unused but kept for compatibility)

    Returns:
        dict: Validated configuration
    """
    if config_file is None:
        # This will create/update the default config if needed
        config_path = create_default_config_file()
    else:
        config_path = str(config_file) if HAS_PATHLIB and isinstance(config_file, Path) else config_file
        if not os.path.exists(config_path):
            raise FileNotFoundError("Config file not found: {}".format(config_path))

    # Load the config file
    with open(config_path, 'r') as f:
        config = json.load(f)

    # Validate and fill in any missing values
    validate_config(config)

    return config