# TODO

## Refactoring

### Consolidate duplicate default parameter definitions
`ConfigVars` (class body in `ampclasslib/ac_util.py`) and `DEFAULT_CONFIG` (dict in `ampclasslib/config_params.py`) define the same 35 parameters with the same values. The two copies can silently drift when a parameter is added or changed in one place but not the other — as happened with `fb_dist_cut` and `tid_max_inner_cn` in June 2026, which required manual edits to both files.

The fix is to derive `ConfigVars` defaults from `DEFAULT_CONFIG` rather than hardcoding them twice. Something like:

```python
# in ac_util.py
from ampclasslib.config_params import DEFAULT_CONFIG

class ConfigVars:
    for _k, _v in DEFAULT_CONFIG.items():
        locals()[_k] = _v
```

Care needed: `set_config_vars()` (also in `ac_util.py`) iterates over `config` keys and calls `setattr(ConfigVars, key, ...)`, which works the same regardless of how the class body is populated.
