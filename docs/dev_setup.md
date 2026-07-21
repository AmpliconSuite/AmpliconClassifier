# Developer Setup

Target: Linux or macOS with conda installed. BFBArchitect is a core dependency; Gurobi is an optional accelerator for its integrated BFB detection.

## Assumptions

- AmpliconClassifier has been cloned locally.
- `$AC_SRC` points to this repository.
- `$AA_DATA_REPO` points to a populated AmpliconArchitect data repository.
- For BFBArchitect development, clone its repository locally; ordinary package installs normally resolve BFBArchitect from PyPI.

Example shell setup:

```bash
export AC_SRC=/path/to/AmpliconClassifier
export AA_DATA_REPO=/path/to/AmpliconArchitect/data_repo
```

---

## 1. Create the conda environment

```bash
conda create -n bfb_ac_dev -c conda-forge -c bioconda python=3.10 \
  intervaltree scipy pandas pysam numpy \
  matplotlib-base cnvkit flask -y
```

Then activate and install pip-only packages:

```bash
conda activate bfb_ac_dev
pip install gurobipy pulp
```

---

## 2. Install BFBArchitect (editable, from source)

BFBArchitect is a core AmpliconClassifier dependency and is normally installed automatically from PyPI. For development, install it as an editable checkout instead so local changes to BFBArchitect are picked up immediately:

```bash
pip install -e /path/to/BFBArchitect
```

Verify:

```bash
python -c "from bfbarchitect.reconstruct import reconstruct_bfb_from_graph; print('BFBArchitect OK')"
```

---

## 3. Verify AmpliconClassifier

```bash
cd "$AC_SRC"
conda run -n bfb_ac_dev python amplicon_classifier.py --help
```

Run the unit tests to confirm the environment is intact:

```bash
conda run -n bfb_ac_dev python -m unittest discover -s tests
```

---

## 4. Verify Gurobi

```bash
conda run -n bfb_ac_dev python -c "import gurobipy; m = gurobipy.Model(); print('Gurobi OK')"
```

If the license isn't found, check that `~/gurobi.lic` exists or `$GRB_LICENSE_FILE` is set.

---

## Notes

- Always run AC with `conda run -n bfb_ac_dev` or activate the environment first.
- BFBArchitect is installed as an editable install, so pulling updates to the repo is sufficient (no reinstall needed unless `pyproject.toml` dependencies change).
- If BFBArchitect dependencies change, re-run `pip install -e /path/to/BFBArchitect` after pulling.
