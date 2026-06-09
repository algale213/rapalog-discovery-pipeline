# Rapalog Discovery Pipeline

## Python Compatibility

This project targets Python 3.11 for compatibility with the package stack used in
`Project2_Chem277A_v9.ipynb` (notably RDKit, Vina, Meeko, SHAP, and UMAP).

Recommended:
- Python 3.11

Avoid for this workflow:
- Python 3.13+
- Python 3.14

## Create The Environment

From the project root:

```bash
conda env create -f environment.yml
conda activate rapalog-py311
python -m ipykernel install --user --name rapalog-py311 --display-name "rapalog-py311"
```

If the environment already exists:

```bash
conda env update -f environment.yml --prune
conda activate rapalog-py311
```

## Run Notebook 02

1. Open `notebooks/02_fragement_and_similarity.ipynb`.
2. Select kernel `rapalog-py311`.
3. Run cells in order.

Current defaults in `src/config.py` are intentionally conservative:
- `RUN_FRAGMENT_EXTRACTION = True`
- `RUN_COCONUT_PREP = False`
- `RUN_FULL_SIMILARITY = False`

Enable heavy steps only when ready to load/process full candidate data.

## Notes On Docking Packages

`vina` and `meeko` are installed from `conda-forge` in `environment.yml` to avoid
source-build failures (notably Boost-related errors on pip builds).

If you ever need newer versions, prefer conda first:

```bash
conda install -n rapalog-py311 -c conda-forge vina meeko
```# rapalog-discovery-pipeline
