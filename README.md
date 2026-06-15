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

## Large Dataset Setup (Not Tracked In Git)

GitHub blocks files above 100 MB, so the full COCONUT dataset files are intentionally not committed:
- `data/raw/coconut_csv-04-2026.csv`
- `data/raw/coconut_csv-04-2026.zip`

Expected location:
- Place both files under `data/raw/`.

Expected command check from project root:

```bash
ls -lh data/raw/coconut_csv-04-2026.csv data/raw/coconut_csv-04-2026.zip
```

If these files are missing, Notebook 02 downstream steps that depend on full COCONUT data will not run.


## Configuration Scenarios
| Scenario | Key flags | What it does |
|---|---|---|
| Notebook 02 first-time build | `RUN_FRAGMENT_EXTRACTION=True`, `RUN_COCONUT_PREP=True`, `RUN_FULL_SIMILARITY=True`, `REBUILD_MOLECULES=True`, `REBUILD_MOLECULES_FILTERED=True`, `REBUILD_SIMILARITY_SCORES=True` | Builds molecule and similarity caches from scratch. |
| Notebook 02 cache reuse | `RUN_FRAGMENT_EXTRACTION=True`, `RUN_COCONUT_PREP=True`, `RUN_FULL_SIMILARITY=True`, `REBUILD_MOLECULES=False`, `REBUILD_MOLECULES_FILTERED=False`, `REBUILD_SIMILARITY_SCORES=False` | Loads existing caches and avoids rebuilding them. |
| Notebook 03 first-time build | `P3_K_RANGE=range(2,15)`, `P3_EVAL_SAMPLE_FRAC=0.10`, `P3_K_VALUES=[2,4]`, `P3_TOP_SIMILARITY_FRACTION=0.10`, `P3_UMAP_SAMPLE_SIZE=10000` | Builds default unsupervised diagnostics and plots from scratch. |
| Notebook 03 cache reuse | `SIMILARITY_SCORES_CACHE_FILE="similarity_scores.pkl"`, `P3_SIMILARITY_COL="sim_r1_tv_frag"` | Loads existing similarity-score cache as Phase 3 input. |
| Notebook 03 reproducible rerun | `RANDOM_STATE=42`, `P3_K_VALUES=[2,4]`, `P3_UMAP_N_NEIGHBORS=50`, `P3_UMAP_MIN_DIST=0.5` | Reuses fixed clustering and UMAP settings for stable reruns. |
| Notebook 03 debug / iteration | `P3_EVAL_SAMPLE_FRAC=0.05`, `P3_UMAP_SAMPLE_SIZE=3000`, `P3_K_RANGE=range(2,10)` | Runs a reduced configuration for faster iteration. |
| Notebook 05 first-time build | `DOCKING_RUN_MODE="production"`, `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=False`, `DOCKING_RECOMPUTE_LIGANDS=True`, `DOCKING_SKIP_RUN=False` | Builds production candidate and ligand caches, then runs docking. |
| Notebook 05 cache reuse | `DOCKING_RUN_MODE="production"`, `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=True`, `DOCKING_RECOMPUTE_LIGANDS=False`, `DOCKING_SKIP_RUN=False` | Loads cached candidates and ligands and avoids rebuilding them. |
| Notebook 05 reproducible rerun | `DOCKING_RUN_MODE="production"`, `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=True`, `DOCKING_RECOMPUTE_LIGANDS=False`, `DOCKING_SKIP_RUN=False` | Reuses fixed candidates and ligand cache for stable reruns. |
| Notebook 05 debug ligand prep only (random 25) | `DOCKING_RUN_MODE="debug"`, `DOCKING_RANDOM_TEST_SET_SIZE=25`, `DOCKING_SAMPLING_METHOD="random"`, `DOCKING_RECOMPUTE_LIGANDS=True`, `DOCKING_SKIP_RUN=True` | Prepares ligands only (no docking) using a random debug subset. |
| Notebook 05 debug docking from prepared ligands (random 25) | `DOCKING_RUN_MODE="debug"`, `DOCKING_RANDOM_TEST_SET_SIZE=25`, `DOCKING_SAMPLING_METHOD="random"`, `DOCKING_LOCK_FINAL_TEST_SET=True`, `DOCKING_RECOMPUTE_LIGANDS=False`, `DOCKING_SKIP_RUN=False` | Runs docking on previously prepared random debug ligands. |
| Notebook 05 debug ligand prep only (stratified 25) | `DOCKING_RUN_MODE="debug"`, `DOCKING_N_SAMPLES_TOTAL=25`, `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=False`, `DOCKING_RECOMPUTE_LIGANDS=True`, `DOCKING_SKIP_RUN=True` | Prepares ligands only using stratified sampling from Cell 3 logic. |
| Notebook 05 debug docking from prepared ligands (stratified 25) | `DOCKING_RUN_MODE="debug"`, `DOCKING_N_SAMPLES_TOTAL=25`, `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=True`, `DOCKING_RECOMPUTE_LIGANDS=False`, `DOCKING_SKIP_RUN=False` | Runs docking on previously prepared stratified debug ligands. |
| Notebook 05 debug / iteration baseline | `DOCKING_RUN_MODE="debug"`, `DOCKING_LOCK_FINAL_TEST_SET=False`, `DOCKING_RECOMPUTE_LIGANDS=True`, `DOCKING_SKIP_RUN=True` | Runs a reduced workflow for faster iteration. |
| Shared preprocessing switch | `SKIP_COCONUT_LOADING=True` | Skips COCONUT loading and validation in shared preprocessing utilities. |

### Notebook 03 Parameters Used Directly
| Parameter | Used for |
|---|---|
| `SIMILARITY_SCORES_CACHE_FILE` | Input cache file loaded for similarity-ranked candidates. |
| `RANDOM_STATE` | Seed for K-Means, sampling, and UMAP reproducibility. |
| `P3_K_RANGE` | Candidate k values evaluated in elbow and silhouette diagnostics. |
| `P3_EVAL_SAMPLE_FRAC` | Fraction of rows sampled for K-Means diagnostics. |
| `P3_K_VALUES` | Primary k values used for final K-Means cluster plots. |
| `P3_SIMILARITY_COL` | Similarity column used for top-similarity overlays. |
| `P3_TOP_SIMILARITY_FRACTION` | Fraction used to define top-similarity subset in overlays. |
| `P3_UMAP_SAMPLE_SIZE` | Subsample size used in UMAP embedding. |
| `P3_UMAP_N_NEIGHBORS` | UMAP neighborhood size. |
| `P3_UMAP_MIN_DIST` | UMAP minimum distance parameter. |



## Run Notebook 01

1. Open `notebooks/01_exploratory_analysis.ipynb`.
2. Select kernel `rapalog-py311`.
3. Run cells in order.

Key toggle in `src/config.py`:
- `SKIP_COCONUT_LOADING = False` for full analysis

For a quick dry run:
- set `SKIP_COCONUT_LOADING = True`

## Run Notebook 02

1. Open `notebooks/02_fragment_and_similarity.ipynb`.
2. Select kernel `rapalog-py311`.
3. Run cells in order.

Current defaults in `src/config.py` are intentionally conservative:
- `RUN_FRAGMENT_EXTRACTION = True`
- `RUN_COCONUT_PREP = False`
- `RUN_FULL_SIMILARITY = False`

Enable heavy steps only when ready to load/process full candidate data.

## Run Notebook 03

1. Open `notebooks/03_unsupervised_analysis.ipynb`.
2. Select kernel `rapalog-py311`.
3. Run cells in order.

Recommended defaults in `src/config.py`:
- `P3_K_VALUES = [2, 4]`
- `P3_EVAL_SAMPLE_FRAC = 0.10`
- `P3_UMAP_SAMPLE_SIZE = 10000`

For faster iteration:
- lower `P3_EVAL_SAMPLE_FRAC` (for example `0.05`)
- lower `P3_UMAP_SAMPLE_SIZE` (for example `3000`)
- narrow `P3_K_RANGE` (for example `range(2, 10)`)

## Run Notebook 04

1. Open `notebooks/04_similarity_prediction.ipynb`.
2. Select kernel `rapalog-py311`.
3. Run cells in order.

Notebook 04 depends on:
- `SIMILARITY_SCORES_CACHE_FILE` (produced upstream in Notebook 02)
- `RANDOM_STATE` (used for train/test split reproducibility)

## Run Notebook 05

1. Open `notebooks/05_docking_pipeline.ipynb`.
2. Select kernel `rapalog-py311`.
3. Run cells in order.

Typical production settings in `src/config.py`:
- `DOCKING_RUN_MODE = "production"`
- `DOCKING_SAMPLING_METHOD = "stratified"`
- `DOCKING_LOCK_FINAL_TEST_SET = False` (first run), then `True` for reruns
- `DOCKING_RECOMPUTE_LIGANDS = True` (first run), then `False` for reruns
- `DOCKING_SKIP_RUN = False`

## Notes On Docking Packages

`vina` and `meeko` are installed from `conda-forge` in `environment.yml` to avoid
source-build failures (notably Boost-related errors on pip builds).

If you ever need newer versions, prefer conda first:

```bash
conda install -n rapalog-py311 -c conda-forge vina meeko
```


