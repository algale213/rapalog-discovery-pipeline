# Rapalog Discovery Pipeline

## Project Overview

This project develops a reproducible computational pipeline for identifying natural products with chemical features and predicted FKBP12 binding behavior related to the FKBP12-binding region of rapamycin. The workflow combines molecular filtering, fragment-based similarity analysis, unsupervised learning, molecular docking, and supervised regression.

The project was developed from coursework in computational chemistry, cheminformatics, and machine learning and has been reorganized into modular Python functions and sequential Jupyter notebooks.

## Scientific Motivation

Rapamycin binds FKBP12 before interacting with the FKBP–rapamycin-binding (FRB) domain of mTOR complex 1. This project investigates whether natural products contain chemical features resembling the FKBP12-facing portion of rapamycin and whether selected candidates are predicted to bind favorably within the FKBP12 pocket.

The current workflow evaluates FKBP12 binding only. It does not establish formation of the complete FKBP12–compound–FRB ternary complex or biological modulation of mTOR signaling.

## Pipeline

1. Explore and quality-check the COCONUT natural-product dataset.
2. Extract the FKBP12-facing region of rapamycin from PDB structure 1FKB.
3. Filter compounds using physicochemical and structural criteria.
4. Calculate Morgan fingerprints and Tanimoto/Tversky similarity.
5. Characterize chemical space using principal component analysis, K-means clustering, and UMAP.
6. Select a stratified candidate subset for AutoDock Vina docking.
7. Model molecular similarity and Vina-derived estimated binding affinity using linear and regularized regression.
8. Save figures, tables, and reproducible intermediate results.

## Verified Production Results

- Approximately 740,000 COCONUT entries were evaluated.
- 96,682 compounds passed the project’s physicochemical filters.
- 1,000 compounds were selected for production docking using cluster- and similarity-aware stratified sampling.
- 984 compounds successfully completed ligand preparation and AutoDock Vina docking.
- The production run used a 25 × 25 × 25 Å docking box, exhaustiveness of 8, and up to 9 poses per ligand.
- 213 compounds received Vina scores below −8.0 kcal/mol.
- The strongest-scoring candidates included NPD14003, STL518675, CNP0253821.1, and Pentacecilide B.

Vina scores are computational estimates used for ranking. They are not experimentally measured binding affinities or validated free energies of binding.

## Repository Structure

| Path | Purpose |
|---|---|
| `notebooks/` | Sequential exploratory, similarity, clustering, docking, and regression analyses |
| `src/` | Reusable data-processing, fragment-extraction, docking, and visualization functions |
| `data/raw/` | User-supplied structures and external source data |
| `data/processed/` | Local intermediate checkpoints, excluded from Git |
| `results/` | Selected figures, tables, and production results |
| `docs/` | Methodology, biological interpretation, results, and future work |
| `run_pipeline.py` | Environment validation and notebook execution entry point |

## Important Limitations

This is a computational screening study. High-ranking candidates require independent validation through repeated docking, alternative receptor conformations, molecular dynamics where appropriate, biochemical binding assays, and functional experiments. Similarity to rapamycin or a favorable FKBP12 docking score does not establish rapamycin-like biological activity.

## Python Compatibility

This project targets Python 3.11 for compatibility with the tested RDKit, AutoDock Vina, Meeko, SHAP, and UMAP package stack.

## Create The Environment

From the project root:

```bash
conda env create -f environment.yml
conda activate rapalog-py311
python -m ipykernel install --user --name rapalog-py311 --display-name "rapalog-py311"
```

This repository is meant to stay portable across machines, so use the shared `environment.yml` instead of hardcoding a local interpreter path in the project. After creating the environment, select the `rapalog-py311` kernel in VS Code for notebooks.

If the environment already exists:

```bash
conda env update -f environment.yml --prune
conda activate rapalog-py311
```

## Validation and Tests

Validate the environment and required project files:

```bash
python run_pipeline.py --check --require-coconut


## Large Dataset Setup (Not Tracked In Git)

GitHub blocks files above 100 MB, so the full COCONUT dataset files are intentionally not committed:
- `data/raw/coconut_csv-04-2026.csv`
- `data/raw/coconut_csv-04-2026.zip`

Expected location:
- Place `coconut_csv-04-2026.csv` under `data/raw/`. Retaining the downloaded ZIP archive is optional.

Expected command check from project root:

```bash
ls -lh data/raw/coconut_csv-04-2026.csv
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
| Notebook 05 first-time build | `DOCKING_RUN_MODE="production"`, `DOCKING_N_SAMPLES_TOTAL=1000`,  `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=False`, `DOCKING_RECOMPUTE_LIGANDS=True`, `DOCKING_SKIP_RUN=False` | Builds a stratified set of 1,000 production candidates, prepares the ligands, and runs docking. |
| Notebook 05 cache reuse | `DOCKING_RUN_MODE="production"`, `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=True`, `DOCKING_RECOMPUTE_LIGANDS=False`, `DOCKING_SKIP_RUN=True` | Loads the fixed candidates, prepared ligands, and completed docking checkpoint without rerunning docking. |
| Notebook 05 reproducible rerun | `DOCKING_RUN_MODE="production"`, `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=True`, `DOCKING_RECOMPUTE_LIGANDS=False`, `DOCKING_SKIP_RUN=False` | Reuses fixed candidates and ligand cache for stable reruns. |
| Notebook 05 debug ligand prep only (random 25) | `DOCKING_RUN_MODE="debug"`, `DOCKING_RANDOM_TEST_SET_SIZE=25`, `DOCKING_SAMPLING_METHOD="random"`, `DOCKING_RECOMPUTE_LIGANDS=True`, `DOCKING_SKIP_RUN=True` | Prepares ligands only (no docking) using a random debug subset. |
| Notebook 05 debug docking from prepared ligands (random 25) | `DOCKING_RUN_MODE="debug"`, `DOCKING_RANDOM_TEST_SET_SIZE=25`, `DOCKING_SAMPLING_METHOD="random"`, `DOCKING_LOCK_FINAL_TEST_SET=True`, `DOCKING_RECOMPUTE_LIGANDS=False`, `DOCKING_SKIP_RUN=False` | Runs docking on previously prepared random debug ligands. |
| Notebook 05 debug ligand prep only (stratified 25) | `DOCKING_RUN_MODE="debug"`, `DOCKING_N_SAMPLES_TOTAL=25`, `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=False`, `DOCKING_RECOMPUTE_LIGANDS=True`, `DOCKING_SKIP_RUN=True` | Prepares ligands only using stratified sampling from Cell 3 logic. |
| Notebook 05 debug docking from prepared ligands (stratified 25) | `DOCKING_RUN_MODE="debug"`, `DOCKING_N_SAMPLES_TOTAL=25`, `DOCKING_SAMPLING_METHOD="stratified"`, `DOCKING_LOCK_FINAL_TEST_SET=True`, `DOCKING_RECOMPUTE_LIGANDS=False`, `DOCKING_SKIP_RUN=False` | Runs docking on previously prepared stratified debug ligands. |
| Notebook 05 debug / iteration baseline | `DOCKING_RUN_MODE="debug"`, `DOCKING_LOCK_FINAL_TEST_SET=False`, `DOCKING_RECOMPUTE_LIGANDS=True`, `DOCKING_SKIP_RUN=True` | Runs a reduced workflow for faster iteration. |
| Notebook 06 first-time build (debug) | `DOCKING_RUN_MODE="debug"`, `BINDING_TARGET_COL="vina_score"`, `BINDING_TEST_SIZE=0.20`, `BINDING_P6_K_RANGE=range(2,15)` | Loads debug Phase 5 docking results, builds Phase 6/7 artifacts, trains all models, and saves outputs to `results/notebook_06/debug/`. |
| Notebook 06 first-time build (production) | `DOCKING_RUN_MODE="production"`, `BINDING_TARGET_COL="vina_score"`, `BINDING_TEST_SIZE=0.20`, `BINDING_P6_K_RANGE=range(2,15)` | Loads production Phase 5 docking results and saves outputs to `results/notebook_06/prod/`. |
| Notebook 06 split/regularization tuning | `BINDING_TEST_SIZE`, `BINDING_P7_ALPHA_LOGSPACE_START`, `BINDING_P7_ALPHA_LOGSPACE_STOP`, `BINDING_P7_ALPHA_LOGSPACE_NUM`, `BINDING_P7_ELASTIC_L1_RATIOS` | Adjusts holdout split and ElasticNet/Ridge/Lasso hyperparameter search without editing notebook code. |
| Notebook 06 PCA/K diagnostics tuning | `BINDING_P6_PCA_CUMVAR_THRESHOLD`, `BINDING_P6_K_RANGE` | Adjusts PCA threshold line and KMeans evaluation range in Phase 6 diagnostics. |
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

### Notebook 06 Parameters Used Directly
| Parameter | Used for |
|---|---|
| `DOCKING_RUN_MODE` | Selects debug/production folder suffix and input source from Notebook 05 results. |
| `BINDING_TARGET_COL` | Target column used for binding-affinity regression (`vina_score`). |
| `BINDING_SIMILARITY_COL` | Similarity feature column used in Phase 6/7 modeling. |
| `BINDING_TEST_SIZE` | Holdout fraction for stratified train/test split in Phase 7. |
| `RANDOM_STATE` | Reproducibility seed for split logic and CV models. |
| `BINDING_P6_PCA_CUMVAR_THRESHOLD` | Reference line threshold in PCA cumulative variance plot. |
| `BINDING_P6_K_RANGE` | Evaluated k-range for elbow and silhouette diagnostics. |
| `BINDING_P7_ALPHA_LOGSPACE_START` | Lower log10 bound for regularization alpha grid. |
| `BINDING_P7_ALPHA_LOGSPACE_STOP` | Upper log10 bound for regularization alpha grid. |
| `BINDING_P7_ALPHA_LOGSPACE_NUM` | Number of alpha values in regularization grid. |
| `BINDING_P7_ELASTIC_L1_RATIOS` | Elastic Net L1/L2 mixture candidates for cross-validation. |



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
- `DOCKING_N_SAMPLES_TOTAL = 1000`
- `DOCKING_SAMPLING_METHOD = "stratified"`
- `DOCKING_LOCK_FINAL_TEST_SET = False` (first run), then `True` for reruns
- `DOCKING_RECOMPUTE_LIGANDS = True` (first run), then `False` for reruns
- `DOCKING_SKIP_RUN = False`

## Run Notebook 06

1. Open `notebooks/06_binding_affinity_prediction.ipynb`.
2. Select kernel `rapalog-py311`.
3. Run cells in order.

Notebook 06 depends on Phase 5 output:
- `results/notebook_05/debug/docking_results_full.csv` when `DOCKING_RUN_MODE = "debug"`
- `results/notebook_05/prod/docking_results_full.csv` when `DOCKING_RUN_MODE = "production"`

Typical settings in `src/config.py`:
- `BINDING_TARGET_COL = "vina_score"`
- `BINDING_SIMILARITY_COL = DOCKING_SIMILARITY_COL`
- `BINDING_TEST_SIZE = 0.20`
- `BINDING_P6_K_RANGE = range(2, 15)`
- `BINDING_P6_PCA_CUMVAR_THRESHOLD = 0.85`

## Notes On Docking Packages

`vina` and `meeko` are installed from `conda-forge` in `environment.yml` to avoid
source-build failures (notably Boost-related errors on pip builds).

If you ever need newer versions, prefer conda first:

```bash
conda install -n rapalog-py311 -c conda-forge vina meeko
```


