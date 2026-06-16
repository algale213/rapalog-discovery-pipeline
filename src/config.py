"""Centralized runtime flags and constants for the rapalog pipeline."""
# Shared / cross-notebook flags
# Used by shared src/ modules and multiple notebooks
SKIP_COCONUT_LOADING = True

# Pipeline phase controls
RUN_EDA = True
RUN_SIMILARITY = True
RUN_CLUSTERING = True
RUN_SIMILARITY_MODEL = True
RUN_DOCKING = True
RUN_BINDING_MODEL = True

# Notebook 02: fragment + similarity workflow
# Scenario: build or reuse the molecule and similarity caches that rank docking candidates
RUN_FRAGMENT_EXTRACTION = True
RUN_COCONUT_PREP = True
RUN_FULL_SIMILARITY = True

# Notebook 02 cache rebuild flags
# REBUILD_MOLECULES           -> raw RDKit molecule cache
# REBUILD_MOLECULES_FILTERED  -> filtered molecule cache
# REBUILD_SIMILARITY_SCORES   -> similarity score cache
REBUILD_MOLECULES = False
REBUILD_MOLECULES_FILTERED = False
REBUILD_SIMILARITY_SCORES = False

# Notebook 02 cache file names (resolved relative to project root)
MOLECULES_CACHE_FILE = "rdkit_coconut_molecules.pkl"
FILTERED_MOLECULES_CACHE_FILE = "rdkit_coconut_molecules_filtered.pkl"
SIMILARITY_SCORES_CACHE_FILE = "similarity_scores.pkl"

# Notebook 02 prefilter thresholds
RAPAMYCIN_PREFILTERS = {
	"mw_min": 350,
	"mw_max": 700,
	"hba_min": 3,
	"hba_max": 8,
	"hbd_min": 1,
	"hbd_max": 4,
	"logp_min": 2,
	"logp_max": 5,
	"rings_min": 2,
	"rings_max": 5,
	"tpsa_max": 150,
}

# Notebook 02 similarity / ranking defaults
TOP_N_CANDIDATES = 5000
RANDOM_STATE = 42

# Notebook 03: unsupervised analysis high-value controls
P3_K_RANGE = range(2, 15)
P3_EVAL_SAMPLE_FRAC = 0.10
P3_K_VALUES = [2, 4]
P3_SIMILARITY_COL = "sim_r1_tv_frag"
P3_TOP_SIMILARITY_FRACTION = 0.10
P3_UMAP_SAMPLE_SIZE = 10000
P3_UMAP_N_NEIGHBORS = 50
P3_UMAP_MIN_DIST = 0.5

# Notebook 05: docking workflow
# debug -> smaller / faster iteration
# production -> full run with fixed candidate set and full ligand prep
DOCKING_N_SAMPLES_TOTAL = 25
DOCKING_RUN_MODE = "debug"  # Run in one of: "debug" or "production"
DOCKING_SIMILARITY_COL = "sim_r1_tv_frag"	# Selects similarity scoring method/hyperparameters 
DOCKING_SAMPLING_METHOD = "random" if DOCKING_RUN_MODE == "debug" else "stratified"  # "random" or "stratified"
DOCKING_RANDOM_TEST_SET_SIZE = 25
DOCKING_LOCK_FINAL_TEST_SET = False  # True reuses the cached candidate set; False rebuilds it

# True regenerates ligand PDBQT cache; False reuses it if present.
DOCKING_RECOMPUTE_LIGANDS = False
DOCKING_LIGAND_PREP_PROGRESS_EVERY = 5
DOCKING_SKIP_RUN = False

# Notebook 05 receptor inputs / outputs
DOCKING_PDB_FILE = "1FKB.pdb"
DOCKING_CLEAN_PDB_FILE = "1FKB_clean.pdb"
DOCKING_RECEPTOR_PDBQT_FILE = "fkbp12.pdbqt"

# Notebook 05 mode-specific cache filenames
DOCKING_CANDIDATES_FILE = "docking_candidates.pkl"
DOCKING_LIGANDS_PICKLE_FILE = "ligand_pdbqt_data.pkl"
DOCKING_CHECKPOINT_FILE = "docking_checkpoint.pkl"
DOCKING_CANDIDATES_FILE_DEBUG = "docking_candidates_debug.pkl"
DOCKING_CANDIDATES_FILE_PRODUCTION = "docking_candidates_prod.pkl"
DOCKING_LIGANDS_PICKLE_FILE_DEBUG = "ligand_pdbqt_data_debug.pkl"
DOCKING_LIGANDS_PICKLE_FILE_PRODUCTION = "ligand_pdbqt_data_prod.pkl"
DOCKING_CHECKPOINT_FILE_DEBUG = "docking_checkpoint_debug.pkl"
DOCKING_CHECKPOINT_FILE_PRODUCTION = "docking_checkpoint_prod.pkl"

# Notebook 05 docking parameters
DOCKING_BOX_PADDING = 5.0
DOCKING_FIXED_BOX_SIZE = [25.0, 25.0, 25.0]
DOCKING_EXHAUSTIVENESS = 8
DOCKING_NUM_POSES = 9
DOCKING_CHECKPOINT_EVERY = 25

DOCKING_OBABEL_PATH = "/opt/homebrew/bin/obabel"
DOCKING_REMOVE_HETATM_RESNAMES = ["RAP", "HOH"]

# Notebook 06: binding-affinity prediction workflow
BINDING_TARGET_COL = "vina_score"
BINDING_SIMILARITY_COL = DOCKING_SIMILARITY_COL
BINDING_TEST_SIZE = 0.20

# Phase 6 defaults
BINDING_P6_K_RANGE = range(2, 15)
BINDING_P6_PCA_CUMVAR_THRESHOLD = 0.85
BINDING_P6_UMAP_N_NEIGHBORS = 50
BINDING_P6_UMAP_MIN_DIST = 0.5

# Phase 7 regularization defaults
# These define np.logspace(start, stop, num) used for ElasticNet/Ridge/Lasso CV.
BINDING_P7_ALPHA_LOGSPACE_START = -4
BINDING_P7_ALPHA_LOGSPACE_STOP = 2
BINDING_P7_ALPHA_LOGSPACE_NUM = 50
BINDING_P7_ELASTIC_L1_RATIOS = [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]

# Future analysis defaults
RUN_PCA = True
RUN_UMAP = True
N_CLUSTERS = 4

