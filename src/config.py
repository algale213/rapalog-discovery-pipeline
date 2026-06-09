"""Centralized runtime flags and constants for the rapalog pipeline."""

# Global project flags
SKIP_COCONUT_LOADING = False

# Phase controls
RUN_EDA = True
RUN_SIMILARITY = True
RUN_CLUSTERING = True
RUN_SIMILARITY_MODEL = True
RUN_DOCKING = True
RUN_BINDING_MODEL = True

# Notebook 02 (fragment + similarity) controls
RUN_FRAGMENT_EXTRACTION = True
RUN_COCONUT_PREP = True
RUN_FULL_SIMILARITY = True

# Cache controls
REBUILD_MOLECULES = False
REBUILD_MOLECULES_FILTERED = False
REBUILD_SIMILARITY_SCORES = False

# Cache file names (resolved relative to project root)
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

# Similarity/ranking defaults
TOP_N_CANDIDATES = 5000
RANDOM_STATE = 42

# Future analysis defaults
RUN_PCA = True
RUN_UMAP = True
N_CLUSTERS = 4

