import os
import pickle
from typing import Dict, Iterable, Optional, Sequence, Tuple

import pandas as pd
from rdkit import DataStructs

from src.fingerprints import get_fp_from_mole


def compute_similarity_scores(
	molecules: Sequence[Tuple[str, object, str, str]],
	binding_fragment,
) -> pd.DataFrame:
	"""Compute Tanimoto and Tversky similarity scores for candidate molecules."""
	binding_fp_r1 = get_fp_from_mole(binding_fragment, 1)
	binding_fp_r2 = get_fp_from_mole(binding_fragment, 2)

	results = []
	for smi, mol, name, cid in molecules:
		try:
			fp_r1 = get_fp_from_mole(mol, 1)
			fp_r2 = get_fp_from_mole(mol, 2)

			tan_r1 = DataStructs.TanimotoSimilarity(binding_fp_r1, fp_r1)
			tan_r2 = DataStructs.TanimotoSimilarity(binding_fp_r2, fp_r2)
			sim_tanimoto = max(tan_r1, tan_r2)

			tv_sym_r1 = DataStructs.TverskySimilarity(binding_fp_r1, fp_r1, 0.5, 0.5)
			tv_sym_r2 = DataStructs.TverskySimilarity(binding_fp_r2, fp_r2, 0.5, 0.5)
			sim_tv_sym = max(tv_sym_r1, tv_sym_r2)

			tv_frag_r1 = DataStructs.TverskySimilarity(binding_fp_r1, fp_r1, 0.1, 0.9)
			tv_frag_r2 = DataStructs.TverskySimilarity(binding_fp_r2, fp_r2, 0.1, 0.9)
			sim_tv_frag = max(tv_frag_r1, tv_frag_r2)

			results.append(
				{
					"smiles": smi,
					"name": name,
					"coconut_id": cid,
					"sim_tanimoto": sim_tanimoto,
					"sim_r1_tan": tan_r1,
					"sim_r2_tan": tan_r2,
					"sim_tv_sym": sim_tv_sym,
					"sim_r1_tv_sym": tv_sym_r1,
					"sim_r2_tv_sym": tv_sym_r2,
					"sim_tv_frag": sim_tv_frag,
					"sim_r1_tv_frag": tv_frag_r1,
					"sim_r2_tv_frag": tv_frag_r2,
				}
			)
		except Exception:
			continue

	return pd.DataFrame(results)


def rank_candidates(
	results_df: pd.DataFrame,
	by: str = "sim_tv_frag",
	ascending: bool = False,
	top_n: Optional[int] = None,
) -> pd.DataFrame:
	"""Rank candidates by selected similarity column and add entry index."""
	ranked = results_df.sort_values(by=by, ascending=ascending).reset_index(drop=True).copy()
	ranked["entry_n"] = range(1, len(ranked) + 1)
	if top_n is not None:
		return ranked.head(top_n).copy()
	return ranked


def compare_similarity_metrics(results_df: pd.DataFrame) -> Dict[str, float]:
	"""Return metric comparison summary (correlation and score ranges)."""
	corr_all = results_df["sim_r1_tan"].corr(results_df["sim_r1_tv_frag"])
	top_5000 = rank_candidates(results_df, by="sim_tv_frag", top_n=5000)
	corr_top = top_5000["sim_r1_tan"].corr(top_5000["sim_r1_tv_frag"])

	return {
		"corr_all": float(corr_all),
		"corr_top_5000": float(corr_top),
		"max_tanimoto": float(results_df["sim_tanimoto"].max()),
		"max_tversky_fragment": float(results_df["sim_tv_frag"].max()),
	}


def load_or_compute_similarity_scores(
	molecules: Sequence[Tuple[str, object, str, str]],
	binding_fragment,
	pickle_path: str,
	rebuild: bool = False,
) -> pd.DataFrame:
	"""Load cached similarity scores or compute and cache them."""
	if (not rebuild) and os.path.exists(pickle_path):
		with open(pickle_path, "rb") as f:
			cached = pickle.load(f)
		if isinstance(cached, pd.DataFrame):
			return cached
		return pd.DataFrame(cached)

	results_df = compute_similarity_scores(molecules, binding_fragment)
	with open(pickle_path, "wb") as f:
		pickle.dump(results_df, f, protocol=pickle.HIGHEST_PROTOCOL)
	return results_df
