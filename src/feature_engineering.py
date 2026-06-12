"""Feature engineering helpers for molecular clustering workflows."""

from typing import Callable, Dict, Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Fragments, rdMolDescriptors


def extract_features_safe(smiles: str) -> Optional[Dict[str, float]]:
	"""Extract a small, robust set of chemistry descriptors from a SMILES string."""
	mol = Chem.MolFromSmiles(smiles)
	if mol is None:
		return None

	try:
		Chem.SanitizeMol(mol)
		return {
			"hbd": rdMolDescriptors.CalcNumHBD(mol),
			"hba": rdMolDescriptors.CalcNumHBA(mol),
			"rings": rdMolDescriptors.CalcNumRings(mol),
			"aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
			"tpsa": Descriptors.TPSA(mol),
			"logp": Descriptors.MolLogP(mol),
			"mw": Descriptors.MolWt(mol),
			"amide": Fragments.fr_amide(mol),
		}
	except Exception:
		return None


def build_feature_df(
	df_input: pd.DataFrame,
	target_col: str = "sim_r1_tv_frag",
	extractor: Callable[[str], Optional[Dict[str, float]]] = extract_features_safe,
) -> pd.DataFrame:
	"""Build a clean feature dataframe with optional similarity target column."""
	rows = []

	for _, row in df_input.iterrows():
		feats = extractor(row["smiles"])
		if feats is None:
			continue

		feats.update(
			{
				"name": row.get("name", np.nan),
				"smiles": row.get("smiles", np.nan),
				target_col: row.get(target_col, np.nan),
			}
		)
		rows.append(feats)

	return pd.DataFrame(rows)
