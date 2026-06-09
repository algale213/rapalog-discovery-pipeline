from rdkit.Chem import AllChem


def get_fp_from_mole(mol, radius: int, n_bits: int = 2048, use_features: bool = True):
	"""Create Morgan fingerprint bit vector from RDKit molecule."""
	return AllChem.GetMorganFingerprintAsBitVect(
		mol,
		radius,
		n_bits,
		useFeatures=use_features,
	)


def get_fp_from_mol(mol, radius: int, n_bits: int = 2048, use_features: bool = True):
	"""Alias retained for compatibility with prior notebook naming."""
	return get_fp_from_mole(mol, radius, n_bits=n_bits, use_features=use_features)
