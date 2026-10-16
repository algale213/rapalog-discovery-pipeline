from collections import namedtuple
import os
import pickle
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdmolops

from src import config


Interaction = namedtuple(
	"Interaction",
	["lig_serial", "lig_atom_name", "prot_res_name", "prot_res_id", "prot_atom_name"],
)


def load_fkbp12_structure(
	cif_path: str,
	structure_id: str = "FKBP12_RAPA",
	quiet: bool = True,
):
	"""Load FKBP12-rapamycin structure from CIF and return structure + model 0."""
	parser = MMCIFParser(QUIET=quiet)
	structure = parser.get_structure(structure_id, cif_path)
	model = structure[0]
	return structure, model


def identify_contact_atoms(
	structure,
	model=None,
	ligand_resname: str = "RAP",
	distance: float = 3.5,
	excluded_resnames: Sequence[str] = ("RAP", "HOH"),
):
	"""Identify ligand atoms that contact protein atoms within distance cutoff."""
	if model is None:
		model = structure[0]

	atoms = list(structure.get_atoms())
	ns = NeighborSearch(atoms)
	contact_atoms: Set = set()

	for chain in model:
		for residue in chain:
			if residue.resname == ligand_resname:
				for lig_atom in residue:
					neighbors = ns.search(lig_atom.coord, distance)
					for prot_atom in neighbors:
						if prot_atom.parent.resname not in excluded_resnames:
							contact_atoms.add(lig_atom)
							break

	return contact_atoms, ns


def build_interaction_map(
	contact_atoms: Iterable,
	neighbor_search: NeighborSearch,
	distance: float = 3.5,
	excluded_resnames: Sequence[str] = ("RAP", "HOH"),
) -> List[Interaction]:
	"""Build atom-level ligand-protein interaction entries for contact atoms."""
	interaction_map: List[Interaction] = []

	for lig_atom in contact_atoms:
		neighbors = neighbor_search.search(lig_atom.coord, distance)
		for prot_atom in neighbors:
			if prot_atom.parent.resname in excluded_resnames:
				continue
			interaction_map.append(
				Interaction(
					lig_serial=lig_atom.serial_number,
					lig_atom_name=lig_atom.name,
					prot_res_name=prot_atom.parent.resname,
					prot_res_id=prot_atom.parent.id[1],
					prot_atom_name=prot_atom.name,
				)
			)

	return interaction_map


def classify_interactions(
	interaction_map: Iterable[Interaction],
	hydrophobic_residues: Optional[Set[str]] = None,
	polar_residues: Optional[Set[str]] = None,
) -> Dict[str, List[Interaction]]:
	"""Classify interactions into hydrogen-bond, hydrophobic, and polar groups."""
	if hydrophobic_residues is None:
		hydrophobic_residues = {"PHE", "TRP", "LEU", "ILE", "VAL", "MET", "ALA", "PRO"}
	if polar_residues is None:
		polar_residues = {"SER", "THR", "ASN", "GLN", "TYR", "HIS"}

	hbond_interactions: List[Interaction] = []
	hydrophobic_interactions: List[Interaction] = []
	polar_interactions: List[Interaction] = []

	for entry in interaction_map:
		lig_atom = entry.lig_atom_name.strip()
		prot_atom = entry.prot_atom_name.strip()
		prot_res = entry.prot_res_name.strip()

		if lig_atom and prot_atom and lig_atom[0] in {"O", "N"} and prot_atom[0] in {"O", "N"}:
			hbond_interactions.append(entry)
		elif lig_atom and lig_atom[0] == "C" and prot_res in hydrophobic_residues:
			hydrophobic_interactions.append(entry)
		elif lig_atom and lig_atom[0] in {"O", "N"} and prot_res in polar_residues:
			polar_interactions.append(entry)

	return {
		"hbond_interactions": hbond_interactions,
		"hydrophobic_interactions": hydrophobic_interactions,
		"polar_interactions": polar_interactions,
	}


def extract_pharmacophore_atoms(
	contact_atoms: Iterable,
	hbond_interactions: Iterable[Interaction],
	hydrophobic_interactions: Iterable[Interaction],
	polar_interactions: Iterable[Interaction],
):
	"""Extract unique ligand atoms involved in classified interactions."""
	pharmacophore_atoms: Set = set()
	contact_atoms_list = list(contact_atoms)

	for group in [hbond_interactions, hydrophobic_interactions, polar_interactions]:
		for entry in group:
			lig_serial = entry.lig_serial
			for atom in contact_atoms_list:
				if atom.serial_number == lig_serial:
					pharmacophore_atoms.add(atom)

	return pharmacophore_atoms


def load_full_structure_rdkit(pdb_path: str, remove_hs: bool = False):
	"""Load full FKBP12-rapamycin PDB structure in RDKit."""
	return Chem.MolFromPDBFile(pdb_path, removeHs=remove_hs)


def identify_rapamycin_atom_indices(full_structure, ligand_resname: str = "RAP") -> List[int]:
	"""Identify RDKit atom indices for rapamycin atoms in full structure."""
	rapa_atom_indices: List[int] = []

	for atom in full_structure.GetAtoms():
		info = atom.GetPDBResidueInfo()
		if info is not None and info.GetResidueName().strip() == ligand_resname:
			rapa_atom_indices.append(atom.GetIdx())

	return rapa_atom_indices


def map_pharmacophore_to_rdkit_indices(
	full_structure,
	pharmacophore_atoms: Iterable,
	ligand_resname: str = "RAP",
) -> List[int]:
	"""Map pharmacophore atom names to RDKit indices in full structure."""
	rdkit_pharmacophore_indices: List[int] = []

	for atom in full_structure.GetAtoms():
		info = atom.GetPDBResidueInfo()
		if info is None:
			continue
		if info.GetResidueName().strip() != ligand_resname:
			continue
		rdkit_name = info.GetName().strip()

		for bio_atom in pharmacophore_atoms:
			if bio_atom.name.strip() == rdkit_name:
				rdkit_pharmacophore_indices.append(atom.GetIdx())
				break

	return rdkit_pharmacophore_indices


def extract_binding_fragment(
	full_structure,
	rdkit_pharmacophore_indices: Sequence[int],
	sanitize: bool = True,
):
	"""Extract and optionally sanitize rapamycin binding fragment from full structure."""
	binding_fragment = Chem.PathToSubmol(full_structure, rdkit_pharmacophore_indices)
	if sanitize:
		rdmolops.SanitizeMol(binding_fragment)
	return binding_fragment


def load_coconut_subset(csv_path: str) -> pd.DataFrame:
	"""Load COCONUT subset needed for downstream similarity workflows."""
	coconut = pd.read_csv(csv_path)
	coconut = coconut[["canonical_smiles", "name", "identifier"]].dropna(subset=["canonical_smiles"])
	coconut = coconut.rename(columns={"canonical_smiles": "smiles"})
	return coconut


def convert_smiles_to_molecules(coconut_df: pd.DataFrame) -> Tuple[List[Tuple], int]:
	"""Convert SMILES strings to RDKit molecules and return valid tuples + invalid count."""
	work_df = coconut_df.copy()
	work_df["mol"] = work_df["smiles"].apply(Chem.MolFromSmiles)
	valid = work_df[work_df["mol"].notnull()]
	molecules = list(zip(valid["smiles"], valid["mol"], valid["name"], valid["identifier"]))
	invalid_mols = len(work_df) - len(valid)
	return molecules, invalid_mols


def load_or_build_molecules(
	coconut_df: pd.DataFrame,
	pickle_path: str,
	rebuild: bool = False,
) -> Tuple[List[Tuple], int]:
	"""Load cached molecule tuples or rebuild them from SMILES and cache to disk."""
	if (not rebuild) and os.path.exists(pickle_path):
		with open(pickle_path, "rb") as f:
			molecules = pickle.load(f)
		return molecules, 0

	molecules, invalid_mols = convert_smiles_to_molecules(coconut_df)
	with open(pickle_path, "wb") as f:
		pickle.dump(molecules, f, protocol=pickle.HIGHEST_PROTOCOL)
	return molecules, invalid_mols


def prefilter_molecules(
	molecules: Sequence[Tuple[str, Chem.Mol, str, str]],
	filters: Optional[Dict[str, float]] = None,
) -> List[Tuple[str, Chem.Mol, str, str]]:
	"""Prefilter molecules by rapamycin-like chemistry property ranges."""
	if filters is None:
		filters = config.COCONUT_PREFILTERS

	filtered: List[Tuple[str, Chem.Mol, str, str]] = []

	for smi, mol, name, cid in molecules:
		try:
			mw = Descriptors.MolWt(mol)
			hba = Descriptors.NumHAcceptors(mol)
			hbd = Descriptors.NumHDonors(mol)
			logp = Crippen.MolLogP(mol)
			rings = Descriptors.RingCount(mol)
			tpsa = Descriptors.TPSA(mol)

			if not (filters["mw_min"] < mw < filters["mw_max"]):
				continue
			if not (filters["hba_min"] <= hba <= filters["hba_max"]):
				continue
			if not (filters["hbd_min"] <= hbd <= filters["hbd_max"]):
				continue
			if not (filters["logp_min"] < logp < filters["logp_max"]):
				continue
			if not (filters["rings_min"] <= rings <= filters["rings_max"]):
				continue
			if tpsa > filters["tpsa_max"]:
				continue

			filtered.append((smi, mol, name, cid))
		except Exception:
			continue

	return filtered
