"""Utilities for Phase 5 docking (AutoDock Vina workflow)."""

from __future__ import annotations

import pickle
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src import config


def _project_root() -> Path:
	"""Return project root based on this module location."""
	return Path(__file__).resolve().parents[1]


def _resolve_path(path_like: str | Path, base_dir: str | Path | None = None) -> Path:
	"""Resolve a path relative to base_dir (or project root)."""
	path = Path(path_like)
	if path.is_absolute():
		return path
	base = Path(base_dir) if base_dir is not None else _project_root()
	return base / path


def define_binding_box(
	pdb_file: str | Path = config.DOCKING_PDB_FILE,
	*,
	ligand_resname: str = "RAP",
	box_padding: float = config.DOCKING_BOX_PADDING,
	plot: bool = True,
	plot_output_path: str | Path | None = None,
) -> dict[str, np.ndarray]:
	"""Extract ligand coordinates from PDB and compute docking box center/size."""
	try:
		from Bio.PDB import PDBParser
	except ImportError as exc:
		raise ImportError("Biopython is required for define_binding_box().") from exc

	pdb_path = _resolve_path(pdb_file)
	if not pdb_path.exists():
		raise FileNotFoundError(f"PDB file not found: {pdb_path}")

	parser = PDBParser(QUIET=True)
	structure = parser.get_structure("receptor", str(pdb_path))

	ligand_coords_list: list[np.ndarray] = []
	for model in structure:
		for chain in model:
			for residue in chain:
				if residue.get_resname().strip() == ligand_resname:
					for atom in residue:
						ligand_coords_list.append(atom.coord)

	if not ligand_coords_list:
		raise ValueError(
			f"No atoms found for ligand residue '{ligand_resname}' in {pdb_path.name}"
		)

	ligand_coords = np.array(ligand_coords_list)
	center_coords = ligand_coords.mean(axis=0)
	min_coords = ligand_coords.min(axis=0)
	max_coords = ligand_coords.max(axis=0)
	box_size_coords = (max_coords - min_coords) + float(box_padding)

	box_min = center_coords - box_size_coords / 2
	box_max = center_coords + box_size_coords / 2
	contains_min = bool((min_coords >= box_min).all())
	contains_max = bool((max_coords <= box_max).all())

	if plot:
		fig = plt.figure(figsize=(8, 6))
		ax = fig.add_subplot(111, projection="3d")
		ax.scatter(ligand_coords[:, 0], ligand_coords[:, 1], ligand_coords[:, 2], s=10)
		ax.scatter(center_coords[0], center_coords[1], center_coords[2], color="red", s=100)

		ax.set_title("Ligand + Docking Center", pad=18)
		ax.set_xlabel("x-dimension (A)", labelpad=10)
		ax.set_ylabel("y-dimension (A)", labelpad=10)
		ax.set_zlabel("")
		ax.text2D(-0.02, 0.55, "z-dimension (A)", transform=ax.transAxes, rotation=90, va="center")
		ax.view_init(elev=22, azim=55)
		ax.set_box_aspect((1.1, 1.1, 0.9))
		fig.subplots_adjust(left=0.10, right=0.92, bottom=0.08, top=0.90)

		if plot_output_path is not None:
			output_path = _resolve_path(plot_output_path)
			output_path.parent.mkdir(parents=True, exist_ok=True)
			fig.savefig(output_path, dpi=300, bbox_inches="tight")
		plt.show()

	return {
		"ligand_coords": ligand_coords,
		"center_coords": center_coords,
		"min_coords": min_coords,
		"max_coords": max_coords,
		"box_size_coords": box_size_coords,
		"box_min": box_min,
		"box_max": box_max,
		"contains_min": contains_min,
		"contains_max": contains_max,
	}


def prepare_ligand_pdbqt(
	docking_candidates_df: pd.DataFrame,
	*,
	ligands_pickle_path: str | Path = config.DOCKING_LIGANDS_PICKLE_FILE,
	recompute: bool = config.DOCKING_RECOMPUTE_LIGANDS,
	random_seed: int = config.RANDOM_STATE,
	obabel_path: str | Path = config.DOCKING_OBABEL_PATH,
	progress_every: int = config.DOCKING_LIGAND_PREP_PROGRESS_EVERY,
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
	"""Generate or load ligand PDBQT strings from candidate SMILES."""
	if "name" not in docking_candidates_df.columns or "smiles" not in docking_candidates_df.columns:
		raise KeyError("docking_candidates_df must contain 'name' and 'smiles' columns.")

	try:
		from rdkit import Chem
	except ImportError as exc:
		raise ImportError("RDKit is required for prepare_ligand_pdbqt().") from exc

	ligands_path = _resolve_path(ligands_pickle_path)
	ligands_path.parent.mkdir(parents=True, exist_ok=True)

	def _single_smiles_to_pdbqt_with_obabel(smiles: str) -> str:
		obabel = _resolve_path(obabel_path)
		if not obabel.exists():
			raise FileNotFoundError(f"Open Babel executable not found: {obabel}")

		cmd = [
			str(obabel),
			f"-:{smiles}",
			"--gen3d",
			"-h",
			"-opdbqt",
		]
		proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
		if proc.returncode != 0:
			raise RuntimeError(
				"Open Babel ligand conversion failed. "
				f"stderr: {proc.stderr.strip() or '<empty>'}"
			)

		pdbqt_text = proc.stdout
		if not pdbqt_text.strip():
			raise RuntimeError("Open Babel returned empty ligand PDBQT output.")
		return pdbqt_text

	def _single_smiles_to_pdbqt(smiles: str) -> str:
		mol = Chem.MolFromSmiles(smiles)
		if mol is None:
			raise ValueError("Invalid SMILES")
		return _single_smiles_to_pdbqt_with_obabel(smiles)

	if recompute:
		ligand_pdbqt_data: list[tuple[str, str]] = []
		failed: list[tuple[str, str]] = []
		total = len(docking_candidates_df)

		for i, (_, row) in enumerate(docking_candidates_df.iterrows(), start=1):
			name = str(row["name"])
			smiles = str(row["smiles"])
			try:
				pdbqt_str = _single_smiles_to_pdbqt(smiles)
				ligand_pdbqt_data.append((name, pdbqt_str))
			except Exception as exc:  # noqa: BLE001 - preserve notebook behavior
				failed.append((name, str(exc)))

			if progress_every > 0 and (i % progress_every == 0 or i == total):
				print(
					f"Ligand prep progress: {i}/{total} | "
					f"success={len(ligand_pdbqt_data)} | failed={len(failed)}"
				)

		with open(ligands_path, "wb") as f:
			pickle.dump(ligand_pdbqt_data, f, protocol=pickle.HIGHEST_PROTOCOL)
		return ligand_pdbqt_data, failed

	if not ligands_path.exists():
		raise FileNotFoundError(
			f"{ligands_path} not found. Set recompute=True to create it."
		)

	with open(ligands_path, "rb") as f:
		ligand_pdbqt_data = pickle.load(f)

	return ligand_pdbqt_data, []


def prepare_receptor(
	*,
	input_pdb: str | Path = config.DOCKING_PDB_FILE,
	clean_pdb: str | Path = config.DOCKING_CLEAN_PDB_FILE,
	output_pdbqt: str | Path = config.DOCKING_RECEPTOR_PDBQT_FILE,
	obabel_path: str | Path = config.DOCKING_OBABEL_PATH,
	remove_resnames: list[str] | tuple[str, ...] = tuple(config.DOCKING_REMOVE_HETATM_RESNAMES),
) -> dict[str, Any]:
	"""Clean receptor PDB and convert to PDBQT via Open Babel."""
	input_path = _resolve_path(input_pdb)
	clean_path = _resolve_path(clean_pdb)
	output_path = _resolve_path(output_pdbqt)
	obabel = _resolve_path(obabel_path)

	if not input_path.exists():
		raise FileNotFoundError(f"Input PDB file not found: {input_path}")
	if not obabel.exists():
		raise FileNotFoundError(f"Open Babel executable not found: {obabel}")

	with open(input_path) as f:
		lines = f.readlines()

	remove_set = {r.upper() for r in remove_resnames}
	clean_lines: list[str] = []
	removed_count = 0

	for line in lines:
		if line.startswith("HETATM"):
			resname = line[17:20].strip().upper()
			if resname in remove_set:
				removed_count += 1
				continue
		clean_lines.append(line)

	clean_path.parent.mkdir(parents=True, exist_ok=True)
	output_path.parent.mkdir(parents=True, exist_ok=True)
	with open(clean_path, "w") as f:
		f.writelines(clean_lines)

	cmd = [
		str(obabel),
		str(clean_path),
		"-O",
		str(output_path),
		"-xr",
		"-xh",
	]

	proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
	if proc.returncode != 0:
		raise RuntimeError(
			"Open Babel conversion failed. "
			f"stderr: {proc.stderr.strip() or '<empty>'}"
		)

	first_line = ""
	if output_path.exists():
		with open(output_path) as f:
			first_line = f.readline().strip()

	return {
		"clean_pdb": clean_path,
		"receptor_pdbqt": output_path,
		"removed_hetatm_count": removed_count,
		"stdout": proc.stdout,
		"stderr": proc.stderr,
		"validation_first_line": first_line,
	}


def _dock_one_via_subprocess(
	receptor_path: str,
	center_coords: list,
	box_size: list,
	exhaustiveness: int,
	n_poses: int,
	pdbqt_str: str,
	timeout: int = 300,
) -> "float | None":
	"""Dock one ligand in an isolated subprocess.

	Each ligand runs in its own fresh Python process.  A C-level segfault or
	memory corruption in the Vina extension kills only that process; the
	Jupyter kernel keeps running and the next ligand is attempted normally.
	The subprocess inherits stdout/stderr so Vina's progress output still
	appears in the notebook cell.
	"""
	import base64
	import json
	import sys as _sys
	import tempfile

	pdbqt_b64 = base64.b64encode(pdbqt_str.encode("utf-8")).decode("ascii")
	with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as tf:
		result_file = tf.name

	script = "\n".join([
		"import base64, json",
		"from vina import Vina",
		f"pdbqt_str = base64.b64decode({pdbqt_b64!r}.encode()).decode('utf-8')",
		f"receptor_path = {receptor_path!r}",
		f"center_coords = {center_coords!r}",
		f"box_size = {box_size!r}",
		f"exhaustiveness = {exhaustiveness!r}",
		f"n_poses = {n_poses!r}",
		f"result_file = {result_file!r}",
		"v = Vina(sf_name='vina')",
		"v.set_receptor(receptor_path)",
		"v.compute_vina_maps(center=center_coords, box_size=box_size)",
		"v.set_ligand_from_string(pdbqt_str)",
		"v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)",
		"score = float(v.score()[0])",
		"with open(result_file, 'w') as _f: json.dump({'score': score}, _f)",
	])

	try:
		subprocess.run(
			[_sys.executable, "-c", script],
			timeout=timeout,
			check=False,
		)
	except subprocess.TimeoutExpired:
		Path(result_file).unlink(missing_ok=True)
		return None

	try:
		with open(result_file) as _f:
			return float(json.load(_f)["score"])
	except Exception:  # noqa: BLE001
		return None
	finally:
		Path(result_file).unlink(missing_ok=True)


def run_vina(
	ligand_pdbqt_data: list[tuple[str, str]],
	center_coords: np.ndarray | list[float],
	*,
	receptor_pdbqt: str | Path = config.DOCKING_RECEPTOR_PDBQT_FILE,
	checkpoint_file: str | Path = config.DOCKING_CHECKPOINT_FILE,
	skip_docking: bool = config.DOCKING_SKIP_RUN,
	exhaustiveness: int = config.DOCKING_EXHAUSTIVENESS,
	n_poses: int = config.DOCKING_NUM_POSES,
	box_size: list[float] | tuple[float, float, float] = tuple(config.DOCKING_FIXED_BOX_SIZE),
	checkpoint_every: int = config.DOCKING_CHECKPOINT_EVERY,
	timeout_per_ligand: int = 300,
) -> list[tuple[str, float]]:
	"""Run Vina docking; each ligand is isolated in its own subprocess to prevent kernel crashes."""
	receptor_path = _resolve_path(receptor_pdbqt)
	checkpoint_path = _resolve_path(checkpoint_file)

	if skip_docking:
		if checkpoint_path.exists():
			with open(checkpoint_path, "rb") as f:
				return pickle.load(f)
		return []

	try:
		from vina import Vina  # noqa: F401 – verify vina is importable before the loop
	except ImportError as exc:
		raise ImportError("vina package is required for run_vina().") from exc

	if not receptor_path.exists():
		raise FileNotFoundError(f"Receptor PDBQT not found: {receptor_path}")

	# Resume from an existing checkpoint so prior successful results are not lost.
	docking_output: list[tuple[str, float]] = []
	if checkpoint_path.exists():
		with open(checkpoint_path, "rb") as f:
			docking_output = pickle.load(f)
		already_docked = {name for name, _ in docking_output}
		remaining = [(n, p) for n, p in ligand_pdbqt_data if n not in already_docked]
		if already_docked:
			print(
				f"Resuming from checkpoint: {len(docking_output)} already docked, "
				f"{len(remaining)} remaining."
			)
	else:
		remaining = list(ligand_pdbqt_data)

	center_list: list = np.array(center_coords).tolist()
	box_list: list = list(box_size)

	for i, (name, pdbqt_str) in enumerate(remaining):
		score = _dock_one_via_subprocess(
			receptor_path=str(receptor_path),
			center_coords=center_list,
			box_size=box_list,
			exhaustiveness=int(exhaustiveness),
			n_poses=int(n_poses),
			pdbqt_str=pdbqt_str,
			timeout=timeout_per_ligand,
		)
		if score is not None:
			docking_output.append((name, score))
		else:
			print(f"Ligand {name}: subprocess failed or timed out, skipping.")

		processed = len(docking_output)
		if checkpoint_every > 0 and processed > 0 and processed % checkpoint_every == 0:
			checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
			with open(checkpoint_path, "wb") as f:
				pickle.dump(docking_output, f)
			print(f"Checkpoint: {processed} docked; {i + 1}/{len(remaining)} processed.")

	checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
	with open(checkpoint_path, "wb") as f:
		pickle.dump(docking_output, f)

	return docking_output


def analyze_docking_results(
	docking_output: list[tuple[str, float]],
	docking_candidates_df: pd.DataFrame,
	*,
	results_csv_path: str | Path | None = None,
) -> dict[str, pd.DataFrame | int]:
	"""Build merged docking results table and optional CSV export."""
	docking_results_df = pd.DataFrame(docking_output, columns=["name", "vina_score"])
	if len(docking_results_df) > 0:
		docking_results_df = docking_results_df.sort_values("vina_score").reset_index(drop=True)

	docking_full_df = docking_results_df.merge(docking_candidates_df, on="name", how="inner")
	dups = int(docking_candidates_df["name"].duplicated().sum()) if "name" in docking_candidates_df.columns else 0

	if results_csv_path is not None:
		csv_path = _resolve_path(results_csv_path)
		csv_path.parent.mkdir(parents=True, exist_ok=True)
		docking_full_df.to_csv(csv_path, index=False)

	return {
		"docking_results_df": docking_results_df,
		"docking_full_df": docking_full_df,
		"candidate_duplicate_names": dups,
	}
