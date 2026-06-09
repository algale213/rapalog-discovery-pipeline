from pathlib import Path
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors

import pandas as pd

from src import config


# Load Coconut Database of Natural Products
def project_root() -> Path:
    """Return the project root based on this file location."""
    return Path(__file__).resolve().parents[1]


def raw_data_dir() -> Path:
    """Return the project's raw data directory."""
    return project_root() / "data" / "raw"


def coconut_csv_path(filename: str = "coconut_csv-04-2026.csv") -> Path:
    """Return the expected path to the COCONUT CSV in data/raw."""
    return raw_data_dir() / filename


def load_coconut_database(
    filename: str = "coconut_csv-04-2026.csv",
    skip_loading: Optional[bool] = None,
    preview_rows: int = 5,
) -> Optional[pd.DataFrame]:
    """Load the COCONUT CSV from data/raw and return a DataFrame."""
    if skip_loading is None:
        skip_loading = config.SKIP_COCONUT_LOADING

    if skip_loading:
        print("Skipping COCONUT loading (SKIP_COCONUT_LOADING=True).")
        return None

    csv_path = coconut_csv_path(filename)
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Expected COCONUT file not found: {csv_path}. "
            "Place the file under data/raw."
        )

    df = pd.read_csv(csv_path)
    print(f"Loaded COCONUT data from: {csv_path}")
    print(df.columns)
    print(df.head(preview_rows))
    return df


# Function to validate RDKit converts SMILES string into a molecule object on subset of data
def verify_smiles_to_molecule(
    df: pd.DataFrame,
    smiles_column: str = "canonical_smiles",
    sample_size: int = 100,
    skip_loading: Optional[bool] = None,
) -> Optional[int]:
    """Verify RDKit can parse SMILES strings into molecule objects."""
    if skip_loading is None:
        skip_loading = config.SKIP_COCONUT_LOADING

    if skip_loading:
        print("Skipping SMILES validation (SKIP_COCONUT_LOADING=True).")
        return None

    if smiles_column not in df.columns:
        raise KeyError(f"Column '{smiles_column}' not found in dataframe.")

    valid_count = 0
    sample = df[smiles_column].head(sample_size)

    for smi in sample:
        if pd.isna(smi):
            continue
        mol = Chem.MolFromSmiles(str(smi))
        if mol is not None:
            valid_count += 1

    print(f"Valid molecules: {valid_count}/{len(sample)}")
    return valid_count


# Function to quantify amount of missing SMILES data in full dataframe
def analyze_smiles_quality(
    df: pd.DataFrame,
    smiles_column: str = "canonical_smiles",
    skip_loading: Optional[bool] = None,
) -> Optional[dict]:
    """Compute missing/invalid/valid SMILES counts across the full dataframe."""
    if skip_loading is None:
        skip_loading = config.SKIP_COCONUT_LOADING

    if skip_loading:
        print("Skipping full SMILES quality analysis (SKIP_COCONUT_LOADING=True).")
        return None

    if smiles_column not in df.columns:
        raise KeyError(f"Column '{smiles_column}' not found in dataframe.")

    total = len(df)
    missing = int(df[smiles_column].isna().sum())

    invalid = 0
    for smi in df[smiles_column].dropna():
        if Chem.MolFromSmiles(str(smi)) is None:
            invalid += 1

    valid = total - missing - invalid
    percent_valid = (valid / total * 100) if total > 0 else 0.0

    summary = {
        "total": total,
        "missing": missing,
        "invalid": invalid,
        "valid": valid,
        "percent_valid": percent_valid,
    }

    print(f"Total molecules: {summary['total']}")
    print(f"Missing SMILES: {summary['missing']}")
    print(f"Invalid SMILES: {summary['invalid']}")
    print(f"Valid molecules: {summary['valid']}")
    print(f"% Valid: {summary['percent_valid']:.2f}%")

    return summary


# Function to calculate molecular weights from SMILES
def calculate_molecular_weights(
    df: pd.DataFrame,
    smiles_column: str = "canonical_smiles",
    mw_column: str = "MW",
    skip_loading: Optional[bool] = None,
) -> Optional[pd.DataFrame]:
    """Compute RDKit molecular weights and return a dataframe with MW column."""
    if skip_loading is None:
        skip_loading = config.SKIP_COCONUT_LOADING

    if skip_loading:
        print("Skipping molecular weight calculation (SKIP_COCONUT_LOADING=True).")
        return None

    if smiles_column not in df.columns:
        raise KeyError(f"Column '{smiles_column}' not found in dataframe.")

    def safe_mw(smi):
        if pd.isna(smi):
            return None
        mol = Chem.MolFromSmiles(str(smi))
        if mol is not None:
            return Descriptors.MolWt(mol)
        return None

    print("Computing molecular weights...")

    df_with_mw = df.copy()
    df_with_mw[mw_column] = df_with_mw[smiles_column].apply(safe_mw)

    invalid_count = int(df_with_mw[mw_column].isna().sum())
    print(f"Invalid molecules (failed parsing): {invalid_count}")

    cleaned_df = df_with_mw.dropna(subset=[mw_column]).copy()
    print(f"Dataset size after cleaning: {len(cleaned_df)}")

    return cleaned_df
