import numpy as np
import pandas as pd
import pytest
from rdkit import Chem, DataStructs

from src.clustering import run_kmeans
from src.feature_engineering import build_feature_df, extract_features_safe
from src.fingerprints import get_fp_from_mol, get_fp_from_mole
from src.similarity import compute_similarity_scores, rank_candidates


def test_extract_features_for_ethanol():
    features = extract_features_safe("CCO")

    assert features is not None
    assert features["hbd"] == 1
    assert features["hba"] == 1
    assert features["rings"] == 0
    assert features["mw"] == pytest.approx(46.069, abs=0.01)


def test_extract_features_rejects_invalid_smiles():
    assert extract_features_safe("not-a-smiles") is None


def test_build_feature_df_skips_invalid_smiles():
    source = pd.DataFrame(
        {
            "name": ["ethanol", "invalid"],
            "smiles": ["CCO", "not-a-smiles"],
            "sim_r1_tv_frag": [0.25, 0.10],
        }
    )

    result = build_feature_df(source)

    assert len(result) == 1
    assert result.iloc[0]["name"] == "ethanol"
    assert result.iloc[0]["sim_r1_tv_frag"] == pytest.approx(0.25)


def test_fingerprint_aliases_are_equivalent():
    molecule = Chem.MolFromSmiles("CCO")

    fp_original = get_fp_from_mole(molecule, radius=2, n_bits=128)
    fp_alias = get_fp_from_mol(molecule, radius=2, n_bits=128)

    assert fp_original.GetNumBits() == 128
    assert DataStructs.TanimotoSimilarity(fp_original, fp_alias) == pytest.approx(1.0)


def test_identical_molecule_has_maximum_similarity():
    molecule = Chem.MolFromSmiles("CCO")
    candidates = [("CCO", molecule, "ethanol", "test-1")]

    result = compute_similarity_scores(candidates, molecule)

    assert len(result) == 1
    assert result.iloc[0]["coconut_id"] == "test-1"
    assert result.iloc[0]["sim_tanimoto"] == pytest.approx(1.0)
    assert result.iloc[0]["sim_tv_sym"] == pytest.approx(1.0)
    assert result.iloc[0]["sim_tv_frag"] == pytest.approx(1.0)


def test_rank_candidates_orders_scores_and_adds_rank():
    source = pd.DataFrame(
        {
            "name": ["low", "high", "middle"],
            "sim_tv_frag": [0.10, 0.80, 0.40],
        }
    )

    ranked = rank_candidates(source, by="sim_tv_frag")

    assert ranked["name"].tolist() == ["high", "middle", "low"]
    assert ranked["entry_n"].tolist() == [1, 2, 3]


def test_kmeans_is_reproducible_with_fixed_seed():
    values = np.array(
        [
            [0.0, 0.0],
            [0.1, 0.1],
            [5.0, 5.0],
            [5.1, 5.1],
        ]
    )

    first = run_kmeans(values, k=2, random_state=42)
    second = run_kmeans(values, k=2, random_state=42)

    assert np.array_equal(first, second)
    assert len(set(first)) == 2