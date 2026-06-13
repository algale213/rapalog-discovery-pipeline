"""Utilities for Phase 4 similarity prediction modeling.

This module provides a clean API for preparing a similarity regression dataset,
training a baseline linear model, evaluating performance, and saving diagnostic
plots/feature-importance charts.
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

from src import config
from src.feature_engineering import build_feature_df


DEFAULT_FEATURE_COLS = [
	"hbd",
	"hba",
	"rings",
	"aromatic_rings",
	"tpsa",
	"logp",
	"mw",
	"amide",
]


def _project_root() -> Path:
	"""Return project root based on this module's location."""
	return Path(__file__).resolve().parents[1]


def _results_dir_for_notebook04(results_dir: str | Path | None = None) -> Path:
	"""Resolve and create a results directory for notebook 04 outputs."""
	if results_dir is None:
		resolved = _project_root() / "results" / "notebook_04"
	else:
		resolved = Path(results_dir)
	resolved.mkdir(parents=True, exist_ok=True)
	return resolved


def _build_phase4_model_df(
	target_col: str = "sim_r1_tv_frag",
	feature_cols: list[str] | None = None,
) -> pd.DataFrame:
	"""Construct the Phase 4 modeling dataframe from cached similarity data."""
	if feature_cols is None:
		feature_cols = DEFAULT_FEATURE_COLS

	raw_dir = _project_root() / "data" / "raw"
	similarity_cache_path = raw_dir / config.SIMILARITY_SCORES_CACHE_FILE
	if not similarity_cache_path.exists():
		raise FileNotFoundError(f"Missing similarity cache: {similarity_cache_path}")

	with open(similarity_cache_path, "rb") as f:
		cached_scores = pickle.load(f)

	results_df = (
		cached_scores
		if isinstance(cached_scores, pd.DataFrame)
		else pd.DataFrame(cached_scores)
	)
	results_df = results_df.drop_duplicates(subset=["smiles"]).copy()
	if "name" in results_df.columns and "coconut_id" in results_df.columns:
		results_df["name"] = results_df["name"].fillna(results_df["coconut_id"])

	ranked = (
		results_df.sort_values(by="sim_tv_frag", ascending=False)
		.reset_index(drop=True)
		.copy()
	)
	ranked["entry_n"] = np.arange(1, len(ranked) + 1)

	feat_all = build_feature_df(ranked, target_col=target_col)
	model_df = feat_all[["name", "smiles"] + feature_cols + [target_col]].dropna().copy()
	return model_df


def prepare_similarity_dataset(
	model_df: pd.DataFrame | None = None,
	*,
	target_col: str = "sim_r1_tv_frag",
	test_size: float = 0.2,
	random_state: int = config.RANDOM_STATE,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series]:
	"""Prepare scaled train/test data for similarity-score linear regression.

	Returns:
		X_train_scaled, X_test_scaled, y_train, y_test
	"""
	if model_df is None:
		model_df = _build_phase4_model_df(target_col=target_col)

	drop_cols = ["name", "smiles", target_col, "PC1", "PC2"]
	drop_cols = [c for c in drop_cols if c in model_df.columns]
	X = model_df.drop(columns=drop_cols).copy()
	y = model_df[target_col].copy()

	X_train, X_test, y_train, y_test = train_test_split(
		X,
		y,
		test_size=test_size,
		random_state=random_state,
	)

	scaler = StandardScaler()
	X_train_scaled = pd.DataFrame(
		scaler.fit_transform(X_train),
		columns=X_train.columns,
		index=X_train.index,
	)
	X_test_scaled = pd.DataFrame(
		scaler.transform(X_test),
		columns=X_test.columns,
		index=X_test.index,
	)

	return X_train_scaled, X_test_scaled, y_train, y_test


def train_similarity_model(
	X_train: pd.DataFrame,
	y_train: pd.Series,
) -> LinearRegression:
	"""Train a baseline linear regression model for similarity prediction."""
	model = LinearRegression()
	model.fit(X_train, y_train)
	return model


def evaluate_similarity_model(
	model: LinearRegression,
	X_test: pd.DataFrame,
	y_test: pd.Series,
	*,
	results_dir: str | Path | None = None,
	metrics_filename: str = "similarity_model_metrics.csv",
) -> dict[str, Any]:
	"""Evaluate similarity model and persist metrics to CSV."""
	y_pred = model.predict(X_test)

	mse = mean_squared_error(y_test, y_pred)
	rmse = float(np.sqrt(mse))
	mae = mean_absolute_error(y_test, y_pred)
	r2 = r2_score(y_test, y_pred)
	score_range = float(y_test.max() - y_test.min())
	nmae_range = float(mae / score_range) if score_range > 0 else np.nan

	metrics = {
		"r2": float(r2),
		"mse": float(mse),
		"rmse": rmse,
		"mae": float(mae),
		"range_normalized_mae": nmae_range,
		"range_normalized_mae_pct": float(nmae_range * 100) if np.isfinite(nmae_range) else np.nan,
		"n_test": int(len(y_test)),
		"n_unique_y_test": int(pd.Series(y_test).nunique()),
	}

	out_dir = _results_dir_for_notebook04(results_dir)
	metrics_df = pd.DataFrame([metrics])
	metrics_df.to_csv(out_dir / metrics_filename, index=False)

	metrics["y_pred"] = y_pred
	return metrics


def plot_similarity_diagnostics(
	y_test: pd.Series,
	y_pred: np.ndarray,
	*,
	results_dir: str | Path | None = None,
) -> None:
	"""Save diagnostic plots for similarity regression predictions."""
	out_dir = _results_dir_for_notebook04(results_dir)

	# Predicted vs actual scatter
	plt.figure(figsize=(7, 5))
	plt.scatter(y_test, y_pred, edgecolors="none", alpha=0.1)
	plt.xlabel("Actual")
	plt.ylabel("Predicted")
	plt.title("Predicted vs Actual Similarity Score")
	plt.tight_layout()
	plt.savefig(out_dir / "similarity_pred_vs_actual.png", dpi=300, bbox_inches="tight")
	plt.show()

	# Distribution of actual similarity scores
	plt.figure(figsize=(8, 5))
	plt.hist(y_test, bins=100, alpha=0.8)
	plt.xlabel("Actual Similarity Score")
	plt.ylabel("Frequency")
	plt.title("Distribution of Tversky Similarity Scores")
	plt.tight_layout()
	plt.savefig(out_dir / "similarity_actual_distribution.png", dpi=300, bbox_inches="tight")
	plt.show()

	# Residuals histogram
	residuals = y_test - y_pred
	plt.figure(figsize=(7, 5))
	plt.hist(residuals, bins=50)
	plt.title("Residuals")
	plt.xlabel("Residual")
	plt.ylabel("Frequency")
	plt.tight_layout()
	plt.savefig(out_dir / "similarity_residual_histogram.png", dpi=300, bbox_inches="tight")
	plt.show()


def feature_importance_similarity(
	model: LinearRegression,
	feature_names: list[str],
	*,
	results_dir: str | Path | None = None,
	coefficients_filename: str = "similarity_model_coefficients.csv",
) -> pd.DataFrame:
	"""Save coefficient table + feature-weight bar chart for similarity model."""
	out_dir = _results_dir_for_notebook04(results_dir)

	coef_df = pd.DataFrame(
		{
			"feature": feature_names,
			"coefficient": model.coef_,
			"abs_coefficient": np.abs(model.coef_),
		}
	).sort_values("abs_coefficient", ascending=False)

	coef_df.to_csv(out_dir / coefficients_filename, index=False)

	coef_plot_df = coef_df.sort_values("coefficient")
	plt.figure(figsize=(8, 5))
	plt.barh(coef_plot_df["feature"], coef_plot_df["coefficient"])
	plt.axvline(0, linewidth=1)
	plt.xlabel("Standardized coefficient")
	plt.ylabel("Feature")
	plt.title("Baseline Linear Regression Feature Weights")
	plt.tight_layout()
	plt.savefig(out_dir / "similarity_feature_weights.png", dpi=300, bbox_inches="tight")
	plt.show()

	return coef_df
