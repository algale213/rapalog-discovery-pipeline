"""Utilities for Phase 6/7 binding-affinity modeling.

This module centralizes dataset preparation, linear-model training,
regularized-model training, model comparison, feature-importance plotting,
and subgroup error analysis for the binding-affinity workflow.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.linear_model import (
	ElasticNetCV,
	LassoCV,
	LinearRegression,
	RidgeCV,
)

from src import config
from src.clustering import run_kmeans


DEFAULT_CHEMICAL_FEATURE_COLS = [
	"hbd",
	"hba",
	"rings",
	"aromatic_rings",
	"tpsa",
	"logp",
	"mw",
	"amide",
]


def _safe_group_sample_size(group_sizes: pd.Series, requested_n: int) -> int:
	"""Return a valid per-group sample size (>=1 and <= smallest group)."""
	if group_sizes.empty:
		return 0
	return max(1, min(int(requested_n), int(group_sizes.min())))


def _regression_metrics(y_true: pd.Series, y_pred: np.ndarray) -> dict[str, float]:
	"""Compute standard regression metrics."""
	rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
	mae = float(mean_absolute_error(y_true, y_pred))
	r2 = float(r2_score(y_true, y_pred))
	return {
		"r2": r2,
		"rmse": rmse,
		"mae": mae,
	}


def prepare_binding_dataset(
	docking_full_df: pd.DataFrame,
	*,
	target_col: str = "vina_score",
	similarity_col: str = config.DOCKING_SIMILARITY_COL,
	chemical_feature_cols: list[str] | None = None,
	test_size: float = 0.2,
	random_state: int = config.RANDOM_STATE,
) -> dict[str, Any]:
	"""Prepare Phase 6/7 data products from docking results.

	Returns a dictionary containing:
	- Phase 6 artifacts: feature matrix, scaled matrix, PCA/loadings, k-means labels
	- Phase 7 artifacts: modeling dataframe, stratified train/test split, scaled X/y
	"""
	if chemical_feature_cols is None:
		chemical_feature_cols = DEFAULT_CHEMICAL_FEATURE_COLS

	id_cols = ["name", "smiles"]
	optional_sampling_cols = ["cluster_k2", "quartile"]
	feature_cols = list(chemical_feature_cols) + [similarity_col]

	required_for_model = id_cols + feature_cols + [target_col]
	missing = [c for c in required_for_model if c not in docking_full_df.columns]
	if missing:
		raise ValueError(f"Missing required columns from docking_full_df: {missing}")

	model_df = docking_full_df.copy()
	model_df = model_df.dropna(subset=required_for_model).reset_index(drop=True)
	model_df["original_index"] = np.arange(len(model_df))

	# Ensure sampling columns exist for stratified split.
	if "cluster_k2" not in model_df.columns:
		X_tmp = model_df[feature_cols]
		X_tmp_scaled = StandardScaler().fit_transform(X_tmp)
		model_df["cluster_k2"] = run_kmeans(X_tmp_scaled, 2, random_state=random_state)

	if "quartile" not in model_df.columns:
		model_df["quartile"] = pd.qcut(
			model_df[similarity_col],
			q=4,
			labels=[3, 2, 1, 0],
			duplicates="drop",
		)

	# Build Phase 6 feature space.
	phase6_drop_cols = [
		"name",
		"smiles",
		"vina_score",
		"cluster_k2",
		"resid_diff",
		"quartile",
		"original_index",
	]
	phase6_drop_cols = [c for c in phase6_drop_cols if c in model_df.columns]

	X_chemsim_features = model_df.drop(columns=phase6_drop_cols).copy()
	if similarity_col in X_chemsim_features.columns and similarity_col != "similarity":
		X_chemsim_features = X_chemsim_features.rename(columns={similarity_col: "similarity"})

	scaler_phase6 = StandardScaler()
	X_chemsim_scaled = scaler_phase6.fit_transform(X_chemsim_features)

	pca_chemsim = PCA(n_components=2, random_state=random_state)
	X_chemsim_pca = pca_chemsim.fit_transform(X_chemsim_scaled)

	phase6_df = model_df.copy()
	phase6_df["PC1"] = X_chemsim_pca[:, 0]
	phase6_df["PC2"] = X_chemsim_pca[:, 1]

	loadings_chemsim = pd.DataFrame(
		pca_chemsim.components_.T,
		columns=["PC1", "PC2"],
		index=X_chemsim_features.columns,
	)

	labels_chemsim_k2 = run_kmeans(X_chemsim_scaled, 2, random_state=random_state)
	labels_chemsim_k4 = run_kmeans(X_chemsim_scaled, 4, random_state=random_state)
	phase6_df["cluster_k2_phase6"] = labels_chemsim_k2
	phase6_df["cluster_k4_phase6"] = labels_chemsim_k4

	# Build Phase 7 train/test split, preserving cluster/quartile structure.
	n_total = len(model_df)
	regression_df_c0 = model_df[model_df["cluster_k2"] == 0].copy()
	regression_df_c1 = model_df[model_df["cluster_k2"] == 1].copy()

	ratio_c0 = len(regression_df_c0) / n_total if n_total else 0.0
	ratio_c1 = len(regression_df_c1) / n_total if n_total else 0.0

	n_test_c0_per_q = int(round((n_total * ratio_c0 * test_size) / 4.0))
	n_test_c1_per_q = int(round((n_total * ratio_c1 * test_size) / 4.0))

	c0_group_sizes = regression_df_c0.groupby("quartile").size()
	c1_group_sizes = regression_df_c1.groupby("quartile").size()
	n0 = _safe_group_sample_size(c0_group_sizes, n_test_c0_per_q)
	n1 = _safe_group_sample_size(c1_group_sizes, n_test_c1_per_q)

	test_c0 = (
		regression_df_c0.groupby("quartile", group_keys=False)
		.sample(n=n0, random_state=random_state)
		if n0 > 0
		else regression_df_c0.iloc[0:0]
	)
	test_c1 = (
		regression_df_c1.groupby("quartile", group_keys=False)
		.sample(n=n1, random_state=random_state)
		if n1 > 0
		else regression_df_c1.iloc[0:0]
	)

	test_df = pd.concat([test_c0, test_c1]).reset_index(drop=True)
	train_df = model_df[
		~model_df["original_index"].isin(test_df["original_index"])
	].reset_index(drop=True)

	X_train = train_df[feature_cols].copy()
	y_train = train_df[target_col].copy()
	X_test = test_df[feature_cols].copy()
	y_test = test_df[target_col].copy()

	scaler_phase7 = StandardScaler()
	X_train_scaled = scaler_phase7.fit_transform(X_train)
	X_test_scaled = scaler_phase7.transform(X_test)

	return {
		"feature_cols": feature_cols,
		"target_col": target_col,
		"similarity_col": similarity_col,
		"model_df": model_df,
		"phase6_df": phase6_df,
		"X_chemsim_features": X_chemsim_features,
		"X_chemsim_scaled": X_chemsim_scaled,
		"X_chemsim_pca": X_chemsim_pca,
		"pca_chemsim": pca_chemsim,
		"loadings_chemsim": loadings_chemsim,
		"explained_var_chemsim": pca_chemsim.explained_variance_ratio_,
		"labels_chemsim_k2": labels_chemsim_k2,
		"labels_chemsim_k4": labels_chemsim_k4,
		"train_df": train_df,
		"test_df": test_df,
		"X_train": X_train,
		"y_train": y_train,
		"X_test": X_test,
		"y_test": y_test,
		"X_train_scaled": X_train_scaled,
		"X_test_scaled": X_test_scaled,
		"scaler_phase6": scaler_phase6,
		"scaler_phase7": scaler_phase7,
	}


def train_baseline_regression(
	X_train_scaled: np.ndarray,
	y_train: pd.Series,
	X_test_scaled: np.ndarray,
	y_test: pd.Series,
) -> dict[str, Any]:
	"""Train baseline linear regression and return predictions + metrics."""
	model = LinearRegression()
	model.fit(X_train_scaled, y_train)

	y_train_pred = model.predict(X_train_scaled)
	y_test_pred = model.predict(X_test_scaled)

	train_metrics = _regression_metrics(y_train, y_train_pred)
	test_metrics = _regression_metrics(y_test, y_test_pred)

	return {
		"name": "Baseline Linear Regression",
		"model": model,
		"y_train_pred": y_train_pred,
		"y_test_pred": y_test_pred,
		"train_metrics": train_metrics,
		"test_metrics": test_metrics,
	}


def train_elastic_net(
	X_train_scaled: np.ndarray,
	y_train: pd.Series,
	X_test_scaled: np.ndarray,
	y_test: pd.Series,
	*,
	alphas: np.ndarray | None = None,
	l1_ratios: list[float] | None = None,
	random_state: int = config.RANDOM_STATE,
) -> dict[str, Any]:
	"""Train Elastic Net with CV and return predictions + metrics."""
	if alphas is None:
		alphas = np.logspace(-4, 2, 50)
	if l1_ratios is None:
		l1_ratios = [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]

	model = ElasticNetCV(
		alphas=alphas,
		l1_ratio=l1_ratios,
		cv=5,
		max_iter=10000,
		n_jobs=-1,
		random_state=random_state,
	)
	model.fit(X_train_scaled, y_train)

	y_train_pred = model.predict(X_train_scaled)
	y_test_pred = model.predict(X_test_scaled)

	train_metrics = _regression_metrics(y_train, y_train_pred)
	test_metrics = _regression_metrics(y_test, y_test_pred)

	return {
		"name": "Elastic Net",
		"model": model,
		"best_alpha": float(model.alpha_),
		"best_l1_ratio": float(model.l1_ratio_),
		"y_train_pred": y_train_pred,
		"y_test_pred": y_test_pred,
		"train_metrics": train_metrics,
		"test_metrics": test_metrics,
	}


def train_ridge(
	X_train_scaled: np.ndarray,
	y_train: pd.Series,
	X_test_scaled: np.ndarray,
	y_test: pd.Series,
	*,
	alphas: np.ndarray | None = None,
) -> dict[str, Any]:
	"""Train RidgeCV model and return predictions + metrics."""
	if alphas is None:
		alphas = np.logspace(-4, 2, 50)

	model = RidgeCV(alphas=alphas, cv=5)
	model.fit(X_train_scaled, y_train)

	y_train_pred = model.predict(X_train_scaled)
	y_test_pred = model.predict(X_test_scaled)

	train_metrics = _regression_metrics(y_train, y_train_pred)
	test_metrics = _regression_metrics(y_test, y_test_pred)

	return {
		"name": "Ridge",
		"model": model,
		"best_alpha": float(model.alpha_),
		"y_train_pred": y_train_pred,
		"y_test_pred": y_test_pred,
		"train_metrics": train_metrics,
		"test_metrics": test_metrics,
	}


def train_lasso(
	X_train_scaled: np.ndarray,
	y_train: pd.Series,
	X_test_scaled: np.ndarray,
	y_test: pd.Series,
	*,
	alphas: np.ndarray | None = None,
	random_state: int = config.RANDOM_STATE,
) -> dict[str, Any]:
	"""Train LassoCV model and return predictions + metrics."""
	if alphas is None:
		alphas = np.logspace(-4, 2, 50)

	model = LassoCV(
		alphas=alphas,
		cv=5,
		max_iter=10000,
		random_state=random_state,
		n_jobs=-1,
	)
	model.fit(X_train_scaled, y_train)

	y_train_pred = model.predict(X_train_scaled)
	y_test_pred = model.predict(X_test_scaled)

	train_metrics = _regression_metrics(y_train, y_train_pred)
	test_metrics = _regression_metrics(y_test, y_test_pred)

	return {
		"name": "Lasso",
		"model": model,
		"best_alpha": float(model.alpha_),
		"y_train_pred": y_train_pred,
		"y_test_pred": y_test_pred,
		"train_metrics": train_metrics,
		"test_metrics": test_metrics,
	}


def compare_models(model_results: list[dict[str, Any]]) -> pd.DataFrame:
	"""Create a comparison dataframe across trained models."""
	rows: list[dict[str, Any]] = []
	for result in model_results:
		rows.append(
			{
				"model": result["name"],
				"train_R2": result["train_metrics"]["r2"],
				"test_R2": result["test_metrics"]["r2"],
				"train_RMSE": result["train_metrics"]["rmse"],
				"test_RMSE": result["test_metrics"]["rmse"],
				"train_MAE": result["train_metrics"]["mae"],
				"test_MAE": result["test_metrics"]["mae"],
			}
		)
	return pd.DataFrame(rows)


def plot_feature_importance(
	model_result: dict[str, Any],
	feature_names: list[str],
	*,
	title: str | None = None,
	figsize: tuple[int, int] = (8, 5),
	show: bool = True,
	output_path: str | Path | None = None,
) -> pd.DataFrame:
	"""Build coefficient table and optional horizontal bar plot."""
	model = model_result["model"]
	if not hasattr(model, "coef_"):
		raise ValueError("Model does not expose coef_.")

	coef_df = pd.DataFrame(
		{
			"feature": feature_names,
			"coefficient": model.coef_,
			"abs_coefficient": np.abs(model.coef_),
		}
	).sort_values("abs_coefficient", ascending=False)

	plot_df = coef_df.sort_values("coefficient")
	plt.figure(figsize=figsize)
	plt.barh(plot_df["feature"], plot_df["coefficient"])
	plt.axvline(0, linewidth=1)
	plt.xlabel("Standardized coefficient")
	plt.ylabel("Feature")
	plt.title(title or f"{model_result['name']} Feature Weights")
	plt.tight_layout()

	if output_path is not None:
		out = Path(output_path)
		out.parent.mkdir(parents=True, exist_ok=True)
		plt.savefig(out, dpi=300, bbox_inches="tight")

	if show:
		plt.show()
	else:
		plt.close()

	return coef_df


def error_analysis(
	test_df: pd.DataFrame,
	y_true: pd.Series,
	y_pred: np.ndarray,
	*,
	cluster_col: str = "cluster_k2",
	quartile_col: str = "quartile",
	similarity_col: str = config.DOCKING_SIMILARITY_COL,
) -> dict[str, pd.DataFrame]:
	"""Compute subgroup error summaries and return plotting dataframe."""
	eval_df = test_df.copy()
	eval_df["y_true"] = np.asarray(y_true)
	eval_df["y_pred"] = y_pred
	eval_df["residual"] = eval_df["y_true"] - eval_df["y_pred"]
	eval_df["abs_error"] = eval_df["residual"].abs()
	eval_df["squared_error"] = eval_df["residual"] ** 2

	group_cols = [c for c in [cluster_col, quartile_col] if c in eval_df.columns]
	if not group_cols:
		raise ValueError("No grouping columns found for error analysis.")

	summary_df = (
		eval_df.groupby(group_cols)
		.agg(
			n=("y_true", "size"),
			mean_actual_vina=("y_true", "mean"),
			mean_predicted_vina=("y_pred", "mean"),
			mean_residual=("residual", "mean"),
			MAE=("abs_error", "mean"),
			RMSE=("squared_error", lambda x: np.sqrt(x.mean())),
		)
		.reset_index()
	)

	cols = ["y_true", "y_pred", "residual", similarity_col]
	cols = [c for c in cols if c in eval_df.columns]
	residual_df = eval_df[cols].copy()

	return {
		"eval_df": eval_df,
		"summary_df": summary_df,
		"residual_df": residual_df,
	}
