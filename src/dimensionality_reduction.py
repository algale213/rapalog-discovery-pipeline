"""Dimensionality reduction helpers used in notebook Phase 3 workflows."""

from typing import Any, Dict, Optional, Sequence

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

try:
	import umap  # type: ignore
except Exception:  # pragma: no cover
	umap = None


def run_pca(
	X_scaled: np.ndarray,
	model_df: pd.DataFrame,
	feature_cols: Sequence[str],
	n_components: int = 2,
	random_state: int = 42,
) -> Dict[str, Any]:
	"""Fit PCA and return transformed coordinates plus loadings."""
	pca = PCA(n_components=n_components, random_state=random_state)
	X_pca = pca.fit_transform(X_scaled)

	# Add PC1 and PC2 to model
	model_df_out = model_df.copy()
	model_df_out["PC1"] = X_pca[:, 0]
	model_df_out["PC2"] = X_pca[:, 1]

	loadings = pd.DataFrame(
		pca.components_.T,
		columns=["PC1", "PC2"],
		index=list(feature_cols),
	)

	return {
		"pca": pca,
		"X_pca": X_pca,
		"model_df": model_df_out,
		"var_explained": pca.explained_variance_ratio_,
		"loadings": loadings,
	}


def plot_variance_explained(
	X_scaled: np.ndarray,
	threshold: float = 0.85,
) -> Dict[str, np.ndarray]:
	"""Plot explained variance by principal component with cumulative curve."""
	pca_full = PCA()
	pca_full.fit(X_scaled)

	explained_var = pca_full.explained_variance_ratio_
	cumulative_var = np.cumsum(explained_var)
	components = np.arange(1, len(explained_var) + 1)

	plt.figure(figsize=(9, 5))
	bars = plt.bar(components, explained_var, color="steelblue", alpha=0.9, label="Explained variance")

	for bar, var in zip(bars, explained_var):
		plt.text(
			bar.get_x() + bar.get_width() / 2,
			bar.get_height() + 0.005,
			f"{var * 100:.1f}%",
			ha="center",
			va="bottom",
			fontsize=9,
		)

	plt.plot(components, cumulative_var, color="darkorange", marker="o", linewidth=2, label="Cumulative variance")
	plt.axhline(threshold, linestyle="--", color="crimson", label=f"{int(threshold * 100)}% threshold")
	plt.xticks(components)
	plt.xlabel("Principal Component")
	plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
	plt.ylabel("Variance Explained")
	plt.title("PCA Variance Explained by Component")
	plt.legend()
	plt.tight_layout()
	plt.show()

	return {
		"explained_var": explained_var,
		"cumulative_var": cumulative_var,
	}


def run_umap(
	X_scaled: np.ndarray,
	labels: Optional[np.ndarray] = None,
	sample_size: int = 10000,
	n_neighbors: int = 50,
	min_dist: float = 0.5,
	random_state: int = 42,
) -> Dict[str, Any]:
	"""Subsample, fit UMAP, and optionally return sampled labels."""
	if umap is None:
		raise ImportError("UMAP is not installed. Install umap-learn to use run_umap().")

	n = len(X_scaled)
	size = min(sample_size, n)
	rng = np.random.default_rng(random_state)
	idx = rng.choice(n, size=size, replace=False)
	X_umap_in = X_scaled[idx]

	umap_model = umap.UMAP(
		n_neighbors=n_neighbors,
		min_dist=min_dist,
		n_components=2,
		random_state=random_state,
	)
	X_umap = umap_model.fit_transform(X_umap_in)

	labels_subset = None if labels is None else labels[idx]

	return {
		"embedding": X_umap,
		"indices": idx,
		"labels": labels_subset,
		"model": umap_model,
	}
