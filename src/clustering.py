"""Clustering helpers for unsupervised chemistry analysis."""

from typing import Iterable, List, Sequence, Tuple

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def evaluate_k(X: np.ndarray, k_range: Iterable[int] = range(2, 15)) -> Tuple[List[float], List[float]]:
	"""Return inertia and silhouette scores across a range of k values."""
	inertia: List[float] = []
	silhouette_scores: List[float] = []

	for k in k_range:
		kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
		labels = kmeans.fit_predict(X)
		inertia.append(float(kmeans.inertia_))
		silhouette_scores.append(float(silhouette_score(X, labels)))

	return inertia, silhouette_scores


def run_kmeans(X: np.ndarray, k: int, random_state: int = 42) -> np.ndarray:
	"""Run K-Means and return cluster labels."""
	kmeans = KMeans(
		n_clusters=k,
		random_state=random_state,
		n_init=10,
	)
	return kmeans.fit_predict(X)
