"""Strongly connected component observables."""
from __future__ import annotations
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


def scc_decomposition(graph: csr_matrix):
    """Return vertex labels and component sizes."""
    n_components, labels = connected_components(
        graph.tocsr(), directed=True, connection="strong"
    )
    sizes = np.bincount(labels, minlength=n_components)
    return labels, sizes


def largest_scc_fraction(graph: csr_matrix) -> float:
    """Fraction of vertices in the largest SCC."""
    if graph.shape[0] == 0:
        return 0.0
    _, sizes = scc_decomposition(graph)
    return float(sizes.max() / graph.shape[0])


def largest_scc_members(graph: csr_matrix) -> np.ndarray:
    """Indices of one largest SCC, with deterministic label tie-breaking."""
    labels, sizes = scc_decomposition(graph)
    label = int(np.flatnonzero(sizes == sizes.max())[0])
    return np.flatnonzero(labels == label)


def second_largest_scc_size(graph: csr_matrix) -> int:
    _, sizes = scc_decomposition(graph)
    order = np.sort(sizes)[::-1]
    return int(order[1]) if len(order) > 1 else 0
