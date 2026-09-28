"""Filtered-graph recurrence, core induction, and surface observables."""
from __future__ import annotations
import numpy as np
from scipy.sparse import csr_matrix
from .scc import scc_decomposition


def strict_core_mask(points: np.ndarray, L: float, width: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    return (
        (points[:, 0] > width)
        & (points[:, 0] < L - width)
        & (points[:, 1] > width)
        & (points[:, 1] < L - width)
    )


def induce_vertices(graph: csr_matrix, mask: np.ndarray):
    """Induce a graph on mask and return (subgraph, original indices)."""
    indices = np.flatnonzero(np.asarray(mask, dtype=bool))
    return graph[indices][:, indices].tocsr(), indices


def internal_recurrence(graph: csr_matrix, mask: np.ndarray) -> dict:
    """Participation and concentration among nontrivial SCCs in an induced core."""
    core, indices = induce_vertices(graph, mask)
    labels, sizes = scc_decomposition(core)
    nontrivial = sizes >= 2
    recurrent_mass = int(sizes[nontrivial].sum())
    largest = int(sizes[nontrivial].max()) if nontrivial.any() else 0
    n = int(core.shape[0])
    participation = recurrent_mass / n if n else 0.0
    concentration = largest / recurrent_mass if recurrent_mass else 0.0
    return {
        "n_core": n,
        "recurrent_mass": recurrent_mass,
        "largest_recurrent_component": largest,
        "f_internal_rec": float(participation),
        "eta_internal": float(concentration),
        "largest_fraction_of_core": float(largest / n if n else 0.0),
        "identity_residual": float(
            largest / n - participation * concentration if n else 0.0
        ),
        "indices": indices,
        "labels": labels,
        "sizes": sizes,
    }


def downstream_depth(points: np.ndarray, L: float, side: str = "bottom") -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if side == "bottom":
        return points[:, 1]
    if side == "top":
        return L - points[:, 1]
    if side == "left":
        return points[:, 0]
    if side == "right":
        return L - points[:, 0]
    raise ValueError("side must be bottom, top, left, or right")


def component_d90(points: np.ndarray, members: np.ndarray, L: float, side="bottom") -> float:
    depth = np.sort(downstream_depth(points, L, side)[np.asarray(members, dtype=int)])
    if not len(depth):
        return float("nan")
    return float(depth[int(np.ceil(0.9 * len(depth))) - 1])


def occupancy_profile(
    points: np.ndarray,
    members: np.ndarray,
    L: float,
    bin_edges: np.ndarray,
    side: str = "bottom",
):
    depth = downstream_depth(points, L, side)
    member_mask = np.zeros(len(points), dtype=bool)
    member_mask[np.asarray(members, dtype=int)] = True
    total, _ = np.histogram(depth, bins=bin_edges)
    occupied, _ = np.histogram(depth[member_mask], bins=bin_edges)
    profile = np.divide(
        occupied, total, out=np.zeros_like(occupied, dtype=float), where=total > 0
    )
    return profile, total
