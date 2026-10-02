"""Fixed-multiset spatial reassignment and fidelity observables."""
from __future__ import annotations
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.stats import ks_2samp, wasserstein_distance
from .scc import largest_scc_fraction


def smooth_periodic_field(
    points: np.ndarray,
    L: float,
    correlation_length: float,
    seed: int,
    grid_size: int = 256,
    shift=(0.0, 0.0),
):
    """Evaluate a periodic Gaussian-smoothed random field at point coordinates."""
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((grid_size, grid_size))
    grid_spacing = L / grid_size
    sigma_grid = correlation_length / np.sqrt(2.0) / grid_spacing
    field = gaussian_filter(noise, sigma=sigma_grid, mode="wrap")
    x = ((points[:, 0] + shift[0]) % L) / grid_spacing
    y = ((points[:, 1] + shift[1]) % L) / grid_spacing
    return map_coordinates(field, [y, x], order=1, mode="wrap")


def reassign_fixed_multiset(
    points: np.ndarray,
    orientations: np.ndarray,
    L: float,
    correlation_length: float,
    seed: int,
    grid_size: int = 256,
    shift=(0.0, 0.0),
):
    """Assign sorted orientations by the rank of an independent smooth field."""
    values = smooth_periodic_field(
        points, L, correlation_length, seed, grid_size=grid_size, shift=shift
    )
    order = np.argsort(values, kind="mergesort")
    reassigned = np.empty_like(np.asarray(orientations, dtype=float))
    reassigned[order] = np.sort(np.asarray(orientations, dtype=float))
    return reassigned, values


def edge_observables(graph, points, theta):
    src = np.repeat(np.arange(graph.shape[0]), np.diff(graph.indptr))
    dst = graph.indices
    delta = points[dst] - points[src]
    lengths = np.linalg.norm(delta, axis=1)
    unit = delta / lengths[:, None]
    direction = np.column_stack((np.sin(theta[src]), -np.cos(theta[src])))
    alpha = np.arctan2(delta[:, 0], -delta[:, 1]) - theta[src]
    alpha = np.arctan2(np.sin(alpha), np.cos(alpha))
    return {
        "length": lengths,
        "source_relative_angle": alpha,
        "vertical_displacement": delta[:, 1],
        "source_relative_alignment": np.einsum("ij,ij->i", unit, direction),
        "local_alignment": float(np.mean(np.einsum("ij,ij->i", unit, direction))),
        "downward_alignment": float(np.mean(-unit[:, 1])),
        "upstream_fraction": float(np.mean(delta[:, 1] > 0)),
        "S": largest_scc_fraction(graph),
    }


def marginal_fidelity(original: dict, reassigned: dict, R: float = 1.0) -> dict:
    """Canonical length/angle/displacement distances; angle is wrapped to [-pi, pi]."""
    fields = {
        "length": 2 * R,
        "source_relative_angle": np.pi,
        "vertical_displacement": 2 * R,
    }
    out = {}
    for name, scale in fields.items():
        a = np.asarray(original[name])
        b = np.asarray(reassigned[name])
        out[f"{name}_ks"] = float(ks_2samp(a, b).statistic)
        out[f"{name}_wasserstein_normalized"] = float(
            wasserstein_distance(a, b) / scale
        )
    return out
