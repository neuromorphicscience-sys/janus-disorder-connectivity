"""Fixed-N selection kernel and clean-boundary feedback predictions."""
from __future__ import annotations
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import ndtr
from scipy.stats import binom


def circular_cap_area(projection: np.ndarray, radius: float) -> np.ndarray:
    s = np.asarray(projection, dtype=float)
    u = np.clip(s / radius, -1.0, 1.0)
    return (
        radius * radius * np.arccos(u)
        - s * np.sqrt(np.maximum(0.0, radius * radius - s * s))
    )


def upstream_probability(alpha: np.ndarray, W: float) -> np.ndarray:
    alpha = np.asarray(alpha, dtype=float)
    if W == 0:
        return (np.cos(alpha) < 0).astype(float)
    answer = np.zeros_like(alpha)
    for winding in range(-5, 6):
        answer += ndtr((1.5 * np.pi - alpha + 2 * np.pi * winding) / W)
        answer -= ndtr((0.5 * np.pi - alpha + 2 * np.pi * winding) / W)
    return answer


def fixed_n_bulk_predictions(
    N: int,
    L: float,
    q: int,
    candidate_radius: float,
    W_values=(0.805, 0.815, 0.86, 0.90),
    n_radial: int = 384,
    n_transverse: int = 96,
):
    """Integrate the exact fixed-N single-edge rank kernel in a bulk disk."""
    z, wz = leggauss(n_radial)
    v, wv = leggauss(n_transverse)
    beta = (z + 1) * np.pi / 2
    wb = wz * np.pi / 2
    s = candidate_radius * np.cos(beta)
    halfwidth = candidate_radius * np.sin(beta)
    transverse = halfwidth[:, None] * v[None, :]
    alpha = np.arctan2(transverse, s[:, None])
    rank_kernel = binom.cdf(
        q - 1,
        N - 2,
        circular_cap_area(s, candidate_radius) / (L * L),
    )
    weights = (
        ((N - 1) / (L * L))
        * (wb * halfwidth * rank_kernel)[:, None]
        * (halfwidth[:, None] * wv[None, :])
    )
    degree = float(weights.sum())
    alignment = float((weights * np.cos(alpha)).sum() / degree)
    rows = []
    for W in W_values:
        rows.append(
            {
                "W": float(W),
                "N": int(N),
                "L": float(L),
                "mean_outdegree": degree,
                "local_alignment_fixedN": alignment,
                "global_drift_fixedN": alignment * np.exp(-W * W / 2),
                "bulk_up_fraction_fixedN": float(
                    np.sum(weights * upstream_probability(alpha, W)) / degree
                ),
            }
        )
    return rows


def _halfwidth_integral(limit: float, depth: float, radius: float) -> float:
    limit = min(max(limit, 0.0), radius)
    depth = min(max(depth, 0.0), radius)
    z0 = np.sqrt(max(0.0, radius * radius - limit * limit))
    if depth <= z0:
        return limit * depth

    def fint(d):
        d = np.clip(d, 0, radius)
        return 0.5 * (
            d * np.sqrt(np.maximum(0, radius * radius - d * d))
            + radius * radius * np.arcsin(d / radius)
        )

    return limit * z0 + fint(depth) - fint(z0)


def downstream_candidate_area(x: float, depth: float, L: float, radius: float) -> float:
    return _halfwidth_integral(x, depth, radius) + _halfwidth_integral(
        L - x, depth, radius
    )


def _fixed_n_deficit(area, N, L, q):
    j = np.arange(q)
    return np.sum(
        (q - j) * binom.pmf(j, N - 1, np.asarray(area)[..., None] / (L * L)),
        axis=-1,
    )


def clean_boundary_depth_profile(
    depths: np.ndarray,
    N: int,
    L: float,
    q: int,
    candidate_radius: float,
    lateral_nodes: int = 64,
):
    """Average selected upstream edges per source over lateral position."""
    xg, xw = leggauss(lateral_nodes)
    xs = candidate_radius * (xg + 1) / 2
    weights = xw * candidate_radius / 2
    values = []
    for depth in np.asarray(depths, dtype=float):
        area = np.array(
            [
                downstream_candidate_area(x, depth, L, candidate_radius)
                for x in xs
            ]
        )
        half = np.array(
            [
                downstream_candidate_area(x, candidate_radius, L, candidate_radius)
                for x in xs
            ]
        )
        side = _fixed_n_deficit(area, N, L, q) - _fixed_n_deficit(
            area + half, N, L, q
        )
        center_area = downstream_candidate_area(
            L / 2, depth, L, candidate_radius
        )
        half_center = np.pi * candidate_radius**2 / 2
        center = _fixed_n_deficit(center_area, N, L, q) - _fixed_n_deficit(
            center_area + half_center, N, L, q
        )
        values.append(
            ((L - 2 * candidate_radius) * center + 2 * np.dot(weights, side))
            / L
        )
    return np.asarray(values)


def clean_feedback_high_density(L: float, q: int, candidate_radius: float) -> float:
    """High-density square law including side and corner correction."""
    return float(
        q
        * (q + 1)
        / 2
        * (L / (2 * candidate_radius) + 2 * np.log(2) - 1)
    )
