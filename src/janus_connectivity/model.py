"""Fixed-N point and quenched-orientation model."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .graph_build import build_directed_graph


@dataclass(frozen=True)
class ModelParameters:
    L: float = 128.0
    sigma: float = 24.0
    q: int = 7
    R: float = 1.0
    W: float = 0.805

    @property
    def N(self) -> int:
        return round(self.sigma * self.L * self.L)


def wrap_angles(theta: np.ndarray) -> np.ndarray:
    wrapped = (np.asarray(theta, dtype=float) + np.pi) % (2 * np.pi) - np.pi
    return np.where(wrapped <= -np.pi, np.pi, wrapped)


def generate_points(N: int, L: float, seed: int) -> np.ndarray:
    """Generate N independent uniform points in the open square [0,L]^2."""
    if N < 0 or L <= 0:
        raise ValueError("N must be nonnegative and L must be positive")
    return np.random.default_rng(seed).random((N, 2), dtype=np.float64) * float(L)


def generate_orientations(N: int, W: float, seed: int) -> np.ndarray:
    """Generate theta_i = W z_i modulo 2 pi, z_i standard normal."""
    if N < 0 or W < 0:
        raise ValueError("N and W must be nonnegative")
    return wrap_angles(float(W) * np.random.default_rng(seed).standard_normal(N))


def generate_realization(
    parameters: ModelParameters,
    geometry_seed: int,
    orientation_seed: int,
):
    """Generate positions, orientations, and the projection-ranked graph."""
    points = generate_points(parameters.N, parameters.L, geometry_seed)
    theta = generate_orientations(parameters.N, parameters.W, orientation_seed)
    graph, candidate_counts, visits = build_directed_graph(
        points, theta, q=parameters.q, R=parameters.R
    )
    return points, theta, graph, candidate_counts, visits
