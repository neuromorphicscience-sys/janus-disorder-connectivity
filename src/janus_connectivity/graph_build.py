"""Exact top-q projection graph construction with cell-list neighbor search."""
from __future__ import annotations
import numpy as np
from numba import njit
from scipy.sparse import csr_matrix


def prepare_cells(points: np.ndarray, R: float = 1.0):
    """Build cell-list arrays for candidate radius a=2R."""
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("finite N x 2 points required")
    if not np.isfinite(R) or R <= 0:
        raise ValueError("positive R required")
    cell = np.floor(points / (2 * R)).astype(np.int64)
    origin = cell.min(axis=0) if len(cell) else np.zeros(2, dtype=np.int64)
    cell -= origin
    shape = cell.max(axis=0) + 1 if len(cell) else np.ones(2, dtype=np.int64)
    ids = cell[:, 0] + shape[0] * cell[:, 1]
    order = np.argsort(ids, kind="stable").astype(np.int32)
    count = np.bincount(ids, minlength=int(shape.prod()))
    ptr = np.r_[0, np.cumsum(count)].astype(np.int64)
    return cell, order, ptr, int(shape[0]), int(shape[1])


@njit(cache=True, fastmath=False)
def _select(points, cosine, sine, q, radius2, cell, order, ptr, nx, ny):
    n = len(points)
    chosen = np.full((n, q), -1, np.int32)
    counts = np.zeros(n, np.int32)
    keys = np.empty(q, np.float64)
    targets = np.empty(q, np.int32)
    visits = 0
    for i in range(n):
        k = 0
        cx, cy = cell[i, 0], cell[i, 1]
        px, py = points[i, 0], points[i, 1]
        co, si = cosine[i], sine[i]
        for yy in range(max(0, cy - 1), min(ny, cy + 2)):
            for xx in range(max(0, cx - 1), min(nx, cx + 2)):
                c = xx + nx * yy
                for pos in range(ptr[c], ptr[c + 1]):
                    j = order[pos]
                    if j == i:
                        continue
                    visits += 1
                    dx = points[j, 0] - px
                    dy = points[j, 1] - py
                    if dx * dx + dy * dy >= radius2:
                        continue
                    counts[i] += 1
                    if q == 0:
                        continue
                    # Smallest key is largest projection on e_i=(sin theta,-cos theta).
                    key = points[j, 1] * co - points[j, 0] * si
                    rank = k
                    while rank > 0 and (
                        key < keys[rank - 1]
                        or (key == keys[rank - 1] and j < targets[rank - 1])
                    ):
                        rank -= 1
                    if rank >= q:
                        continue
                    for z in range(min(k, q - 1), rank, -1):
                        keys[z] = keys[z - 1]
                        targets[z] = targets[z - 1]
                    keys[rank] = key
                    targets[rank] = j
                    if k < q:
                        k += 1
        # Canonical target ordering produces a stable CSR representation.
        for j in range(k):
            chosen[i, j] = targets[j]
        for j in range(1, k):
            value = chosen[i, j]
            z = j
            while z > 0 and chosen[i, z - 1] > value:
                chosen[i, z] = chosen[i, z - 1]
                z -= 1
            chosen[i, z] = value
    return chosen, counts, visits


def build_directed_graph(
    points: np.ndarray,
    theta: np.ndarray,
    q: int = 7,
    R: float = 1.0,
    cells=None,
):
    """Construct the directed graph under candidate radius 2R and top-q ranking."""
    points = np.asarray(points, dtype=np.float64)
    theta = np.asarray(theta, dtype=np.float64)
    if theta.shape != (len(points),) or not np.isfinite(theta).all():
        raise ValueError("theta must be a finite vector with one value per point")
    if q < 0 or q != int(q):
        raise ValueError("q must be a nonnegative integer")
    cells = prepare_cells(points, R) if cells is None else cells
    selected, counts, visits = _select(
        points, np.cos(theta), np.sin(theta), int(q), (2 * R) ** 2, *cells
    )
    degree = np.minimum(counts, q)
    indptr = np.r_[0, np.cumsum(degree)].astype(np.int64)
    indices = selected[selected >= 0]
    graph = csr_matrix(
        (np.ones(len(indices), np.uint8), indices, indptr),
        shape=(len(points), len(points)),
    )
    return graph, counts, int(visits)


def brute_force_graph(points: np.ndarray, theta: np.ndarray, q: int = 7, R: float = 1.0):
    """Small-N reference implementation used in tests."""
    points = np.asarray(points, dtype=float)
    theta = np.asarray(theta, dtype=float)
    edges = []
    for i, point in enumerate(points):
        direction = np.array([np.sin(theta[i]), -np.cos(theta[i])])
        candidates = [
            j
            for j in range(len(points))
            if j != i and np.linalg.norm(points[j] - point) < 2 * R
        ]
        candidates.sort(
            key=lambda j: (-float(np.dot(points[j] - point, direction)), j)
        )
        edges.extend((i, j) for j in candidates[:q])
    if edges:
        src, dst = np.asarray(edges, dtype=np.int64).T
    else:
        src = dst = np.empty(0, dtype=np.int64)
    return csr_matrix(
        (np.ones(len(src), np.uint8), (src, dst)),
        shape=(len(points), len(points)),
    )
