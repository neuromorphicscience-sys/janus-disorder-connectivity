import numpy as np
from janus_connectivity.graph_build import build_directed_graph, brute_force_graph


def test_cell_graph_matches_brute_force():
    rng = np.random.default_rng(1234)
    points = rng.random((80, 2)) * 8
    theta = rng.normal(0, 0.8, 80)
    fast, counts, _ = build_directed_graph(points, theta, q=4, R=1.0)
    slow = brute_force_graph(points, theta, q=4, R=1.0)
    assert (fast != slow).nnz == 0
    assert np.all(np.diff(fast.indptr) == np.minimum(counts, 4))


def test_outdegree_budget():
    rng = np.random.default_rng(8)
    points = rng.random((60, 2)) * 5
    theta = rng.normal(size=60)
    graph, _, _ = build_directed_graph(points, theta, q=7, R=1)
    assert np.diff(graph.indptr).max() <= 7
    assert graph.diagonal().sum() == 0
