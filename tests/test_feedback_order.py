import numpy as np
from scipy.sparse import csr_matrix
from janus_connectivity.feedback_order import (
    ORDER1, ORDER2, GT2_OR_INF, classify_graph, make_filtered_graph
)


def graph(n, edges):
    src, dst = np.asarray(edges, dtype=int).T
    return csr_matrix((np.ones(len(edges), np.uint8), (src, dst)), shape=(n, n))


def edge_labels(result):
    return {(int(u), int(v)): int(label) for u, v, label in zip(result["u"], result["v"], result["labels"])}


def test_order1_edge():
    G = graph(2, [(0, 1), (1, 0)])
    result = classify_graph(G, np.array([1.0, 0.0]))
    assert list(result["labels"]) == [ORDER1]


def test_reciprocal_order2_pair():
    G = graph(4, [(0, 3), (3, 1), (1, 2), (2, 0)])
    result = classify_graph(G, np.array([0.0, 1.0, 2.0, 3.0]))
    assert np.all(result["labels"] == ORDER2)
    assert result["witness"][result["witness"][0]] == 0


def test_order2_partner_may_be_order1():
    # i=(0,3) closes only through j=(1,2); j has its own D-only return 2->1.
    G = graph(4, [(0, 3), (1, 2), (3, 1), (2, 1), (2, 0)])
    result = classify_graph(G, np.array([0.0, 1.0, 2.0, 3.0]))
    labels = edge_labels(result)
    assert labels[(0, 3)] == ORDER2
    assert labels[(1, 2)] == ORDER1
    i = np.flatnonzero((result["u"] == 0) & (result["v"] == 3))[0]
    partner = int(result["witness"][i])
    assert result["labels"][partner] == ORDER1


def test_true_order3_is_not_promoted():
    G = graph(6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)])
    result = classify_graph(G, np.array([0.0, 5.0, 1.0, 4.0, 2.0, 3.0]))
    assert np.all(result["labels"] == GT2_OR_INF)


def test_edge_nesting():
    G = graph(6, [(0, 1), (1, 0), (2, 3), (3, 4), (4, 5), (5, 2)])
    result = classify_graph(G, np.array([1.0, 0.0, 0.0, 3.0, 1.0, 2.0]))
    G1 = make_filtered_graph(result, include_order2=False)
    G2 = make_filtered_graph(result, include_order2=True)
    e1 = set(zip(*G1.nonzero()))
    e2 = set(zip(*G2.nonzero()))
    eg = set(zip(*G.nonzero()))
    assert e1 <= e2 <= eg
