import numpy as np
from scipy.sparse import csr_matrix
from janus_connectivity.feedback_order import classify_graph, make_filtered_graph
from janus_connectivity.scc import largest_scc_fraction


def test_g1_g2_g_largest_scc_monotonicity():
    edges = [(0, 1), (1, 0), (2, 3), (3, 4), (4, 5), (5, 2), (1, 2)]
    src, dst = np.asarray(edges).T
    G = csr_matrix((np.ones(len(edges), np.uint8), (src, dst)), shape=(6, 6))
    y = np.array([1.0, 0.0, -1.0, 3.0, 1.0, 2.0])
    result = classify_graph(G, y)
    G1 = make_filtered_graph(result, False)
    G2 = make_filtered_graph(result, True)
    assert largest_scc_fraction(G1) <= largest_scc_fraction(G2)
    assert largest_scc_fraction(G2) <= largest_scc_fraction(G)
