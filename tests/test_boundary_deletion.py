import numpy as np
from scipy.sparse import csr_matrix
from janus_connectivity.recurrence import induce_vertices
from janus_connectivity.scc import scc_decomposition


def largest_size(graph):
    _, sizes = scc_decomposition(graph)
    return int(sizes.max()) if len(sizes) else 0


def test_boundary_deletion_cannot_increase_largest_scc_mass():
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 4)]
    src, dst = np.asarray(edges).T
    graph = csr_matrix((np.ones(len(edges)), (src, dst)), shape=(6, 6))
    subgraph, _ = induce_vertices(graph, np.array([False, True, True, True, True, True]))
    assert largest_size(subgraph) <= largest_size(graph)


def test_filtered_core_remains_subgraph_of_full_core():
    full = csr_matrix(np.array([[0,1,0],[0,0,1],[1,0,0]], dtype=np.uint8))
    filtered = csr_matrix(np.array([[0,1,0],[0,0,1],[0,0,0]], dtype=np.uint8))
    mask = np.array([True, True, True])
    full_core, _ = induce_vertices(full, mask)
    filtered_core, _ = induce_vertices(filtered, mask)
    filtered_edges = set(zip(*filtered_core.nonzero()))
    full_edges = set(zip(*full_core.nonzero()))
    assert filtered_edges <= full_edges
