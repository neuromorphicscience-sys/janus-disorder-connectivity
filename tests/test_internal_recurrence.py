import numpy as np
from scipy.sparse import csr_matrix
from janus_connectivity.recurrence import internal_recurrence


def test_internal_recurrence_identity():
    edges = [(0, 1), (1, 0), (2, 3), (3, 2), (4, 5)]
    src, dst = np.asarray(edges).T
    graph = csr_matrix((np.ones(len(edges)), (src, dst)), shape=(6, 6))
    result = internal_recurrence(graph, np.ones(6, dtype=bool))
    assert result["recurrent_mass"] == 4
    assert result["largest_recurrent_component"] == 2
    assert result["f_internal_rec"] == 4 / 6
    assert result["eta_internal"] == 1 / 2
    assert abs(result["identity_residual"]) < 1e-15
