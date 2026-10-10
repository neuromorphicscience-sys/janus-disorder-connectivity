"""Regression cases for the geometric bounds, independent of archived labels."""
import runpy
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

VALIDATE = runpy.run_path(str(Path(__file__).resolve().parents[1] /
                             "scripts/validate_full_graph_classifier.py"))


def test_all_four_vertex_topologies():
    for _, graph, y in VALIDATE["four_vertex_graphs"]():
        VALIDATE["check_graph"](graph, y)


def test_reverse_ceiling_equality_and_shared_endpoint():
    # i=0->1, j=1->2; v_j=y(v_i)+max_up_dy exactly. D return 2->0.
    graph = csr_matrix((np.ones(3, np.uint8), ([0, 1, 2], [1, 2, 0])), shape=(3, 3))
    result, labels, _ = VALIDATE["check_graph"](graph, np.array([0., 1., 2.]))
    assert result["order2"] == 2
    assert labels.tolist() == [2, 2, 0]


def test_forward_floor_equality_and_order_one_partner():
    # Candidate source 2 is exactly the forward floor. Partner 2->3 has
    # its own D-only return; it must remain eligible to close 0->4.
    edges = [(0, 4), (4, 2), (2, 3), (3, 2), (3, 0), (4, 1)]
    src, dst = zip(*edges)
    graph = csr_matrix((np.ones(len(edges), np.uint8), (src, dst)), shape=(5, 5))
    record, _, _ = VALIDATE["check_graph"](graph, np.arange(5, dtype=float))
    assert record["order1"] == 1 and record["order2"] == 1


def test_negative_coordinate_ceiling_cancellation():
    # h=fl(0.2-(-1.0)) rounds down. nextafter(fl(-1.0+h)) is still
    # below 0.2, so rounding only the final sum would miss edge 0->1.
    graph = csr_matrix((np.ones(3, np.uint8), ([0, 1, 2], [1, 2, 0])), shape=(3, 3))
    record, labels, _ = VALIDATE["check_graph"](graph, np.array([-1.5, -1., .2]))
    assert record["order2"] == 2 and labels.tolist() == [2, 2, 0]


def test_progress_coordinate_translation_and_scaling():
    graph = csr_matrix((np.ones(3, np.uint8), ([0, 1, 2], [1, 2, 0])), shape=(3, 3))
    for scale in (2.**-20, 1., 2.**20):
        for shift in (-3., -1., 0., 3.):
            _, labels, _ = VALIDATE["check_graph"](
                graph, scale * (np.array([-.5, 0., 1.2]) + shift))
            assert labels.tolist() == [2, 2, 0]


def test_complete_model_suite_is_fixed():
    cases = list(VALIDATE["model_cases"]())
    assert len(cases) == 26
    assert {case["N"] for case in cases} == {96, 384, 864, 1536}
