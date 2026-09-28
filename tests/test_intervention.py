import numpy as np
from janus_connectivity.intervention import reassign_fixed_multiset


def test_fixed_multiset_reassignment():
    rng = np.random.default_rng(4)
    points = rng.random((100, 2)) * 16
    theta = rng.normal(size=100)
    reassigned, field = reassign_fixed_multiset(points, theta, 16, 4, seed=18, grid_size=64)
    assert np.array_equal(np.sort(reassigned), np.sort(theta))
    assert np.unique(field).size > 90
