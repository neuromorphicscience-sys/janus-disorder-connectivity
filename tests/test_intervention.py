import numpy as np
from janus_connectivity.intervention import reassign_fixed_multiset


def test_fixed_multiset_reassignment():
    rng = np.random.default_rng(4)
    points = rng.random((100, 2)) * 16
    theta = rng.normal(size=100)
    reassigned, field = reassign_fixed_multiset(points, theta, 16, 4, seed=18, grid_size=64)
    assert np.array_equal(np.sort(reassigned), np.sort(theta))
    assert np.unique(field).size > 90


def test_fidelity_uses_wrapped_source_angle_and_pi_normalizer():
    from janus_connectivity.intervention import marginal_fidelity
    a={"length":np.array([1.,1.]),"source_relative_angle":np.array([0.,0.]),"vertical_displacement":np.array([-1.,-1.])}
    b={**a,"source_relative_angle":np.array([.1,.1])}
    result=marginal_fidelity(a,b)
    assert abs(result["source_relative_angle_wasserstein_normalized"]-.1/np.pi)<1e-14
    assert result["source_relative_angle_ks"]==1
