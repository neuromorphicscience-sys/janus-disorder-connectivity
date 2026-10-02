"""Integrity checks on the complete archived SI cohorts; no graph generation."""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]/'data/supplemental_final'
read=lambda name:pd.read_csv(ROOT/(name+'.csv'),engine='python')

def test_ell_common_fidelity_and_missing_rejections():
    d=read('source_figS5')
    assert len(d)==192 and d.graph_id.nunique()==48
    metrics=[x for x in d if x.endswith('_KS') or x.endswith('_Wasserstein')]
    passed=(d[metrics]<=.03).all(axis=1)&(d.I0_relative_change<=.03)&(d.changed_edge_fraction>=.25)&d.outdegree_exact&d.orientation_multiset_exact
    assert np.array_equal(passed,d.accepted)
    assert passed.sum()==186
    assert d[~passed].S_reassigned.isna().all()
    assert d[passed].S_reassigned.lt(d[passed].S_original).all()

def test_complete_internal_cohort_and_product():
    d=read('source_figS8')
    assert len(d)==1528 and d.edge_hash.nunique()==1528 and d.W.nunique()==24
    for k in [1,2]:
        np.testing.assert_allclose(d[f'S_internal_rec_core_{k}'],d[f'f_internal_rec_{k}']*d[f'eta_internal_{k}'],rtol=0,atol=1e-14)

def test_complete_deletion_population_and_widths():
    d=read('source_figS9')
    assert len(d)==576 and d.edge_hash.nunique()==192
    assert d.groupby('edge_hash').width_over_R.apply(lambda x:set(x)=={2,4,8}).all()
    assert d[d.width_over_R.eq(4)].groupby('W').size().tolist()==[96,48,24,24]
    np.testing.assert_allclose(d.delta_full_G2,d.S_full_core-d.S2_core,rtol=0,atol=1e-14)
    assert d.delta_full_G2.gt(0).all()

def test_complete_size_samples_and_covariant_direction_coordinates():
    d=read('source_figS3')
    assert len(d)==4634 and d.edge_sha256.nunique()==4634
    assert d.groupby(['L','W']).size().max()==150
    d=read('source_figS7_direction')
    assert len(d)==72 and d.graph_pair_id.nunique()==24
    assert set(d.axis_policy)=={'COVARIANT_AXIS'}
    for k in [1,2]:
        np.testing.assert_allclose(d[f'DLOC_G{k}'],d[f'G{k}_downstream_2R_fraction']-d[f'G{k}_upstream_2R_fraction'],rtol=0,atol=1e-14)
