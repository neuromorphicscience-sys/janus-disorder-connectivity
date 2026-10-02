"""Prespecified fidelity gates and descriptive crossing semantics."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from replay_final_hardening import accepts_fidelity
from summarize_supplemental import bracket_crossing


def record():
    r=dict(outdegree_exact=True,orientation_multiset_exact=True,changed_edge_fraction=.25,I0_relative_change=.03)
    for x in ['length','relative_angle','dy']:
        r[x+'_KS']=.03;r[x+'_Wasserstein']=.03
    return r


def test_common_gate_includes_limits_and_requires_changes():
    r=record();assert accepts_fidelity(r)
    r['changed_edge_fraction']=.249999;assert not accepts_fidelity(r)
    r=record();r['changed_edge_fraction']=.99;assert accepts_fidelity(r)


def test_each_fidelity_constraint_and_multiset_are_required():
    for key in ['length_KS','relative_angle_KS','dy_KS','length_Wasserstein','relative_angle_Wasserstein','dy_Wasserstein','I0_relative_change']:
        r=record();r[key]=.030001;assert not accepts_fidelity(r)
    for key in ['outdegree_exact','orientation_multiset_exact']:
        r=record();r[key]=False;assert not accepts_fidelity(r)


def test_crossing_is_interpolated_and_never_extrapolated():
    value,lo,hi,status=bracket_crossing([.7,.8,.9],[.1,.4,.9],.5)
    assert abs(value-.82)<1e-14 and (lo,hi,status)==(.8,.9,'BRACKETED')
    assert np.isnan(bracket_crossing([.7,.8],[.1,.4],.5)[0])


def test_multiple_upward_brackets_are_ambiguous():
    result=bracket_crossing([1,2,3,4],[.1,.8,.2,.9],.5)
    assert np.isnan(result[0]) and result[-1]=='NA_MULTIPLE_CROSSINGS'


def test_frozen_counts_and_rejected_outcome_policy():
    root=Path(__file__).resolve().parents[1]/'data/supplemental'
    e=pd.read_csv(root/'02_ell_sensitivity/GRAPHLEVEL.csv')
    assert len(e)==192 and e.groupby('ell_over_R').size().eq(48).all()
    assert e.loc[~e.accepted,'S_reassigned'].isna().all()
    assert e.loc[~e.accepted,'log10_S_ratio'].isna().all()
    for row in e.to_dict('records'):assert row['accepted']==accepts_fidelity(row)
    assert e[e.ell_over_R.eq(4)].accepted.all()
    d=pd.read_csv(root/'03_boundary_deletion_ensemble/GRAPHLEVEL.csv')
    assert len(d)==576 and not d.duplicated(['edge_hash','width_over_R']).any()
    for width,g in d.groupby('width_over_R'):
        assert g.groupby('W').size().to_dict()=={.805:96,.815:48,.86:24,.9:24}
    assert d.S1_core.le(d.S2_core).all() and d.S2_core.le(d.S_full_core).all()
