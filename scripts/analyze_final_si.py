"""Recompute descriptive SI statistics without generating any base networks."""
from pathlib import Path
import argparse,json
import pandas as pd
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
read=lambda n:pd.read_csv(a.data/(n+'.csv'),engine='python')
def summ(d,groups,columns):
 out=[]
 for key,g in d.groupby(groups,dropna=False):
  if not isinstance(key,tuple):key=(key,)
  r=dict(zip(groups,key));r['n']=len(g)
  for col in columns:
   x=g[col]
   for suffix,fun in [('median',lambda x:x.median()),('Q1',lambda x:x.quantile(.25)),('Q3',lambda x:x.quantile(.75)),('min',lambda x:x.min()),('max',lambda x:x.max())]:r[col+'_'+suffix]=fun(x)
  out.append(r)
 return pd.DataFrame(out)
ell=read('source_figS5');dist=[x for x in ell if x.endswith('_KS') or x.endswith('_Wasserstein')]
passed=(ell[dist]<=.03).all(axis=1)&(ell.I0_relative_change<=.03)&(ell.changed_edge_fraction>=.25)&ell.outdegree_exact&ell.orientation_multiset_exact
assert (passed==ell.accepted).all();assert ell[~passed].S_reassigned.isna().all()
ell['S_ratio']=ell.S_reassigned/ell.S_original
ell_summary=summ(ell[passed],['ell_over_R'],['log10_S_ratio','S_ratio','I0_relative_change','changed_edge_fraction'])
ell_summary['attempted_n']=ell_summary.ell_over_R.map(ell.groupby('ell_over_R').size());ell_summary['rejected_n']=ell_summary.attempted_n-ell_summary.n
ell_summary.to_csv(a.output/'ell_summary.csv',index=False)
fidelity_rows=[]
metrics=[(f'{key}_{kind}',('none' if kind=='KS' else norm),.03,'<=') for key,norm in [('length','2R'),('relative_angle','pi'),('dy','2R')] for kind in ['KS','Wasserstein']]+[('I0_relative_change','max(I0_original,0.01)',.03,'<='),('changed_edge_fraction','original edge count',.25,'>=')]
for scale,g in ell.groupby('ell_over_R'):
 for key,norm,limit,relation in metrics:
  x=g[key];fidelity_rows.append(dict(ell_over_R=scale,metric=key,normalization=norm,threshold=limit,relation=relation,median=x.median(),Q1=x.quantile(.25),Q3=x.quantile(.75),minimum=x.min(),maximum=x.max(),pass_count=int((x<=limit).sum() if relation=='<=' else (x>=limit).sum()),total_count=len(x)))
pd.DataFrame(fidelity_rows).to_csv(a.output/'source_tableS2_recomputed.csv',index=False)
fs=read('source_figS3');fss=summ(fs,['L','W'],['S']);fss['variance']=fs.groupby(['L','W']).S.var().to_numpy();fss.to_csv(a.output/'finite_size_summary.csv',index=False)
delete=read('source_figS9');assert len(delete)==576;assert delete.edge_hash.nunique()==192
assert np.allclose(delete.delta_full_G2,delete.S_full_core-delete.S2_core,rtol=0,atol=1e-14)
summ(delete,['W','width_over_R'],['S_full_core','S1_core','S2_core','delta_full_G2','retained_fraction']).to_csv(a.output/'deletion_summary.csv',index=False)
internal=read('source_figS8');assert len(internal)==1528 and internal.W.nunique()==24
for k in [1,2]:
 assert np.allclose(internal[f'S_internal_rec_core_{k}'],internal[f'f_internal_rec_{k}']*internal[f'eta_internal_{k}'],rtol=0,atol=1e-14)
 internal[f'retention_{k}']=internal[f'internal_recurrent_vertices_{k}']/internal[f'parent_recurrent_vertices_{k}']
summ(internal,['W'],[f'{prefix}{k}' for prefix in ['f_internal_rec_','eta_internal_','S_internal_rec_core_','internal_nontrivial_scc_count_','retention_'] for k in [1,2]]).to_csv(a.output/'internal_summary.csv',index=False)
retention_rows=[]
for w,g in internal.groupby('W'):
 for k in [1,2]:
  r=dict(W=w,n=len(g),k=k)
  for name,x in [('parent_fraction',g[f'f_parent_rec_{k}']),('internal_fraction',g[f'f_internal_rec_{k}']),('retention',g[f'retention_{k}'])]:
   r.update({name+'_median':x.median(),name+'_Q1':x.quantile(.25),name+'_Q3':x.quantile(.75)})
  retention_rows.append(r)
pd.DataFrame(retention_rows).to_csv(a.output/'source_tableS3_recomputed.csv',index=False)
direction=read('source_figS7_direction');assert len(direction)==72
assert set(direction.axis_policy)=={'COVARIANT_AXIS'}
for n in [1,2]:assert np.allclose(direction[f'DLOC_G{n}'],direction[f'G{n}_downstream_2R_fraction']-direction[f'G{n}_upstream_2R_fraction'])
validation=read('source_figS6');assert len(validation)==45 and validation.classification_exact.all()
report={'new_base_network_realizations':0,'new_surrogate_fields':0,'G3_production':0,'ell_attempts':len(ell),'ell_accepted':int(passed.sum()),'ell_rejected':int((~passed).sum()),'all_accepted_suppress_connectivity':bool((ell[passed].S_ratio<1).all()),'common_fidelity_gates_pass':True,'internal_realizations':len(internal),'internal_W_conditions':internal.W.nunique(),'deletion_rows':len(delete),'deletion_unique_graphs':delete.edge_hash.nunique(),'all_deletion_differences_positive':bool((delete.delta_full_G2>0).all()),'direction_pairs':24,'archived_independent_classifier_cases':45,'contradictions':[]}
(a.output/'analysis_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
