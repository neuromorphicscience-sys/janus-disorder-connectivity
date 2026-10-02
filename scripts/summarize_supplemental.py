"""Recompute supplemental descriptive summaries from the supplied graph-level data."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
PROJECT=Path(__file__).resolve().parents[1]
INPUT=PROJECT/'data/supplemental'
OUT=PROJECT/'outputs/supplemental_summaries'
KEYS=['distribution','correlation_family','control_name','direction','q','sigma','L','ell_x','ell_y','control']

def summarize_values(data,keys):
    rows=[]
    for key,g in data.groupby(keys,dropna=False):
        key=key if isinstance(key,tuple) else (key,);x=g.S.to_numpy()
        row=dict(zip(keys,key));row.update(n=len(x),S_median=np.median(x),S_Q1=np.quantile(x,.25),S_Q3=np.quantile(x,.75),S_min=x.min(),S_max=x.max(),S_variance=np.var(x,ddof=1) if len(x)>1 else np.nan);rows.append(row)
    return pd.DataFrame(rows)

def bracket_crossing(x,y,level):
    """Return a unique adjacent upward bracket, else NaN. No extrapolation."""
    x=np.asarray(x);y=np.asarray(y);brackets=[]
    for i in range(len(x)-1):
        if y[i]<=level<=y[i+1] and y[i+1]>y[i]:
            brackets.append((x[i]+(level-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]),x[i],x[i+1]))
    if len(brackets)==1:return (*brackets[0],'BRACKETED')
    return (np.nan,np.nan,np.nan,'NA_NOT_BRACKETED' if not brackets else 'NA_MULTIPLE_CROSSINGS')

def compare(actual,expected,name):
    assert list(actual.columns)==list(expected.columns),name+' columns'
    assert len(actual)==len(expected),name+' row count'
    for col in actual:
        if pd.api.types.is_numeric_dtype(actual[col]) and pd.api.types.is_numeric_dtype(expected[col]):
            np.testing.assert_allclose(actual[col],expected[col],atol=2e-12,rtol=2e-10,equal_nan=True,err_msg=name+': '+col)
        else:
            assert actual[col].fillna('').astype(str).tolist()==expected[col].fillna('').astype(str).tolist(),name+': '+col

def save_check(d,folder,name):
    dest=OUT/folder;dest.mkdir(parents=True,exist_ok=True);d.to_csv(dest/name,index=False)
    compare(d,pd.read_csv(INPUT/folder/name),name)

def main():
    r=pd.read_csv(INPUT/'04_parameter_robustness/ROBUSTNESS_GRAPHLEVEL.csv')
    save_check(summarize_values(r,KEYS),'04_parameter_robustness','ROBUSTNESS_SUMMARY.csv')
    f=pd.read_csv(INPUT/'05_finite_size_diagnostics/FINITE_SIZE_GRAPHLEVEL.csv')
    s=summarize_values(f,['L','W']);rng=np.random.default_rng(2026100201)
    for i,row in s.iterrows():
        x=f[f.L.eq(row.L)&f.W.eq(row.W)].S.to_numpy();med=np.median(x[rng.integers(0,len(x),(2000,len(x)))],axis=1)
        s.loc[i,'median_bootstrap_low']=np.quantile(med,.025);s.loc[i,'median_bootstrap_high']=np.quantile(med,.975)
    save_check(s,'05_finite_size_diagnostics','FINITE_SIZE_SUMMARY.csv')
    cross=[]
    for L,g in s.groupby('L'):
        g=g.sort_values('W');row={'L':int(L),'n_graphs':int(g.n.sum()),'n_W':len(g)}
        for lev in [.1,.5,.9]:
            value,lo,hi,status=bracket_crossing(g.W,g.S_median,lev);key=f'W_{lev:g}'
            row[key]=value;row[key+'_status']=status;row[key+'_bracket_low']=lo;row[key+'_bracket_high']=hi
        row['Delta_W_10_90']=row['W_0.9']-row['W_0.1'];near=g[g.W.between(.75,.85)];peak=near.loc[near.S_variance.idxmax()]
        row['sampled_variance_max_W']=float(peak.W);row['sampled_variance_max']=float(peak.S_variance)
        row['variance_peak_interpretation']='Sampled maximum only; unequal grids and n prevent a resolved peak estimate'
        row['median_downward_steps']=int(np.sum(np.diff(g.S_median)<0));cross.append(row)
    save_check(pd.DataFrame(cross),'05_finite_size_diagnostics','CROSSOVER_DIAGNOSTICS.csv')
    e=pd.read_csv(INPUT/'02_ell_sensitivity/GRAPHLEVEL.csv');rows=[]
    for k,g in e.groupby('ell_over_R'):
        a=g[g.accepted];x=a.log10_S_ratio.to_numpy()
        rows.append(dict(ell_over_R=int(k),attempted_n=len(g),accepted_n=len(a),rejected_n=len(g)-len(a),log10_ratio_median=np.median(x),log10_ratio_Q1=np.quantile(x,.25),log10_ratio_Q3=np.quantile(x,.75),largest_accepted_S_ratio=10**x.max(),KS_max_accepted=a[[c for c in a if c.endswith('_KS')]].to_numpy().max(),Wasserstein_max_accepted=a[[c for c in a if c.endswith('_Wasserstein')]].to_numpy().max(),I0_relative_change_max_accepted=a.I0_relative_change.max(),changed_fraction_min_accepted=a.changed_edge_fraction.min(),status='SUPPORTED' if (x<0).all() else 'NOT SUPPORTED'))
    save_check(pd.DataFrame(rows),'02_ell_sensitivity','SUMMARY.csv')
    d=pd.read_csv(INPUT/'03_boundary_deletion_ensemble/GRAPHLEVEL.csv');rows=[]
    for (w,width),g in d.groupby(['W','width_over_R']):
        row=dict(W=w,width_over_R=width,n=len(g),n_positive_difference=int(g.delta_full_G2.gt(0).sum()),n_equal=int(g.delta_full_G2.eq(0).sum()))
        for name in ['S_full_core','S1_core','S2_core','delta_full_G2','ratio_full_G2','retained_fraction']:
            x=g[name].to_numpy()
            for key,val in [('median',np.median(x)),('Q1',np.quantile(x,.25)),('Q3',np.quantile(x,.75)),('min',x.min()),('max',x.max())]:row[name+'_'+key]=float(val)
        rows.append(row)
    save_check(pd.DataFrame(rows),'03_boundary_deletion_ensemble','SUMMARY.csv')
    print('PASS: parameter, bootstrap, crossover, intervention and deletion summaries reproduce from graph-level tables')
if __name__=='__main__':main()
