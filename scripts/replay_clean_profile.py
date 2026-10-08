"""Replay exactly the existing 24 clean graphs; source-depth mean comparison."""
from pathlib import Path
import sys,os,json,hashlib
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[name]='1'
R=Path(__file__).resolve().parents[1];P=R;F=R/'outputs/final_bridge'
sys.path.insert(0,str(R/'src'))
import numpy as np,pandas as pd
from janus_connectivity.model import generate_points,generate_orientations
from janus_connectivity.graph_build import build_directed_graph
from concurrent.futures import ProcessPoolExecutor,as_completed
def job(row):
    dest=F/f"04_clean_profile/records/{row['graph_id']}.csv"
    if dest.exists():return str(dest)
    n=round(row['sigma']*row['L']**2);p=generate_points(n,row['L'],int(row['geometry_seed']));th=generate_orientations(n,0,int(row['disorder_seed']));g,_,_=build_directed_graph(p,th,row['q'])
    assert hashlib.sha256(g.indptr.tobytes()+g.indices.tobytes()).hexdigest()==row['expected_hash']
    source=np.repeat(np.arange(n),np.diff(g.indptr));up=source[p[g.indices,1]>p[source,1]]
    assert len(up)==int(row['F_count'])
    bins=np.linspace(0,.22,12);nc,_=np.histogram(p[:,1],bins);ec,_=np.histogram(p[up,1],bins)
    pd.DataFrame(dict(graph_id=row['graph_id'],edge_hash=row['expected_hash'],d_min_R=bins[:-1],d_max_R=bins[1:],d_center_R=(bins[1:]+bins[:-1])/2,source_node_n=nc,upward_edge_n=ec,mean_Z_up=np.divide(ec,nc,out=np.zeros_like(ec,dtype=float),where=nc>0),lateral_average='all x in [0,L]')).to_csv(dest,index=False)
    return str(dest)
def main():
    (F/'04_clean_profile/records').mkdir(parents=True,exist_ok=True)
    pop=json.loads((R/'data/final_bridge/frozen_replay_inputs.json').read_text())
    counts=pd.read_csv(R/'data/figure4/clean_feedback_counts.csv').set_index('graph_id')
    parents={r['graph_id']:r for r in pop['populations']['primary']}
    rows=[]
    for gid,c in counts.iterrows():
        row=dict(parents[gid]);assert row['expected_hash']==c.edge_hash and row['control']==0
        row['F_count']=int(c.F_count);rows.append(row)
    assert len(rows)==24
    pd.DataFrame(rows).to_csv(F/'04_clean_profile/CLEAN_REPLAY_MANIFEST.csv',index=False)
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(job,row) for row in rows]):print('DONE',future.result(),flush=True)
    d=pd.concat([pd.read_csv(f) for f in sorted((F/'04_clean_profile/records').glob('*.csv'))]);d.to_csv(F/'04_clean_profile/CLEAN_PROFILE_SOURCE.csv',index=False)
    theory=pd.read_csv(R/'data/figure4/boundary_depth_prediction.csv');out=[]
    for depth,x in d.groupby('d_center_R'):
        lo=x.d_min_R.iloc[0];hi=x.d_max_R.iloc[0];grid=np.linspace(lo,hi,201)
        pred=np.trapezoid(np.interp(grid,theory.d_over_R,theory.mean_Z_up_fixedN_x_averaged),grid)/(hi-lo)
        out.append(dict(d_center_R=depth,d_min_R=lo,d_max_R=hi,n=len(x),mean_Z_up=x.mean_Z_up.mean(),SEM_Z_up=x.mean_Z_up.std(ddof=1)/np.sqrt(len(x)),Q1=x.mean_Z_up.quantile(.25),Q3=x.mean_Z_up.quantile(.75),theory_bin_average=pred))
    pd.DataFrame(out).to_csv(F/'04_clean_profile/CLEAN_PROFILE_SUMMARY.csv',index=False)
    (F/'04_clean_profile/CLEAN_PROFILE_REPORT.md').write_text('# Clean-source profile\n\nAll 24 archived W=0 graphs were replayed with exact edge hashes and feedback counts. Each bin averages outgoing upward counts over sources at depth d=y, including the full lateral interval 0<=x<=L, matching the fixed-N theory. Ensemble mean and SEM (24 independent clean geometries) are used because the prediction is an ensemble mean. Theory is averaged within each 0.02R bin, preventing a midpoint-only comparison. No fitted scale or new realization is used. Raw graph/bin node and edge counts are supplied.\n')
if __name__=='__main__':main()
