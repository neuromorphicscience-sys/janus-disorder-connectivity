"""Prespecified replay and derived analyses; no new base seeds or order-three work."""
from __future__ import annotations
import argparse, hashlib, json, os, sys, time
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
sys.dont_write_bytecode=True
PROJECT=Path(__file__).resolve().parents[1]
INPUT=PROJECT/'data/supplemental'
ROOT=PROJECT/'outputs/final_hardening'
os.environ['NUMBA_CACHE_DIR']=str(PROJECT/'.numba_cache')
sys.path.insert(0,str(PROJECT/'src'))
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.stats import ks_2samp, wasserstein_distance
from scipy.special import logsumexp
from scipy.optimize import minimize_scalar
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from numba import njit
from janus_connectivity.model import generate_points,generate_orientations
from janus_connectivity.graph_build import build_directed_graph,prepare_cells
from janus_connectivity.feedback_order import classify_graph
from concurrent.futures import ProcessPoolExecutor,as_completed

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def edge_hash(g):return hashlib.sha256(g.indptr.tobytes()+g.indices.tobytes()).hexdigest()
def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temp.replace(path)
def scc_mass(g):
    _,labels=connected_components(g,directed=True,connection='strong')
    return int(np.bincount(labels).max())
def replay(row):
    p=generate_points(int(row['N']),float(row['L']),int(row['geometry_seed']))
    th=generate_orientations(len(p),float(row['W']),int(row['disorder_seed']))
    cells=prepare_cells(p)
    g,_,_=build_directed_graph(p,th,int(row['q']),cells=cells)
    assert edge_hash(g)==row['edge_hash'],'Frozen base hash mismatch: '+row['graph_id']
    return p,th,g,cells
def values(g,p,theta):
    src=np.repeat(np.arange(len(p)),np.diff(g.indptr))
    dx=p[g.indices,0]-p[src,0];dy=p[g.indices,1]-p[src,1]
    ang=np.arctan2(dx,-dy)-theta[src]
    return [np.hypot(dx,dy),np.arctan2(np.sin(ang),np.cos(ang)),dy]
def means(v):
    length,angle,dy=v
    return {'A_loc':float(np.cos(angle).mean()),'D':float((-dy/length).mean()),'f_up':float((dy>0).mean())}
def rate(dy):
    h=np.histogram(dy,bins=2048,range=(-2,2))[0]
    x=np.linspace(-2+2/2048,2-2/2048,2048);ok=h>0
    lw=np.log(h[ok]/h.sum())
    opt=minimize_scalar(lambda t:logsumexp(lw+t*x[ok]),bounds=(0,32),method='bounded')
    return float(-opt.fun)
@njit(cache=True)
def changed(ptr,a,b):
    common=0
    for v in range(len(ptr)-1):
        for i in range(ptr[v],ptr[v+1]):
            for j in range(ptr[v],ptr[v+1]):
                if a[i]==b[j]:common+=1;break
    return 1-common/len(a)
def arrange(p,th,ell,rep):
    noise=np.random.default_rng(np.random.SeedSequence([2026091973,rep])).normal(size=(256,256))
    field=gaussian_filter(noise,sigma=ell/np.sqrt(2),mode='wrap')
    v=map_coordinates(field,[p[:,0]+64,p[:,1]+64],order=1,mode='wrap')
    order=np.argsort(v,kind='stable');out=np.empty_like(th);out[order]=np.sort(th)
    assert np.array_equal(np.sort(out),np.sort(th))
    return out
def accepts_fidelity(rec):
    return bool(rec['outdegree_exact'] and rec['orientation_multiset_exact']
        and rec['changed_edge_fraction']>=.25 and rec['I0_relative_change']<=.03
        and all(rec[x+'_KS']<=.03 and rec[x+'_Wasserstein']<=.03
                for x in ('length','relative_angle','dy')))

def ell_job(row):
    start=time.perf_counter();p,th,g,cells=replay(row)
    original=values(g,p,th);m0=means(original);r0=rate(original[2]);s0=scc_mass(g)/len(p)
    assert np.isclose(s0,float(row['S_original']),rtol=0,atol=2e-15)
    out=[]
    for ell in json.loads((INPUT/'02_ell_sensitivity/PROTOCOL.json').read_text())['ell_over_R']:
        name=f"W{float(row['W']):.3f}_r{int(row['realization']):02d}_ell{ell}"
        done=ROOT/'02_ell_sensitivity/records'/f'{name}.json'
        if done.exists():out.append(json.loads(done.read_text()));continue
        target=arrange(p,th,ell,int(row['realization']))
        h,_,_=build_directed_graph(p,target,int(row['q']),cells=cells)
        vh=values(h,p,target);m1=means(vh);r1=rate(vh[2]);delta=abs(r1-r0)/max(r0,.01)
        cf=float(changed(g.indptr,g.indices,h.indices))
        rec={k:row[k] for k in ('graph_id','edge_hash','geometry_seed','disorder_seed','W','realization','L','q','N')}
        rec.update(ell_over_R=ell,field_seed=[2026091973,int(row['realization'])],
                   surrogate_edge_hash=edge_hash(h),changed_edge_fraction=cf,I0_original=r0,I0_reassigned=r1,
                   I0_relative_change=delta,outdegree_exact=bool(np.array_equal(np.diff(g.indptr),np.diff(h.indptr))),
                   orientation_multiset_exact=True,S_original=s0)
        for key in m0:rec[key+'_original']=m0[key];rec[key+'_reassigned']=m1[key]
        for label,a,b,norm in zip(('length','relative_angle','dy'),original,vh,(2,np.pi,2)):
            rec[label+'_KS']=float(ks_2samp(a,b,method='asymp').statistic)
            rec[label+'_Wasserstein']=float(wasserstein_distance(a,b)/norm)
        for quant in (.1,.25,.5,.75,.9,.95):
            qa=float(np.quantile(original[0],quant));qb=float(np.quantile(vh[0],quant))
            rec[f'length_q{quant:g}_original']=qa;rec[f'length_q{quant:g}_reassigned']=qb
            rec[f'length_q{quant:g}_change']=qb-qa
        rec['accepted']=accepts_fidelity(rec)
        # All fidelity metrics are recorded on disk before the surrogate SCC is evaluated.
        dump(ROOT/'02_ell_sensitivity/fidelity'/f'{name}.json',rec)
        rec['fidelity_sha256']=sha(ROOT/'02_ell_sensitivity/fidelity'/f'{name}.json')
        rec['S_reassigned']=scc_mass(h)/len(p) if rec['accepted'] else None
        rec['log10_S_ratio']=float(np.log10(rec['S_reassigned']/s0)) if rec['accepted'] else None
        rec['status']='ACCEPTED' if rec['accepted'] else 'NO_ACCEPTED_INTERVENTION_UNDER_COMMON_FIDELITY_GATE'
        if ell==4:
            assert rec['surrogate_edge_hash']==row['ell4_edge_hash'],'Historical ell4 hash mismatch'
            assert np.isclose(rec['S_reassigned'],float(row['S_ell4']),rtol=0,atol=2e-15)
        rec['elapsed_seconds']=time.perf_counter()-start
        dump(done,rec);out.append(rec)
    return {'graph_id':row['graph_id'],'ell_records':len(out),'seconds':time.perf_counter()-start}
def filtered(g,ix,classification,k):
    keep=np.ones(g.nnz,dtype=bool);keep[ix]=classification<=k
    prefix=np.r_[0,np.cumsum(keep,dtype=np.int32)]
    return csr_matrix((g.data[keep],g.indices[keep],prefix[g.indptr]),shape=g.shape)
def deletion_job(row):
    start=time.perf_counter();done=ROOT/'03_boundary_deletion_ensemble/records'/f"{row['edge_hash']}.json"
    if done.exists():return {'graph_id':row['graph_id'],'reused':True}
    p,th,g,cells=replay(row)
    archive=PROJECT/row['classification_archive'] if row['classification_archive'] else None
    if archive is not None:
        assert sha(archive)==row['classification_sha256']
        with np.load(archive,allow_pickle=False) as a:
            ix=a['feedback_csr_edge_index'];code=a['classification']
    else:
        derived=classify_graph(g,p[:,1]);ix=derived['up_ix'];code=derived['labels']
        for k in (1,2):
            assert np.isclose(scc_mass(filtered(g,ix,code,k))/len(p),float(row[f'prior_S{k}_full']),atol=2e-15,rtol=0)
        archive=ROOT/'03_boundary_deletion_ensemble/classifications'/f"{row['edge_hash']}.npz"
        archive.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(archive,feedback_csr_edge_index=ix,classification=code)
        del derived
    src=np.repeat(np.arange(len(p)),np.diff(g.indptr))
    assert np.array_equal(ix,np.flatnonzero(p[g.indices,1]>p[src,1]))
    result=[]
    for width in (2,4,8):
        core=np.flatnonzero(np.all((p>width)&(p<float(row['L'])-width),axis=1))
        rec={key:row[key] for key in ('graph_id','edge_hash','W','L','q','N','geometry_seed','disorder_seed')}
        rec.update(width_over_R=width,N_retained=int(len(core)),retained_fraction=len(core)/len(p))
        rec['S_full_core']=scc_mass(g[core,:][:,core])/len(core)
        for k in (1,2):
            h=filtered(g,ix,code,k);rec[f'S{k}_core']=scc_mass(h[core,:][:,core])/len(core)
            if width==4 and pd.notna(row[f'prior_S{k}_4R']):
                assert np.isclose(rec[f'S{k}_core'],float(row[f'prior_S{k}_4R']),atol=2e-15,rtol=0)
        assert rec['S1_core']<=rec['S2_core']<=rec['S_full_core']
        rec['delta_full_G2']=rec['S_full_core']-rec['S2_core']
        rec['ratio_full_G2']=rec['S_full_core']/rec['S2_core']
        rec['replay_hash_pass']=True;rec['classification_sha256']=sha(archive)
        result.append(rec)
    dump(done,{'rows':result,'elapsed_seconds':time.perf_counter()-start})
    return {'graph_id':row['graph_id'],'seconds':time.perf_counter()-start}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('task',choices=['ell','deletion']);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--graph-id',help='Replay only this frozen manifest identity');args=ap.parse_args()
    path=INPUT/('02_ell_sensitivity/BASE_MANIFEST.csv' if args.task=='ell' else '03_boundary_deletion_ensemble/BASE_MANIFEST.csv')
    rows=pd.read_csv(path,dtype={'geometry_seed':str,'disorder_seed':str}).to_dict('records')
    for row in rows:
        if args.task=='deletion':row['classification_archive']=''
    if args.graph_id:
        rows=[row for row in rows if row['graph_id']==args.graph_id]
        if len(rows)!=1:raise ValueError('Graph identity must occur exactly once in the fixed manifest')
    fn=ell_job if args.task=='ell' else deletion_job
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i,future in enumerate(as_completed([pool.submit(fn,row) for row in rows]),1):
            print(i,len(rows),future.result(),flush=True)
    records=[]
    folder=ROOT/('02_ell_sensitivity' if args.task=='ell' else '03_boundary_deletion_ensemble')
    for p in sorted((folder/'records').glob('*.json')):
        r=json.loads(p.read_text());records.extend(r['rows'] if args.task=='deletion' else [r])
    pd.DataFrame(records).to_csv(folder/'GRAPHLEVEL.csv',index=False)
if __name__=='__main__':main()
