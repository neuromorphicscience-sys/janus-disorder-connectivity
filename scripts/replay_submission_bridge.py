"""Only the 48 frozen intervention pairs; predeclared observables, no search."""
from pathlib import Path
import os,sys,json,hashlib,time,argparse
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[name]='1'
R=Path(__file__).resolve().parents[1];P=R;F=R/'outputs/final_bridge'
sys.path.insert(0,str(R/'src'))
import numpy as np,pandas as pd
from scipy.ndimage import gaussian_filter,map_coordinates
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from janus_connectivity.model import generate_points,generate_orientations
from janus_connectivity.graph_build import build_directed_graph,prepare_cells
from janus_connectivity.feedback_order import classify_graph,make_filtered_graph
from janus_connectivity.recurrence import internal_recurrence,strict_core_mask
from concurrent.futures import ProcessPoolExecutor,as_completed

def edge_hash(g):return hashlib.sha256(g.indptr.tobytes()+g.indices.tobytes()).hexdigest()
def s(g):
    _,labels=connected_components(g,directed=True,connection='strong')
    return np.bincount(labels).max()/g.shape[0]
def correlation(p,angles):
    # Same prespecified anchor sample for original/reassigned and all W.
    anchors=np.sort(np.random.default_rng(2026100801).choice(len(p),512,replace=False))
    tree=cKDTree(p);neigh=tree.query_ball_point(p[anchors],8)
    out=[]
    for state,theta in angles.items():
        pol2=abs(np.exp(1j*theta).mean())**2
        perm=(len(p)*pol2-1)/(len(p)-1) # exact E over all permutations, i!=j
        counts=np.zeros(4,dtype=np.int64);total=np.zeros(4)
        for i,nb in zip(anchors,neigh):
            j=np.asarray(nb,dtype=np.int64);j=j[j!=i]
            d=np.linalg.norm(p[j]-p[i],axis=1)
            bins=np.searchsorted([1,2,4,8],d,side='left');valid=bins<4
            counts+=np.bincount(bins[valid],minlength=4)
            total+=np.bincount(bins[valid],weights=np.cos(theta[i]-theta[j[valid]]),minlength=4)
        for b,(lo,hi) in enumerate(zip([0,1,2,4],[1,2,4,8])):
            raw=total[b]/counts[b]
            out.append(dict(state=state,r_min_R=lo,r_max_R=hi,anchor_n=len(anchors),ordered_pair_n=int(counts[b]),polarization_squared=pol2,exact_permutation_baseline=perm,raw_cos_mean=raw,C_theta_conn=raw-pol2,C_theta_perm_excess=raw-perm))
    return out
def job(row):
    dest=F/'01_intervention_bridge/records'/f"{row['graph_id']}.json"
    if dest.exists():return row['graph_id'],'cached'
    start=time.perf_counter();n=int(row['N']);L=float(row['L']);q=int(row['q'])
    p=generate_points(n,L,int(row['geometry_seed']));th=generate_orientations(n,float(row['W']),int(row['disorder_seed']))
    noise=np.random.default_rng(np.random.SeedSequence([2026091973,int(row['realization'])])).normal(size=(256,256))
    field=gaussian_filter(noise,sigma=4/np.sqrt(2),mode='wrap')
    v=map_coordinates(field,[p[:,0]+64,p[:,1]+64],order=1,mode='wrap')
    th_r=np.empty_like(th);th_r[np.argsort(v,kind='stable')]=np.sort(th)
    assert np.array_equal(np.sort(th),np.sort(th_r))
    cells=prepare_cells(p);core=strict_core_mask(p,L,4)
    rec={k:row[k] for k in ('graph_id','W','realization','L','N','q','geometry_seed','disorder_seed','edge_hash','ell4_edge_hash')}
    rec['original_graph_id']=row['graph_id'];rec['reassigned_graph_id']=row['graph_id']+'_ell4_reassigned'
    hist=[];outdegrees=[]
    for state,theta,expected_s,expected_hash in [('original',th,row['S_original'],row['edge_hash']),('reassigned',th_r,row['S_ell4'],row['ell4_edge_hash'])]:
        g,_,_=build_directed_graph(p,theta,q,cells=cells)
        assert edge_hash(g)==expected_hash,(row['graph_id'],state,'hash mismatch')
        rec['S_'+state]=float(s(g));assert np.isclose(rec['S_'+state],float(expected_s),rtol=0,atol=2e-15)
        outdegrees.append(np.diff(g.indptr));indeg=np.bincount(g.indices,minlength=n)
        rec['mean_indegree_'+state]=float(indeg.mean());rec['variance_indegree_'+state]=float(indeg.var())
        rec['zero_indegree_fraction_'+state]=float((indeg==0).mean())
        rec['reciprocal_edge_fraction_'+state]=float(g.multiply(g.T).sum()/g.nnz)
        _,weak=connected_components(g,directed=True,connection='weak');rec['WCC_'+state]=float(np.bincount(weak).max()/n)
        for d,c in enumerate(np.bincount(indeg)):
            hist.append(dict(state=state,indegree=d,node_count=int(c),fraction=c/n))
        print(row['graph_id'],state,'classifying',flush=True)
        cls=classify_graph(g,p[:,1])
        archive=F/'01_intervention_bridge/classifications'/f'{expected_hash}.npz'
        np.savez_compressed(archive,feedback_csr_edge_index=cls['up_ix'],classification=cls['labels'])
        rec['classification_sha256_'+state]=hashlib.sha256(archive.read_bytes()).hexdigest()
        for k in (1,2):
            h=make_filtered_graph(cls,include_order2=k==2);rec[f'S{k}_{state}']=float(s(h))
            m=internal_recurrence(h,core)
            for key,target in [('f_internal_rec','f_int_rec'),('eta_internal','eta_int'),('largest_fraction_of_core','S_int'),('recurrent_mass','recurrent_mass'),('n_core','N_core')]:
                rec[f'{target}_{k}_{state}']=None if key=='eta_internal' and m['recurrent_mass']==0 else m[key]
            rec[f'eta_int_zero_convention_{k}_{state}']=m['eta_internal']
            assert abs(m['identity_residual'])<1e-14
        assert rec[f'S1_{state}']<=rec[f'S2_{state}']<=rec[f'S_{state}']
        del g,cls,h
    assert np.array_equal(*outdegrees)
    rec['outdegree_exact']=True;rec['orientation_multiset_exact']=True;rec['both_replay_hashes_verified']=True
    for key in ['S','S1','S2','f_int_rec_1','f_int_rec_2','eta_int_1','eta_int_2','S_int_1','S_int_2','mean_indegree','variance_indegree','zero_indegree_fraction','reciprocal_edge_fraction','WCC']:
        a=rec[key+'_original'];b=rec[key+'_reassigned'];valid=a is not None and b is not None
        rec[key+'_difference']=b-a if valid else None
        rec[key+'_ratio']=b/a if valid and a>0 else None
        rec['log10_'+key+'_ratio']=float(np.log10(b/a)) if valid and a>0 and b>0 else None
    corr=correlation(p,dict(original=th,reassigned=th_r))
    for rows in (hist,corr):
        for v in rows:v.update(graph_id=row['graph_id'],W=row['W'],realization=row['realization'])
    dest.write_text(json.dumps(dict(bridge=rec,indegree=hist,correlation=corr,seconds=time.perf_counter()-start),indent=2)+'\n')
    return row['graph_id'],time.perf_counter()-start
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=3);ap.add_argument('--limit',type=int);a=ap.parse_args()
    manifest=R/'data/final_bridge/INTERVENTION_PAIR_MANIFEST.csv'
    rows=pd.read_csv(manifest,dtype={'geometry_seed':str,'disorder_seed':str}).to_dict('records')
    assert len(rows)==48
    (F/'01_intervention_bridge/records').mkdir(parents=True,exist_ok=True);(F/'01_intervention_bridge/classifications').mkdir(exist_ok=True)
    if a.limit:rows=rows[:a.limit]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(job,row) for row in rows]):print('DONE',future.result(),flush=True)
    objects=[json.loads(p.read_text()) for p in sorted((F/'01_intervention_bridge/records').glob('*.json'))]
    for field,name in [('bridge','INTERVENTION_BRIDGE_SOURCE'),('indegree','INDEGREE_DISTRIBUTIONS'),('correlation','ORIENTATION_CORRELATION_SOURCE')]:
        vals=[]
        for obj in objects:vals.extend([obj[field]] if field=='bridge' else obj[field])
        pd.DataFrame(vals).to_csv(F/f'01_intervention_bridge/{name}.csv',index=False)
if __name__=='__main__':main()
