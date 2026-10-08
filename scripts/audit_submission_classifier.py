"""Reproduce the exact historical subsets and independently audit all edges."""
from pathlib import Path
import os,sys,json,hashlib
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[name]='1'
P=Path(__file__).resolve().parents[1];F=P/'outputs/final_bridge';sys.path.insert(0,str(P/'src'))
import numpy as np,pandas as pd
from scipy.sparse.csgraph import connected_components
from collections import deque
from janus_connectivity.model import generate_points,generate_orientations
from janus_connectivity.graph_build import build_directed_graph
from janus_connectivity.feedback_order import classify_graph
V=P/'data/final_bridge'
def oracle(g,y):
    src=np.repeat(np.arange(len(y)),np.diff(g.indptr));isf=y[g.indices]>y[src];answers=np.zeros(g.nnz,dtype=np.uint8)
    for e in np.flatnonzero(isf):
        target=int(src[e]);start=int(g.indices[e]);dist=np.full(len(y),3,dtype=np.int8);dist[start]=0;q=deque([start])
        while q:
            u=q.popleft()
            for j in range(g.indptr[u],g.indptr[u+1]):
                if j==e:continue
                v=int(g.indices[j]);c=int(dist[u])+int(isf[j])
                if c<=1 and c<dist[v]:
                    dist[v]=c
                    if isf[j]:q.append(v)
                    else:q.appendleft(v)
        answers[e]=dist[target]+1 if dist[target]<=1 else 3
    return answers
def check(g,y,name,parent,rule,indices):
    r=classify_graph(g,y);actual=np.zeros(g.nnz,dtype=np.uint8);actual[r['up_ix']]=r['labels'];want=oracle(g,y)
    agreement=np.array_equal(actual,want)
    pd.DataFrame(dict(source=np.repeat(np.arange(len(y)),np.diff(g.indptr)),target=g.indices,source_y=y[np.repeat(np.arange(len(y)),np.diff(g.indptr))],target_y=y[g.indices],production_label=actual,oracle_label=want)).to_csv(F/f'02_classifier_audit/edges/{name}.csv',index=False)
    _,lab=connected_components(g,directed=True,connection='strong');sizes=np.bincount(lab)
    out=dict(name=name,parent_graph_id=parent['graph_id'],parent_hash=parent['expected_hash'],selection_rule=rule,N=len(y),edges=g.nnz,D0_edges=int((actual==0).sum()),feedback_F_edges=len(r['up_ix']),order1_edges=int((actual==1).sum()),order2_edges=int((actual==2).sum()),gt2_or_infinite_edges=int((actual==3).sum()),nontrivial_SCCs=int((sizes>=2).sum()),edge_by_edge_agreement=agreement,vertex_index_sha256=hashlib.sha256(indices.tobytes()).hexdigest())
    np.save(F/f'02_classifier_audit/vertices/{name}.npy',indices)
    if not agreement:raise AssertionError(out)
    return out
def main():
    for name in ['edges','vertices']:(F/f'02_classifier_audit/{name}').mkdir(parents=True,exist_ok=True)
    pop=json.loads((V/'frozen_replay_inputs.json').read_text());old=json.loads((V/'historical_classifier_baseline.json').read_text());baseline={x['name']:x for x in old['historical_induced_subgraphs']}
    records=[]
    # Decision rule fixed before comparisons: retain historical set if >=20 positive
    # edges at each order pooled and >=1 positive order-2 case at each parent.
    need_rich=(sum(x['F1_count'] for x in baseline.values())<20 or sum(x['F2_exact_count'] for x in baseline.values())<20 or any(not any(x['F2_exact_count']>0 for x in baseline.values() if x['parent_graph_id']==p['graph_id']) for p in pop['benchmark']))
    (F/'02_classifier_audit/PREDECLARED_SELECTION.json').write_text(json.dumps(dict(historical_rule='2 uniform without replacement; 6 weak-graph BFS from random vertex, shuffled neighbors; RNG seeded by first 16 hex hash digits',coverage_rule='>=20 order1 and >=20 order2 edges total and at least one order2-positive subset per historical parent',feedback_rich_needed=need_rich,rich_rule='if needed: 8 evenly indexed archived upward edges per parent, both endpoints plus nearest Euclidean vertices to midpoint until N=128; choices fixed before oracle comparison'),indent=2)+'\n')
    for parent in pop['benchmark']:
        n=round(parent['sigma']*parent['L']**2);p=generate_points(n,parent['L'],int(parent['geometry_seed']));th=generate_orientations(n,parent['control'],int(parent['disorder_seed']));g,_,_=build_directed_graph(p,th,parent['q'])
        assert hashlib.sha256(g.indptr.tobytes()+g.indices.tobytes()).hexdigest()==parent['expected_hash']
        y=p[:,1];rng=np.random.default_rng(int(parent['expected_hash'][:16],16));weak=(g+g.T).tocsr()
        for j in range(8):
            if j<2:ix=rng.choice(n,64,replace=False);rule='uniform_random_64'
            else:
                start=int(rng.integers(n));seen={start};queue=[start];pos=0
                while pos<len(queue) and len(seen)<64:
                    x=queue[pos];pos+=1;nb=weak.indices[weak.indptr[x]:weak.indptr[x+1]].copy();rng.shuffle(nb)
                    for z in nb:
                        if int(z) not in seen:seen.add(int(z));queue.append(int(z))
                        if len(seen)==64:break
                ix=np.array(sorted(seen),dtype=np.int64);rule='weak_neighborhood_BFS_64'
            name=f"{parent['graph_id']}_induced{j}";r=check(g[ix][:,ix].tocsr(),y[ix],name,parent,rule,ix);b=baseline[name]
            assert (r['order1_edges'],r['order2_edges'],r['feedback_F_edges'])==(b['F1_count'],b['F2_exact_count'],b['F_count'])
            records.append(r)
        if need_rich:
            src=np.repeat(np.arange(n),np.diff(g.indptr));fi=np.flatnonzero(y[g.indices]>y[src])
            for j,e in enumerate(fi[np.linspace(0,len(fi)-1,8,dtype=int)]):
                midpoint=(p[src[e]]+p[g.indices[e]])/2;distance=((p-midpoint)**2).sum(axis=1)
                endpoints=[int(src[e]),int(g.indices[e])];near=np.argsort(distance,kind='stable')[:130]
                ix=np.array(endpoints+[int(v) for v in near if int(v) not in endpoints][:126],dtype=np.int64);ix=np.sort(ix)
                records.append(check(g[ix][:,ix].tocsr(),y[ix],f"{parent['graph_id']}_feedbacklocal{j}",parent,'feedback_edge_midpoint_nearest_128',ix))
        pd.DataFrame(records).to_csv(F/'02_classifier_audit/CLASSIFIER_VALIDATION_COVERAGE.csv',index=False)
        print(parent['graph_id'],'audited',flush=True)
    d=pd.DataFrame(records);lines=['# Classifier validation coverage','', 'Exact historical RNG and weak-neighborhood sampling were reproduced from the original code. Labels were compared on every induced edge: D0=0, F1=1, F2=2, >2/infinity=3. The independent oracle uses 0/1 shortest-path searches excluding the tested feedback edge. No graph was selected by disagreement.','']
    for rule,group in d.groupby('selection_rule'):
        lines.append(f'## {rule}: {len(group)} subgraphs')
        for col in ['N','edges','D0_edges','feedback_F_edges','order1_edges','order2_edges','gt2_or_infinite_edges','nontrivial_SCCs']:
            lines.append(f'- {col}: min/median/max={group[col].min():g}/{group[col].median():g}/{group[col].max():g}; total={group[col].sum()}')
    lines+=['',f"All {len(d)} subgraphs agree edge by edge. Total positive labels: order1={d.order1_edges.sum()}, order2={d.order2_edges.sum()}.",f'Feedback-rich supplement required by predeclared coverage rule: {need_rich}.','The five archived toys are retained; their edge-wise independent reproduction is run with the existing public validation script.']
    (F/'02_classifier_audit/CLASSIFIER_VALIDATION_REPORT.md').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':main()
