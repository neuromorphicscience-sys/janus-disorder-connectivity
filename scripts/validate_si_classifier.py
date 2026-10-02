"""Independent 0/1-cost oracle comparison on five unchanged archived toy graphs."""
from pathlib import Path
from collections import deque
import argparse,json,sys
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
PROJECT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(PROJECT/'src'))
from janus_connectivity.feedback_order import classify_graph,make_filtered_graph
ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
nodes=pd.read_csv(a.data/'classifier_toy_nodes.csv',engine='python');edges=pd.read_csv(a.data/'classifier_toy_edges.csv',engine='python');out=[]
for name,nodeset in nodes.groupby('toy'):
 y=nodeset.sort_values('node').progress_y.to_numpy(float);e=edges[edges.toy.eq(name)];pairs=list(zip(e.source,e.target));G=csr_matrix((np.ones(len(e),np.uint8),(e.source,e.target)),shape=(len(y),len(y)));G.sort_indices();r=classify_graph(G,y)
 adj={i:[] for i in range(len(y))}
 for u,v in pairs:adj[u].append(v)
 oracle=[]
 for u,v in zip(r['u'],r['v']):
  q=deque([(int(v),0)]);visited={(int(v),0)};best=2
  while q:
   x,c=q.popleft()
   if x==u:best=c;break
   for z in adj[x]:
    if x==u and z==v:continue
    weight=int(y[z]>y[x]);state=(z,c+weight)
    if state[1]<=1 and state not in visited:
     visited.add(state)
     if weight:q.append(state)
     else:q.appendleft(state)
  oracle.append(best+1 if best<2 else 3)
 assert np.array_equal(r['labels'],oracle),(name,r['labels'],oracle)
 e1=set(zip(*make_filtered_graph(r,include_order2=False).nonzero()));e2=set(zip(*make_filtered_graph(r,include_order2=True).nonzero()));assert e1<=e2<=set(pairs)
 out.append({'toy':name,'feedback_edges':len(oracle),'oracle_matches':True,'edge_nesting':True})
pd.DataFrame(out).to_csv(a.output/'toy_oracle_recheck.csv',index=False);print(json.dumps(out,indent=2))
