"""Reproduce descriptive crossing intervals from frozen shared-stream data."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=ROOT/'data/final_bridge');ap.add_argument('--output',type=Path,default=ROOT/'outputs/final_bridge/05_finite_size_uncertainty');args=ap.parse_args()
F=args.output.parent;P=ROOT;args.output.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(args.data/'FINITE_SIZE_STREAM_SOURCE.csv',dtype={'stream_id':str});assert len(d)==4634 and not d.duplicated(['L','W','stream_id']).any()
def crossing(w,s,level):
    ok=np.isfinite(s);w=w[ok];s=s[ok]
    eq=np.flatnonzero(s==level)
    if len(eq)==1:return float(w[eq[0]])
    if len(eq)>1:return np.nan
    ix=np.flatnonzero((s[:-1]<level)&(s[1:]>level))
    if len(ix)!=1:return np.nan # retain ambiguity, no monotonic smoothing
    i=ix[0];return float(w[i]+(level-s[i])*(w[i+1]-w[i])/(s[i+1]-s[i]))
def weighted_median(values,weights):
    o=np.argsort(values);v=values[o];a=weights[o];c=np.cumsum(a);n=int(c[-1])
    if not n:return np.nan
    lo=v[np.searchsorted(c,(n-1)//2+1)];hi=v[np.searchsorted(c,n//2+1)]
    return (lo+hi)/2
summary=[];replicates=[];counts=[]
B=4000
for L,g in d.groupby('L'):
    streams=sorted(g.stream_id.unique());index={s:i for i,s in enumerate(streams)};groups=[];ws=[]
    for W,z in g.groupby('W'):
        ws.append(W);groups.append((z.S.to_numpy(),np.array([index[s] for s in z.stream_id])))
        counts.append(dict(L=L,W=W,n=len(z),unique_streams_total=len(streams)))
    ws=np.array(ws);med=np.array([np.median(v) for v,i in groups]);point=[crossing(ws,med,p) for p in [.1,.5,.9]]
    rng=np.random.default_rng(np.random.SeedSequence([2026100802,int(L)]));values=[]
    for b in range(B):
        weight=np.bincount(rng.integers(len(streams),size=len(streams)),minlength=len(streams))
        sm=np.array([weighted_median(v,weight[i]) for v,i in groups]);c=[crossing(ws,sm,p) for p in [.1,.5,.9]]
        width=c[2]-c[0] if np.isfinite(c[0]+c[2]) else np.nan
        values.append((c[1],width));replicates.append(dict(L=L,replicate=b,W10=c[0],W_half=c[1],W90=c[2],Delta_W10_90=width))
    arr=np.array(values)
    for j,(metric,central) in enumerate([('W_half',point[1]),('Delta_W10_90',point[2]-point[0])]):
        a=arr[:,j];a=a[np.isfinite(a)];q=np.quantile(a,[.025,.16,.5,.84,.975]);summary.append(dict(L=L,metric=metric,central=central,CI95_low=q[0],CI68_low=q[1],bootstrap_median=q[2],CI68_high=q[3],CI95_high=q[4],successful_fraction=len(a)/B,bootstrap_n=B,unique_stream_n=len(streams),display_decimals=3))
    print('bootstrap',L,'streams',len(streams),flush=True)
pd.DataFrame(summary).to_csv(args.output/'FINITE_SIZE_BOOTSTRAP_SUMMARY.csv',index=False)
pd.DataFrame(replicates).to_csv(args.output/'FINITE_SIZE_BOOTSTRAP_REPLICATES.csv',index=False)
pd.DataFrame(counts).to_csv(args.output/'CONDITION_COUNTS.csv',index=False)
