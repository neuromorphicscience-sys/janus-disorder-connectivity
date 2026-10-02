"""Reproduce nine supplemental figures from explicit archived CSV source tables."""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

MAG='#D81B60'; G1='#004D40'; G2='#1E88E5'; AMBER='#FFC107'; GRAY='#777777'
MARKERS=['o','s','^','D']; STYLES=['-','--','-.',':']
plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans','font.size':8,
 'axes.labelsize':8.5,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,
 'axes.linewidth':.7,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
 'axes.spines.top':True,'axes.spines.right':True,'pdf.fonttype':42,'svg.fonttype':'none',
 'lines.linewidth':1.1,'lines.markersize':3.5,'legend.frameon':False,'text.color':'#111111','axes.labelcolor':'#111111'})
ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
DATA=args.data; OUT=args.output;OUT.mkdir(parents=True,exist_ok=True)
def read(name):return pd.read_csv(DATA/(name+'.csv'),engine='python')
def canvas(rows=2,cols=2,h=110):
 fig,axs=plt.subplots(rows,cols,figsize=(178/25.4,h/25.4),layout='constrained',squeeze=False)
 for i,ax in enumerate(axs.flat):ax.text(-.10,1.04,f'({chr(97+i)})',transform=ax.transAxes,size=10,weight='bold');ax.tick_params(top=True,right=True)
 return fig,axs
def save(fig,n):
 for ext in ['pdf','svg']:fig.savefig(OUT/f'figS{n}.{ext}')
 fig.savefig(OUT/f'figS{n}.png',dpi=600);plt.close(fig);print('Fig. S'+str(n),flush=True)
def err(ax,x,a,color,label=None,marker='o',ls='-',offset=0):
 m=a.median();lo,hi=a.quantile([.25,.75]);ax.errorbar(x+offset,m,yerr=[[m-lo],[hi-m]],color=color,marker=marker,ls=ls,capsize=2,label=label,zorder=5)
def summary(ax,d,x,y='S',color=MAG,label=None,marker='o',ls='-'):
 g=d.groupby(x)[y];m=g.median();lo=g.quantile(.25);hi=g.quantile(.75)
 ax.plot(m.index,m,color=color,marker=marker,ls=ls,label=label)
 ax.fill_between(m.index.to_numpy(float),lo.to_numpy(float),hi.to_numpy(float),color=color,alpha=.12,lw=0)
def ecdf(ax,a,color,label,ls='-'):
 x=np.sort(np.asarray(a));ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',color=color,ls=ls,label=label)
def size_axis(ax):ax.set_xscale('log',base=2);ax.set_xticks([64,128,256,512],labels=['64','128','256','512'])

def s1():
 d=read('source_figS1');fig,axs=canvas(1,2,78)
 for j,q in enumerate([4,5,7,10]):summary(axs[0,0],d[d.q.eq(q)],'control',label=rf'$q={q}$',marker=MARKERS[j],ls=STYLES[j])
 axs[0,0].set(xlabel=r'$W$',ylabel=r'$S$',ylim=(-.02,1.03));axs[0,0].legend(loc='lower left',ncol=2)
 endpoint=[]
 for j,(q,g) in enumerate(d.groupby('q')):
  w=g.control.max();a=g[g.control.eq(w)].S;err(axs[0,1],q,a,MAG);endpoint.append(dict(q=int(q),W=w,n=len(a),median=a.median(),Q1=a.quantile(.25),Q3=a.quantile(.75)))
 axs[0,1].set(xlabel=r'$q$',ylabel=r'$S$ at largest sampled $W$',ylim=(0,1.04),xticks=[4,5,6,7,8,10,12])
 axs[0,1].text(.05,.13,'Different endpoint W for each q',transform=axs[0,1].transAxes,fontsize=7.5)
 pd.DataFrame(endpoint).to_csv(OUT/'q_endpoint_summary.csv',index=False);save(fig,1)
def s2():
 d=read('source_figS2');fig,axs=canvas(3,2,155)
 for ax,(dist,xlab,name) in zip(axs.flat,[('gaussian','W','Gaussian'),('uniform','A','Bounded uniform'),('vonmises',r'\kappa','von Mises'),('laplace','b','Wrapped Laplace')]):
  g=d[d.distribution.eq(dist)&d.correlation_family.eq('iid')];summary(ax,g,'control')
  ax.set(xlabel='$'+xlab+'$',ylabel=r'$S$',ylim=(-.02,1.03));ax.text(.97,.83,name,transform=ax.transAxes,ha='right')
 for ax,anis in zip(axs.flat[4:],[False,True]):
  g=d[d.control_name.eq('alpha')];g=g[g.correlation_family.str.startswith('ANISO') if anis else ~g.correlation_family.str.startswith('ANISO')]
  for j,((ex,ey),a) in enumerate(g.groupby(['ell_x','ell_y'])):summary(ax,a,'control',label=rf'$({ex:g},{ey:g})R$',marker=MARKERS[j],ls=STYLES[j])
  ax.set(xlabel=r'$\alpha$',ylabel=r'$S$',ylim=(-.02,1.03));ax.legend(title=r'$(\ell_x,\ell_y)$',loc='upper right',ncol=2 if not anis else 1)
 save(fig,2)
def s3():
 d=read('source_figS3');fig,axs=canvas(2,2,113);cross=[]
 for j,(L,g) in enumerate(d.groupby('L')):
  ax=axs[0,0];ax.scatter(g.W,g.S,s=4,color=MAG,alpha=.18,marker=MARKERS[j],linewidths=0)
  summary(ax,g,'W',marker=MARKERS[j],ls=STYLES[j])
  summary(axs[0,1],g,'W',label=rf'${L:g}$',marker=MARKERS[j],ls=STYLES[j])
  m=g.groupby('W').S.median().sort_index();r=dict(L=L,n=len(g),n_W=len(m))
  for t in [.1,.5,.9]:
   brackets=[(x1,x2,y1,y2) for x1,x2,y1,y2 in zip(m.index[:-1],m.index[1:],m.iloc[:-1],m.iloc[1:]) if y1<t<=y2]
   assert len(brackets)==1,('Crossing not uniquely bracketed',L,t,brackets)
   x1,x2,y1,y2=brackets[0];r[f'W_{t:g}']=x1+(t-y1)*(x2-x1)/(y2-y1);r[f'bracket_{t:g}']=f'{x1:g},{x2:g}'
  r['width_10_90']=r['W_0.9']-r['W_0.1'];cross.append(r)
 axs[0,0].set(xlabel=r'$W$',ylabel=r'Individual $S$ and median',ylim=(-.02,1.03))
 axs[0,1].set(xlabel=r'$W$',ylabel=r'$S$',xlim=(.75,.835),ylim=(-.02,1.03));axs[0,1].legend(title=r'$L/R$',ncol=2,loc='lower right')
 c=pd.DataFrame(cross);c.to_csv(OUT/'finite_size_crossings.csv',index=False)
 axs[1,0].plot(c.L,c['W_0.5'],color=MAG,marker='o');axs[1,0].set(xlabel=r'$L/R$',ylabel=r'$W_{1/2}$')
 axs[1,1].plot(c.L,c.width_10_90,color=MAG,marker='o');axs[1,1].set(xlabel=r'$L/R$',ylabel=r'$\Delta W_{10-90}$')
 for ax in axs[1]:size_axis(ax)
 save(fig,3)
def paired_diag(ax,d,key):
 x=np.arange(len(d))+1
 for name,marker,color,offset in [('KS','o',GRAY,-.12),('Wasserstein','s','#333333',.12)]:ax.scatter(x+offset,d[key+'_'+name],s=8,marker=marker,c=color,label='KS' if name=='KS' else 'Normalized W1',linewidths=.3)
 ax.axhline(.03,c=AMBER,ls='--',lw=1.2);ax.set(xlabel='Pair index',ylabel='Marginal distance',ylim=(0,.032),xlim=(0,49));ax.set_xticks([1,16,32,48])
 ax.legend(loc='center',bbox_to_anchor=(.65,.65),ncol=1)
def s4():
 d=read('source_figS4').sort_values(['W','realization']);fig,axs=canvas()
 for ax,key in zip(axs.flat,['length','relative_angle','dy']):paired_diag(ax,d,key)
 ax=axs[1,1];x=np.arange(48)+1;ax.scatter(x,d.I0_relative_change,c=GRAY,s=9,marker='o',label='Rate change');ax.axhline(.03,c=AMBER,ls='--');ax.set(xlabel='Pair index',ylabel=r'$|\Delta I_0|/\max(I_{0,o},0.01)$',ylim=(0,.04),xlim=(0,49),xticks=[1,16,32,48])
 twin=ax.twinx();twin.scatter(x,d.changed_edge_fraction,c='#333333',s=9,marker='s');twin.axhline(.25,c=AMBER,ls=':');twin.set(ylabel='Changed-edge fraction',ylim=(0,1));ax.legend(loc='upper center');save(fig,4)
def s5():
 d=read('source_figS5');acc=d[d.accepted];fig,axs=canvas(h=113)
 for i,(ell,g) in enumerate(acc.groupby('ell_over_R')):
  x=i+(np.arange(len(g))-(len(g)-1)/2)/max(len(g),1)*.30
  axs[0,0].scatter(x,g.log10_S_ratio,s=8,c=MAG,alpha=.35,linewidths=0);err(axs[0,0],i,g.log10_S_ratio,MAG)
  axs[0,0].text(i,.98,f'{len(g)}/48',transform=axs[0,0].get_xaxis_transform(),ha='center',va='top')
 for j,(key,c,marker) in enumerate([('A_loc',G1,'o'),('D',G2,'s'),('f_up',GRAY,'^')]):
  for i,(ell,g) in enumerate(acc.groupby('ell_over_R')):err(axs[0,1],i,g[key+'_reassigned']/g[key+'_original'],c,label={'A_loc':r'$A_{\rm loc}$','D':r'$D$','f_up':r'$f_\uparrow$'}[key] if i==0 else None,marker=marker,offset=(j-1)*.13)
 axs[0,1].axhline(1,c=GRAY,ls='--',lw=.6);axs[0,1].set_ylim(.994,1.007);axs[0,1].legend(loc='upper center',ncol=3)
 for i,(ell,g) in enumerate(d.groupby('ell_over_R')):
  for kind,c,marker,off in [('KS',GRAY,'o',-.09),('Wasserstein','#333333','s',.09)]:
   worst=g[[x for x in g if x.endswith('_'+kind)]].max(axis=1);axs[1,0].scatter(np.full(len(g),i+off),worst,c=c,s=8,marker=marker,alpha=.35,label=kind if i==0 else None)
  accepted=g.accepted.to_numpy();axs[1,1].scatter(np.full(len(g),i),g.I0_relative_change,c=GRAY,s=9,marker='o',alpha=.45)
  reject=g[~g.accepted];axs[1,1].scatter(np.full(len(reject),i),reject.I0_relative_change,c='#111111',s=22,marker='x',label='Rejected' if len(reject) else None)
 axs[1,0].axhline(.03,c=AMBER,ls='--');axs[1,0].set(ylabel='Largest marginal distance',ylim=(0,.032));axs[1,0].legend(loc='center',bbox_to_anchor=(.65,.65))
 axs[1,1].axhline(.03,c=AMBER,ls='--');axs[1,1].set(ylabel='Relative rate change',ylim=(0,.05));axs[1,1].legend(loc='upper left')
 twin=axs[1,1].twinx()
 for i,(ell,g) in enumerate(d.groupby('ell_over_R')):twin.scatter(np.full(len(g),i+.13),g.changed_edge_fraction,s=7,c='#333333',marker='s',alpha=.3)
 twin.axhline(.25,c=AMBER,ls=':');twin.set(ylabel='Changed-edge fraction',ylim=(0,1))
 axs[0,0].set(ylabel=r'$\log_{10}(S_r/S_o)$');axs[0,1].set(ylabel='Paired local-statistic ratio')
 for ax in axs.flat:ax.set(xlabel=r'$\ell/R$',xticks=range(4),xticklabels=['1','2','4','8'],xlim=(-.5,3.5))
 save(fig,5)
def s6():
 nodes=read('classifier_toy_nodes');edges=read('classifier_toy_edges');v=read('source_figS6');fig,axs=canvas(h=95)
 for ax,name,col in [(axs[0,0],'B_pure_order1',G1),(axs[0,1],'C_reciprocal_order2',G2)]:
  n=nodes[nodes.toy.eq(name)].sort_values('node');ys=n.progress_y.to_numpy();xs=np.array([0,0]) if len(n)==2 else np.array([0,1,0,1]);pos={i:(xs[i],ys[i]) for i in range(len(n))}
  for r in edges[edges.toy.eq(name)].itertuples():
   u,w=int(r.source),int(r.target);color=col if ys[w]>ys[u] else GRAY
   ax.add_patch(FancyArrowPatch(pos[u],pos[w],arrowstyle='-|>',mutation_scale=10,color=color,linewidth=1.3,shrinkA=8,shrinkB=8,connectionstyle='arc3,rad=.14' if len(n)==2 else 'arc3,rad=0'))
  ax.scatter(xs,ys,s=95,c='white',edgecolors='#111111',zorder=3)
  for i,(x,y) in enumerate(zip(xs,ys)):ax.text(x,y,str(i),ha='center',va='center',fontsize=8,zorder=4)
  ax.set(xlim=(-.6,.6) if len(n)==2 else (-.35,1.35),ylim=(ys.min()-.5,ys.max()+.5));ax.set_axis_off();ax.text(.5,.02,'Order 1' if len(n)==2 else 'Order 2',transform=ax.transAxes,ha='center')
 axs[1,0].set_axis_off();axs[1,0].text(.03,.93,r'$E=D_0\cup F$'+'\n\n'+r'$F_{\leq1}\subseteq F_{\leq2}\subseteq F$'+'\n\n'+r'$G_1\subseteq G_2\subseteq G$',transform=axs[1,0].transAxes,fontsize=11,va='top')
 ax=axs[1,1];toy=v[v.parent_graph_id.isna()];hist=v[v.parent_graph_id.notna()]
 ax.texts[0].set_x(.02)
 ax.plot([0,1],[toy.classification_exact.sum(),hist.classification_exact.sum()],ls='',marker='o',c=GRAY,ms=5)
 ax.set(xticks=[0,1],xticklabels=['Archived toys','Induced graphs'],ylabel='Independent exact matches',ylim=(0,46),xlim=(-.4,1.4));ax.text(0,8,'5/5',ha='center');ax.text(1,42,'40/40',ha='center');save(fig,6)
def s7():
 d=read('source_figS7');d=d[d.W.eq(.805)];s=read('source_figS7_scaling');s=s[s.control.eq(.805)];direction=read('source_figS7_direction');fig,axs=canvas(h=113)
 for j,(L,g) in enumerate(d.groupby('L')):
  for k,c,marker,off in [(1,G1,'o',-.1),(2,G2,'s',.1)]:
   a=g[f'G{k}_d90'];x=j+off+(np.arange(len(a))-(len(a)-1)/2)/len(a)*.10
   axs[0,0].scatter(x,a,c=c,s=7,alpha=.3,marker=marker,linewidths=0);err(axs[0,0],j,a,c,label=rf'$G_{k}$' if j==0 else None,marker=marker,offset=off)
   ecdf(axs[1,1],a,c,rf'$G_{k},\ {L:g}$',STYLES[j])
 for k,c in [(1,G1),(2,G2)]:
  for key,ls,marker in [('ell_major_over_L_','-','o'),('ell_minor_over_L_','--','s')]:summary(axs[0,1],s,'L',key+f'G{k}',c,label=rf'$G_{k}$ '+('major' if ls=='-' else 'minor'),marker=marker,ls=ls)
 for j,name in enumerate(['ORIGINAL','REVERSED','ROTATED_90']):
  g=direction[direction.orientation_transform.eq(name)].sort_values('graph_pair_id')
  for k,c,marker,off in [(1,G1,'o',-.1),(2,G2,'s',.1)]:
   a=g[f'DLOC_G{k}'];axs[1,0].scatter(j+off+(np.arange(len(a))-11.5)/24*.1,a,c=c,s=7,alpha=.3,marker=marker);err(axs[1,0],j,a,c,label=rf'$G_{k}$' if j==0 else None,marker=marker,offset=off)
 axs[0,0].set(xticks=range(4),xticklabels=['64','128','256','512'],xlabel=r'$L/R$',ylabel=r'$d_{90}/R$',yscale='log');axs[0,0].legend(ncol=2)
 axs[0,1].set(xlabel=r'$L/R$',ylabel='Span / L',yscale='log');size_axis(axs[0,1]);axs[0,1].legend(ncol=2,loc='center',bbox_to_anchor=(.5,.80))
 axs[1,0].set(xticks=range(3),xticklabels=['Original','Reversed','Rotated 90°'],ylabel=r'$\mathrm{DLOC}_k$',ylim=(-.04,1.05));axs[1,0].legend(ncol=2,loc='lower left')
 axs[1,1].set(xlabel=r'$d_{90}/R$',ylabel='ECDF',xscale='log',ylim=(0,1.02));axs[1,1].legend(ncol=2,loc='lower right',fontsize=7.5)
 save(fig,7)
def s8():
 d=read('source_figS8');fig,axs=canvas(h=113)
 for ax,prefix,label in zip(axs.flat,['f_internal_rec_','eta_internal_','internal_nontrivial_scc_count_','S_internal_rec_core_'],[r'$f^{\rm int}_{k,rec}$',r'$\eta^{\rm int}_k$','Nontrivial internal SCC count',r'$S_k^{\rm int}$']):
  for k,c,m in [(1,G1,'o'),(2,G2,'s')]:
   assert (d[prefix+str(k)]>0).all(),('Zeros require explicit log handling',prefix,k)
   summary(ax,d,'W',prefix+str(k),c,label=rf'$G_{k}$',marker=m)
  ax.set(xlabel=r'$W$',ylabel=label);ax.legend(ncol=2)
 axs[0,1].set_yscale('log');axs[1,0].set_yscale('log');axs[1,1].set_yscale('log');save(fig,8)
def s9():
 d=read('source_figS9');fig,axs=canvas(h=113);core=d[d.width_over_R.eq(4)]
 for j,(w,g) in enumerate(core.groupby('W')):
  x=j+(np.arange(len(g))-(len(g)-1)/2)/len(g)*.22
  for xx,a,b in zip(x,g.S_full_core,g.S2_core):axs[0,0].plot([xx-.1,xx+.1],[a,b],c='#DDDDDD',lw=.4,zorder=0)
  for y,c,m,off in [('S_full_core',MAG,'o',-.1),('S2_core',G2,'s',.1)]:
   axs[0,0].scatter(x+off,g[y],s=7,c=c,alpha=.4,marker=m,linewidths=0);err(axs[0,0],j,g[y],c,label=r'$G[I]$' if y=='S_full_core' and j==0 else (r'$G_2[I]$' if j==0 else None),marker=m,offset=off)
  axs[0,0].text(j,1.05,f'n={len(g)}',transform=axs[0,0].get_xaxis_transform(),ha='center',va='bottom',fontsize=7.5)
  axs[0,1].scatter(x,g.delta_full_G2,s=7,c=GRAY,alpha=.4);err(axs[0,1],j,g.delta_full_G2,GRAY)
 for j,(w,g) in enumerate(d.groupby('W')):
  for y,c,m,ls in [('S_full_core',MAG,MARKERS[j],STYLES[j]),('S2_core',G2,MARKERS[j],STYLES[j])]:summary(axs[1,0],g,'width_over_R',y,c,label=rf'$W={w:g}$' if y=='S2_core' else None,marker=m,ls=ls)
  summary(axs[1,1],g,'width_over_R','retained_fraction',GRAY,label=rf'$W={w:g}$',marker=MARKERS[j],ls=STYLES[j])
 for ax in axs[0]:ax.set(xlabel=r'$W$',xticks=range(4),xticklabels=['.805','.815','.86','.90'],xlim=(-.5,3.5))
 axs[0,0].set(ylabel='Retained-core SCC fraction',ylim=(-.02,1.03));axs[0,0].legend(loc='center right')
 axs[0,1].set(ylabel=r'$S(G[I])-S(G_2[I])$',ylim=(0,1.03))
 axs[1,0].set(ylabel='Retained-core SCC fraction',ylim=(-.02,1.03));axs[1,0].legend(ncol=2,loc='center',bbox_to_anchor=(.52,.31))
 axs[1,1].set(ylabel='Retained-node fraction',ylim=(0,1.03));axs[1,1].legend(ncol=2,loc='lower left')
 for ax in axs[1]:ax.set(xlabel=r'$\Delta/R$',xticks=[2,4,8])
 save(fig,9)
for f in [s1,s2,s3,s4,s5,s6,s7,s8,s9]:f()
