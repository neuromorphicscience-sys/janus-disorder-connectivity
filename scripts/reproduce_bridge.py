"""Render Fig. S10 from public paired source data, with undefined eta explicit."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=Path('data/final_bridge'));ap.add_argument('--output',type=Path,default=Path('figures/generated/supplemental'));a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8.5,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
    C={'S':'#D81B60','1':'#004D40','2':'#1E88E5'};d=pd.read_csv(a.data/'INTERVENTION_BRIDGE_SOURCE.csv');corr=pd.read_csv(a.data/'ORIENTATION_CORRELATION_SOURCE.csv')
    assert len(d)==48
    fig,axes=plt.subplots(2,2,figsize=(178/25.4,115/25.4));fig.subplots_adjust(left=.10,right=.98,bottom=.12,top=.95,wspace=.35,hspace=.48)
    for ax,label in zip(axes.flat,'abcd'):ax.text(-.18,1.03,label,transform=ax.transAxes,fontweight='bold',fontsize=10)
    ax=axes[0,0]
    for i,(name,color,marker) in enumerate([('S',C['S'],'o'),('S1',C['1'],'s'),('S2',C['2'],'^')]):
        v=d[f'log10_{name}_ratio'].to_numpy();jitter=(np.arange(len(v))%16-7.5)/60
        ax.scatter(i+jitter,v,s=9,marker=marker,facecolors='none',edgecolors=color,alpha=.5,lw=.5)
        q=np.quantile(v,[.25,.5,.75]);ax.errorbar(i,q[1],yerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt=marker,color=color,ms=5,capsize=3,lw=1.3)
    ax.axhline(0,color='.7',lw=.6);ax.set(xticks=[0,1,2],xticklabels=['$S$','$S_1$','$S_2$'],ylabel=r'$\log_{10}(S_{k,r}/S_{k,o})$',ylim=(-4.2,.2),xlim=(-.4,2.4));ax.text(.97,.96,'48 pairs; 16 shared streams',transform=ax.transAxes,ha='right',va='top',fontsize=7.5)
    ax=axes[0,1]
    for k,marker in [(1,'s'),(2,'^')]:
        i=2*(k-1);vo=d[f'f_int_rec_{k}_original'].to_numpy();vr=d[f'f_int_rec_{k}_reassigned'].to_numpy();jitter=(np.arange(len(vo))%16-7.5)/65
        for j in range(len(vo)):ax.plot([i+jitter[j],i+1+jitter[j]],[vo[j],vr[j]],color=C[str(k)],alpha=.16,lw=.55)
        for pos,v in [(i,vo),(i+1,vr)]:
            ax.scatter(pos+jitter,v,s=9,marker=marker,facecolors='none',edgecolors=C[str(k)],lw=.5,alpha=.5)
            q=np.quantile(v,[.25,.5,.75]);ax.errorbar(pos,q[1],yerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt=marker,color=C[str(k)],ms=5,capsize=3,lw=1.3)
    ax.set(xticks=range(4),xticklabels=['$1,o$','$1,r$','$2,o$','$2,r$'],ylabel=r'$f^{\rm int}_{k,\rm rec}$',ylim=(-.008,.15),xlim=(-.4,3.4));ax.text(.03,.96,'strict 4R core\nn=48/order',transform=ax.transAxes,ha='left',va='top',fontsize=7.5)
    ax=axes[1,0]
    for k,marker in [(1,'s'),(2,'^')]:
        pos=2*(k-1);v=d[f'eta_int_{k}_original'].to_numpy();assert d[f'eta_int_{k}_reassigned'].isna().all()
        jitter=(np.arange(len(v))%16-7.5)/60;ax.scatter(pos+jitter,v,s=9,marker=marker,facecolors='none',edgecolors=C[str(k)],alpha=.5,lw=.5)
        q=np.quantile(v,[.25,.5,.75]);ax.errorbar(pos,q[1],yerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt=marker,color=C[str(k)],ms=5,capsize=3,lw=1.3)
        ax.text(pos+1,.005,'undefined',rotation=90,ha='center',va='center',fontsize=7.5,color='.4')
    ax.set(yscale='log',xticks=range(4),xticklabels=['$1,o$','$1,r$','$2,o$','$2,r$'],ylabel=r'$\eta^{\rm int}_k$',ylim=(.0005,.12),xlim=(-.4,3.4));ax.text(.03,.96,'reassigned: no\nrecurrent core',transform=ax.transAxes,ha='left',va='top',fontsize=7.5)
    ax=axes[1,1]
    for state,col,marker,ls in [('original','.4','o','--'),('reassigned',C['S'],'s','-')]:
        x=[];med=[];lo=[];hi=[]
        for (rmin,rmax),g in corr[corr.state==state].groupby(['r_min_R','r_max_R']):
            q=g.C_theta_conn.quantile([.25,.5,.75]);x.append((rmin+rmax)/2);lo.append(q.iloc[0]);med.append(q.iloc[1]);hi.append(q.iloc[2])
        ax.plot(x,med,color=col,marker=marker,ls=ls,ms=4,lw=1.2,label=state);ax.fill_between(x,lo,hi,color=col,alpha=.12)
    ax.axhline(0,color='.7',lw=.6);ax.set(xlabel=r'radial bin center $r/R$',ylabel=r'$C_\theta^{\rm conn}(r)$',xlim=(0,6.5),ylim=(-.04,.53),xticks=[.5,1.5,3,6]);ax.legend(loc='upper right',fontsize=7.5)
    for extension in ['pdf','svg','png']:fig.savefig(a.output/f'figS10.{extension}',dpi=300,facecolor='white')
    plt.close(fig)
if __name__=='__main__':main()
