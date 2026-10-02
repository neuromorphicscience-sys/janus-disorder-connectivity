"""Supplemental figures from complete source tables, using the manuscript palette."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
PROJECT=Path(__file__).resolve().parents[1]
ROOT=PROJECT/'data/supplemental'
OUT=PROJECT/'figures/generated/supplemental';OUT.mkdir(parents=True,exist_ok=True)
COL={'full':'#D81B60','g1':'#004D40','g2':'#1E88E5','ref':'#FFC107','gray':'#777777'}
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':9,'axes.labelsize':10,
    'legend.fontsize':8,'axes.linewidth':.7,'xtick.direction':'in','ytick.direction':'in',
    'xtick.top':True,'ytick.right':True,'pdf.fonttype':42,'svg.fonttype':'none','lines.linewidth':1.1,'lines.markersize':3})

def canvas(rows=2,cols=2,height=105):
    fig,axes=plt.subplots(rows,cols,figsize=(165/25.4,height/25.4),layout='constrained',squeeze=False)
    for i,ax in enumerate(axes.flat):
        ax.text(-.12,1.04,f'({chr(97+i)})',transform=ax.transAxes,weight='bold')
        ax.tick_params(which='both',top=True,right=True)
    return fig,axes
def save(fig,n):
    fig.savefig(OUT/f'figS{n}.pdf');fig.savefig(OUT/f'figS{n}.svg');fig.savefig(OUT/f'figS{n}.png',dpi=250)
    plt.close(fig)
def curve(ax,d,x='control',c=COL['full'],label=None,marker='o',ls='-'):
    d=d.sort_values(x);a=d[x].to_numpy();m=d.S_median.to_numpy()
    ax.plot(a,m,marker=marker,color=c,label=label,ls=ls)
    ax.fill_between(a,d.S_Q1.to_numpy(),d.S_Q3.to_numpy(),color=c,alpha=.12,lw=0)
def medcurve(ax,x,y,color,label,marker='o',ls='-'):
    ax.plot(x,y,marker=marker,color=color,label=label,ls=ls)
def groupsummary(d,x,y):
    g=d.groupby(x)[y];return g.median(),g.quantile(.25),g.quantile(.75)

def robust():
    s=pd.read_csv(ROOT/'04_parameter_robustness/ROBUSTNESS_SUMMARY.csv')
    q=s[s.distribution.eq('gaussian')&s.correlation_family.eq('iid')&s.L.eq(128)&s.sigma.eq(24)]
    fig,axs=canvas(2,4,100)
    for ax,k in zip(axs.flat,[4,5,6,7,8,10,12]):
        d=q[q.q.eq(k)];curve(ax,d);ax.set(xlabel=r'$W$',ylabel=r'$S$',ylim=(-.02,1.02))
        ax.text(.06,.90,rf'$q={k}$',transform=ax.transAxes)
    axs.flat[7].texts[0].remove()
    axs.flat[7].axis('off')
    axs.flat[7].text(0,.83,r'$L/R=128$'+'\n'+r'$\sigma R^2=24$'+'\n\nMedians and IQRs\nAll archived conditions',transform=axs.flat[7].transAxes,fontsize=10)
    save(fig,1)
    fig,axs=canvas(3,2,155)
    cases=[('gaussian','iid','W','Gaussian'),('uniform','iid','A','Bounded uniform'),('vonmises','iid',r'\kappa','von Mises'),('laplace','iid','b','Wrapped Laplace')]
    for ax,(dist,fam,xlab,name) in zip(axs.flat,cases):
        d=s[s.distribution.eq(dist)&s.correlation_family.eq(fam)&s.q.eq(7)&s.L.eq(128)&s.sigma.eq(24)]
        curve(ax,d);ax.set(xlabel=rf'${xlab}$',ylabel=r'$S$',ylim=(-.02,1.02));ax.text(.06,.90,name,transform=ax.transAxes)
    for ax,anisotropic in zip(axs.flat[4:],[False,True]):
        sel=s[s.q.eq(7)&s.L.eq(128)&s.sigma.eq(24)&s.control_name.eq('alpha')]
        sel=sel[sel.correlation_family.str.startswith('ANISO') if anisotropic else ~sel.correlation_family.str.startswith('ANISO')]
        for j,((fam,ex,ey),d) in enumerate(sel.groupby(['correlation_family','ell_x','ell_y'])):
            curve(ax,d,c=COL['full'],label=rf'$({ex:g},{ey:g})R$',marker=['o','s','^','D'][j%4],ls=['-','--','-.',':'][j%4])
        ax.set(xlabel=r'$\alpha$',ylabel=r'$S$',ylim=(-.02,1.02));ax.legend(title=r'$(\ell_x,\ell_y)$',loc='best',frameon=False)
    save(fig,2)

def intervention():
    d=pd.read_csv(ROOT/'02_ell_sensitivity/GRAPHLEVEL.csv');accepted=d[d.accepted]
    fig,axs=canvas(2,2,112)
    ax=axs[0,0]
    for ell,g in accepted.groupby('ell_over_R'):
        jitter=(g.realization.to_numpy()-7.5)/75
        ax.scatter(ell+jitter,g.log10_S_ratio,color=COL['full'],s=10,alpha=.4)
        m=g.log10_S_ratio.median();q=g.log10_S_ratio.quantile([.25,.75]).to_numpy()
        ax.errorbar(ell,m,yerr=[[m-q[0]],[q[1]-m]],color=COL['full'],capsize=3,marker='o',ms=5)
        ax.text(ell,.98,f'{len(g)}/48',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=8)
    ax.axhline(0,color=COL['gray'],lw=.6);ax.set(xlabel=r'$\ell/R$',ylabel=r'$\log_{10}(S_r/S_o)$')
    ax=axs[0,1]
    for key,color,label in [('A_loc',COL['g1'],r'$A_{\rm loc}$'),('D',COL['g2'],r'$D$'),('f_up',COL['gray'],r'$f_{\uparrow}$')]:
        for ell,g in accepted.groupby('ell_over_R'):
            ratio=g[key+'_reassigned']/g[key+'_original'];m=ratio.median();q=ratio.quantile([.25,.75]).to_numpy()
            ax.errorbar(ell,m,yerr=[[m-q[0]],[q[1]-m]],color=color,marker='o',capsize=2,label=label if ell==1 else None)
    ax.axhline(1,color=COL['gray'],lw=.6,ls='--');ax.set(xlabel=r'$\ell/R$',ylabel='Paired local-statistic ratio');ax.legend(frameon=False)
    ax=axs[1,0]
    for kind,col in [('KS',COL['g2']),('Wasserstein',COL['g1'])]:
        for ell,g in d.groupby('ell_over_R'):
            worst=g[[c for c in d if c.endswith('_'+kind)]].max(axis=1)
            ax.scatter(np.full(len(g),ell+(0.12 if kind=='KS' else -.12)),worst,s=8,c=col,alpha=.35,label=kind if ell==1 else None)
    ax.axhline(.03,color=COL['gray'],ls='--',label='Common limit');ax.set(xlabel=r'$\ell/R$',ylabel='Largest marginal distance');ax.legend(frameon=False)
    ax=axs[1,1]
    for ell,g in d.groupby('ell_over_R'):
        ax.scatter(np.full(len(g),ell-.12),g.changed_edge_fraction,s=8,color=COL['gray'],alpha=.35,label='Changed edges' if ell==1 else None)
        ax.scatter(np.full(len(g),ell+.12),g.I0_relative_change,s=8,color=COL['g2'],alpha=.35,label=r'$|\Delta I_0|/I_0$' if ell==1 else None)
    ax.axhline(.25,c=COL['gray'],ls='--',lw=.7);ax.axhline(.03,c=COL['g2'],ls='--',lw=.7)
    ax.set(xlabel=r'$\ell/R$',ylabel='Fraction / relative change');ax.legend(frameon=False)
    for ax in axs.flat:ax.set_xticks([1,2,4,8])
    save(fig,4)
    d=d[d.ell_over_R.eq(4)].sort_values(['W','realization']);fig,axs=canvas(2,2,110)
    ax=axs[0,0]
    for r in d.to_dict('records'):ax.plot([0,1],[r['S_original'],r['S_reassigned']],color='#D8D8D8',lw=.5,zorder=0)
    ax.scatter(np.zeros(len(d)),d.S_original,c=COL['full'],s=13);ax.scatter(np.ones(len(d)),d.S_reassigned,c=COL['gray'],s=13)
    ax.set(yscale='log',ylabel=r'$S$',xticks=[0,1],xticklabels=['Original','Reassigned'])
    ax=axs[0,1]
    for j,label in enumerate(['length','relative_angle','dy']):
        ax.scatter(j-.10+np.arange(len(d))/500,d[label+'_KS'],s=8,c=COL['g2'],alpha=.35,label='KS' if j==0 else None)
        ax.scatter(j+.10+np.arange(len(d))/500,d[label+'_Wasserstein'],s=8,c=COL['g1'],alpha=.35,label='Wasserstein' if j==0 else None)
    ax.axhline(.03,c=COL['gray'],ls='--');ax.set(ylabel='Marginal distance',xticks=[0,1,2],xticklabels=['Length','Rel. angle',r'$\Delta y$']);ax.legend(frameon=False)
    ax=axs[1,0]
    for q in [.1,.25,.5,.75,.9,.95]:
        a=d[f'length_q{q:g}_change'];ax.scatter(np.full(len(a),q),a,s=6,c=COL['gray'],alpha=.25)
        m=a.median();lo,hi=a.quantile([.25,.75]);ax.errorbar(q,m,yerr=[[m-lo],[hi-m]],c=COL['gray'],marker='o',capsize=2)
    ax.axhline(0,c=COL['gray'],lw=.7);ax.set(xlabel='Edge-length quantile probability',ylabel=r'$\Delta\ell_u/R$')
    ax=axs[1,1]
    ax.scatter(d.I0_relative_change,d.log10_S_ratio,c=COL['full'],s=16);ax.set(xlabel=r'$|\Delta I_0|/I_0$',ylabel=r'$\log_{10}(S_r/S_o)$')
    save(fig,3)

def size_and_bulk():
    s=pd.read_csv(ROOT/'05_finite_size_diagnostics/FINITE_SIZE_SUMMARY.csv');c=pd.read_csv(ROOT/'05_finite_size_diagnostics/CROSSOVER_DIAGNOSTICS.csv')
    fig,axs=canvas(2,2,110)
    for j,(L,g) in enumerate(s.groupby('L')):
        g=g.sort_values('W');color=COL['full'];marker=['o','s','^','D'][j];ls=['-','--','-.',':'][j]
        curve(axs[0,0],g,x='W',c=color,label=rf'$L/R={L:g}$',marker=marker,ls=ls)
        medcurve(axs[0,1],g.W,g.S_variance,color,rf'${L:g}$',marker,ls)
    axs[0,0].set(xlabel=r'$W$',ylabel=r'$S$',xlim=(.75,.84),ylim=(-.02,1.02));axs[0,0].legend(frameon=False)
    axs[0,1].set(xlabel=r'$W$',ylabel=r'${\rm Var}(S)$',xlim=(.75,.84))
    for level,color,marker in [(.1,COL['g1'],'v'),(.5,COL['full'],'o'),(.9,COL['g2'],'^')]:
        medcurve(axs[1,0],c.L,c[f'W_{level:g}'],color,rf'$W_{{{level:g}}}$',marker)
    axs[1,0].set(xlabel=r'$L/R$',ylabel='Interpolated crossing');axs[1,0].legend(frameon=False)
    medcurve(axs[1,1],c.L,c.Delta_W_10_90,COL['gray'],None)
    axs[1,1].set(xlabel=r'$L/R$',ylabel=r'$\Delta W_{10-90}$')
    for ax in axs[1]:ax.set_xscale('log',base=2);ax.set_xticks([64,128,256,512],labels=['64','128','256','512'])
    save(fig,5)
    d=pd.read_csv(ROOT/'05_finite_size_diagnostics/SCALING_784_GRAPHLEVEL.csv');d=d[d.control.isin([.7975,.805,.815])]
    fig,axs=canvas(2,2,110)
    for k,color in [(1,COL['g1']),(2,COL['g2'])]:
        ax=axs[0,k-1]
        for j,(w,g) in enumerate(d.groupby('control')):
            med,lo,hi=groupsummary(g,'L',f'L_times_S{k}')
            ax.errorbar(med.index,med,yerr=[med-lo,hi-med],c=color,marker=['o','s','^'][j],ls=['-','--',':'][j],label=rf'$W={w:g}$',capsize=2)
        ax.set(xlabel=r'$L/R$',ylabel=rf'$LS_{k}/R$');ax.legend(frameon=False)
    g=d[d.control.eq(.805)]
    for name,col,label,ls in [('ell_major_over_L_G2',COL['g2'],'Major span','-'),('ell_minor_over_L_G2',COL['g2'],'Minor span','--')]:
        med,lo,hi=groupsummary(g,'L',name);axs[1,0].errorbar(med.index,med,yerr=[med-lo,hi-med],c=col,label=label,marker='o',ls=ls,capsize=2)
    axs[1,0].set(xlabel=r'$L/R$',ylabel='Span / L',yscale='log');axs[1,0].legend(frameon=False)
    for k,color in [(1,COL['g1']),(2,COL['g2'])]:
        med,lo,hi=groupsummary(g,'L',f'bottom_fraction_2R_G{k}');axs[1,1].errorbar(med.index,med,yerr=[med-lo,hi-med],c=color,label=rf'$G_{k}$',marker='o',capsize=2)
    axs[1,1].set(xlabel=r'$L/R$',ylabel=r'Downstream $2R$ strip fraction',ylim=(0,1.05));axs[1,1].legend(frameon=False)
    for ax in axs.flat:ax.set_xscale('log',base=2);ax.set_xticks([64,128,256,512],labels=['64','128','256','512'])
    save(fig,6)
    d=pd.read_csv(ROOT/'06_si/INTERNAL_RECURRENCE_1528.csv');fig,axs=canvas(2,2,110)
    panels=[('f_internal_rec_',r'$f^{\rm int}_{k,rec}$'),('eta_internal_',r'$\eta^{\rm int}_k$'),('internal_nontrivial_scc_count_','Nontrivial internal SCC count')]
    for ax,(prefix,label) in zip(axs.flat,panels):
        for k,col in [(1,COL['g1']),(2,COL['g2'])]:
            med,lo,hi=groupsummary(d,'W',prefix+str(k));ax.plot(med.index,med,c=col,label=rf'$G_{k}$',marker='o');ax.fill_between(med.index,lo,hi,color=col,alpha=.12,lw=0)
        ax.axvspan(.7975,.815,color='#E3E3E3',alpha=.45,zorder=0);ax.set(xlabel=r'$W$',ylabel=label);ax.legend(frameon=False)
    axs[0,1].set_yscale('log');axs[1,0].set_yscale('log')
    for k,col in [(1,COL['g1']),(2,COL['g2'])]:
        d['retention']=d[f'internal_recurrent_vertices_{k}']/d[f'parent_recurrent_vertices_{k}'];med,lo,hi=groupsummary(d,'W','retention');axs[1,1].plot(med.index,med,color=col,marker='o',label=rf'$G_{k}$');axs[1,1].fill_between(med.index,lo,hi,color=col,alpha=.12,lw=0)
    axs[1,1].set(xlabel=r'$W$',ylabel='Internal / parent recurrent mass',ylim=(.85,1.01));axs[1,1].legend(frameon=False)
    save(fig,7)

def deletion():
    d=pd.read_csv(ROOT/'03_boundary_deletion_ensemble/GRAPHLEVEL.csv');r=pd.read_csv(ROOT/'03_boundary_deletion_ensemble/REPRESENTATIVE_PLACEMENT.csv',dtype={'geometry_seed':str,'disorder_seed':str,'classification_sha256':str,'edge_hash':str},engine='python')
    fig,axs=canvas(2,2,118)
    ax=axs[0,0];width4=d[d.width_over_R.eq(4)]
    for j,(w,g) in enumerate(width4.groupby('W')):
        x=j+(np.arange(len(g))-len(g)/2)/max(len(g),1)*.28
        for xx,a,b in zip(x,g.S_full_core,g.S2_core):ax.plot([xx-.1,xx+.1],[a,b],c='#D8D8D8',lw=.5)
        ax.scatter(x-.1,g.S_full_core,c=COL['full'],s=8,alpha=.5,label=r'$G[I]$' if j==0 else None)
        ax.scatter(x+.1,g.S2_core,c=COL['g2'],s=8,alpha=.5,label=r'$G_2[I]$' if j==0 else None)
        ax.text(j,.98,f'n={len(g)}',ha='center',va='top',transform=ax.get_xaxis_transform(),fontsize=8)
    ax.set(ylabel='Retained-core SCC fraction',xticks=range(4),xticklabels=['.805','.815','.86','.90'],xlabel=r'$W$',ylim=(-.02,1.05));ax.legend(loc='center right',frameon=False)
    ax=axs[0,1]
    for j,(w,g) in enumerate(d.groupby('W')):
        for k,col,ls in [('S_full_core',COL['full'],'-'),('S2_core',COL['g2'],'--')]:
            med,lo,hi=groupsummary(g,'width_over_R',k);ax.errorbar(med.index,med,yerr=[med-lo,hi-med],c=col,marker=['o','s','^','D'][j],ls=ls,capsize=2,label=rf'$W={w:g}$' if k=='S2_core' else None)
    ax.set(xlabel=r'Deletion width $w/R$',ylabel='Core fraction',xticks=[2,4,8]);ax.legend(frameon=False,ncol=2,loc='center',bbox_to_anchor=(.55,.34))
    ax=axs[1,0]
    for j,(w,g) in enumerate(width4.groupby('W')):
        x=np.sort(g.delta_full_G2);ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',c=[COL['full'],COL['g1'],COL['g2'],COL['gray']][j],label=rf'$W={w:g}$')
    ax.set(xlabel=r'$S(G[I])-S(G_2[I])$',ylabel='ECDF');ax.legend(frameon=False)
    ax=axs[1,1]
    for key,col,label,marker in [('S_full_core',COL['full'],r'$G[I]$','o'),('S2_core',COL['g2'],r'$G_2[I]$','s')]:
        ax.plot(range(4),r.sort_values('W')[key+'_empirical_percentile'],marker=marker,ls='',c=col,label=label)
    ax.axhline(.5,c=COL['gray'],ls='--',lw=.7);ax.set(ylabel='Representative percentile',xlabel=r'$W$',xticks=range(4),xticklabels=['.805','.815','.86','.90'],ylim=(-.02,1.02));ax.legend(frameon=False)
    save(fig,8)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--existing-only',action='store_true');a=ap.parse_args()
    robust();size_and_bulk()
    if not a.existing_only:intervention();deletion()
