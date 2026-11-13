#!/usr/bin/env python3
"""Re-analysis of processed GEO data (GSE178742 logCPM, GSE168740 FPKM, GSE121885 counts).
Mirrors the in-browser analysis of 2026-11-11 (scripts/differential_expression/browser_analysis.js) so results
can be regenerated from files on disk after scripts/download/download_geo.py.

Methods (processed data only - NOT a unified raw-read pipeline):
 * GSE178742 apiary: colony = experimental unit. Mean logCPM per colony (5 bees), Student t (3 vs 3 colonies, df=4).
   Bee-level Welch t reported only to show the effect of pseudoreplication. ICC(colony) from one-way ANOVA within arm.
 * GSE178742 field: one colony-pooled library per colony; Student t (5v5 foragers, 5v4 nurses). Site confounded.
 * GSE168740: duplicate Gene_ID rows checked (all identical -> de-duplicated); genes with mean FPKM>=1 tested;
   Student t on log2(FPKM+1), 3 v 3 pooled libraries per age.
 * GSE121885: median-of-ratios size factors; genes mean normalised count>=10; OLS on log2(norm+1) with effect-coded
   diet, virus, diet:virus (balanced 2x2, n=6/cell, df=20).
 * BH-FDR across all tested genes in each contrast.
"""
import os, sys, tarfile, gzip, io, argparse
import numpy as np, pandas as pd
from scipy import stats
ap=argparse.ArgumentParser(); ap.add_argument('--project-root',default='.'); a=ap.parse_args(); R=a.project_root
raw=lambda acc,f: os.path.join(R,'data','raw',acc,f)
out=os.path.join(R,'results','differential_expression'); os.makedirs(out,exist_ok=True)
def bh(p):
    p=np.asarray(p,float); q=np.full(p.shape,np.nan); m=~np.isnan(p); pv=p[m]; n=len(pv)
    o=np.argsort(pv); r=pv[o]*n/np.arange(1,n+1); r=np.minimum.accumulate(r[::-1])[::-1]; t=np.empty(n); t[o]=np.minimum(r,1); q[m]=t; return q
meta=pd.read_csv(os.path.join(R,'data','metadata','GSE178742_sample_metadata_verified.tsv'),sep='\t')
summ=[]
for exp in ['Apiary_Foragers','Apiary_Nurses','Field_Foragers','Field_Nurses']:
    f=raw('GSE178742',f'GSE178742_{exp}_Logcpm.txt.gz')
    if not os.path.exists(f): print('missing',f); continue
    X=pd.read_csv(f,sep='\t',index_col=0); X.columns=[c.strip() for c in X.columns]
    md=meta[meta.experiment==exp].set_index('matrix_column').loc[X.columns]
    ex=md.treatment.eq('Exposed').values
    if exp.startswith('Apiary'):
        cm=X.T.groupby(md.colony.values).mean().T; ctrt=md.groupby('colony').treatment.first()
        e=cm.loc[:,ctrt[cm.columns].eq('Exposed').values]; u=cm.loc[:,ctrt[cm.columns].eq('Unexposed').values]
        t,p=stats.ttest_ind(e,u,axis=1); lfc=e.mean(1)-u.mean(1)
        se=np.sqrt(((e.shape[1]-1)*e.var(1,ddof=1)+(u.shape[1]-1)*u.var(1,ddof=1))/(e.shape[1]+u.shape[1]-2)*(1/e.shape[1]+1/u.shape[1]))
        tn,pn=stats.ttest_ind(X.loc[:,ex],X.loc[:,~ex],axis=1,equal_var=False)
        res=pd.DataFrame({'logFC':lfc,'se':se,'t':t,'p':p,'naive_bee_level_p':pn},index=X.index)
        res['naive_fdr']=bh(res.naive_bee_level_p)
    else:
        e=X.loc[:,ex]; u=X.loc[:,~ex]; t,p=stats.ttest_ind(e,u,axis=1); lfc=e.mean(1)-u.mean(1)
        res=pd.DataFrame({'logFC':lfc,'t':t,'p':p},index=X.index)
    res['aveExpr']=X.mean(1); res['fdr']=bh(res.p); res.to_csv(os.path.join(out,f'GSE178742_{exp}_genomewide.tsv.gz'),sep='\t')
    summ.append((exp,len(res),int((res.fdr<0.05).sum()),int((res.p<0.01).sum())))
tf=raw('GSE168740','GSE168740_RAW.tar')
if os.path.exists(tf):
    cols={}
    with tarfile.open(tf) as t:
        for m in t.getmembers():
            d=pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(m).read())),sep='\t',dtype={'Gene_ID':str})
            s=d.columns[2]; g=d.groupby('Gene_ID')[s]
            assert (g.nunique()<=1).all(), f'non-identical duplicate rows in {m.name}'   # never silently drop differing rows
            cols[s]=g.first()
    F=pd.DataFrame(cols)
    for age in ['1d','5d','10d']:
        C=[c for c in F if c.startswith('C') and c.endswith('_'+age)]; T=[c for c in F if c.startswith('T') and c.endswith('_'+age)]
        k=F[C+T].mean(1)>=1; L=np.log2(F.loc[k,C+T]+1); t,p=stats.ttest_ind(L[T],L[C],axis=1)
        res=pd.DataFrame({'log2FC':L[T].mean(1)-L[C].mean(1),'t':t,'p':p,'meanFPKM_C':F.loc[k,C].mean(1),'meanFPKM_T':F.loc[k,T].mean(1)})
        res['fdr']=bh(res.p); res.to_csv(os.path.join(out,f'GSE168740_{age}_genomewide.tsv.gz'),sep='\t'); summ.append(('GSE168740_'+age,len(res),int((res.fdr<0.05).sum()),int((res.p<0.01).sum())))
cf=raw('GSE121885','GSE121885_processedData.txt.gz')
if os.path.exists(cf):
    Cn=pd.read_csv(cf,sep='\t',index_col=0); Cn=Cn[~Cn.index.str.startswith('__')]
    pos=Cn[(Cn>0).all(1)]; lg=np.log(pos); sf=np.exp((lg.sub(lg.mean(1),axis=0)).median(0))
    N=Cn/sf; N=N[N.mean(1)>=10]; Y=np.log2(N+1); n=Y.shape[1]
    diet=np.where(Y.columns.str[1]=='R',1,-1); vir=np.where(Y.columns.str[0]=='V',1,-1)
    D=np.column_stack([np.ones(n),diet,vir,diet*vir]); B=Y.values@D/n; rss=((Y.values-B@D.T)**2).sum(1); se=np.sqrt(rss/(n-4)/n)
    res=pd.DataFrame(index=Y.index); res['meanNorm']=N.mean(1)
    for name,est,s in [('diet',2*B[:,1],2*se),('virus',2*B[:,2],2*se),('interaction',4*B[:,3],4*se),('dietNoVirus',2*(B[:,1]-B[:,3]),2*np.sqrt(2)*se)]:
        res[name+'_log2FC']=est; res[name+'_p']=2*stats.t.sf(np.abs(est/s),n-4); res[name+'_fdr']=bh(res[name+'_p'])
    res.to_csv(os.path.join(out,'GSE121885_factorial_genomewide.tsv.gz'),sep='\t')
    for c in ['diet','virus','interaction','dietNoVirus']: summ.append(('GSE121885_'+c,len(res),int((res[c+'_fdr']<0.05).sum()),int((res[c+'_p']<0.01).sum())))
print(pd.DataFrame(summ,columns=['contrast','tested','FDR<0.05','p<0.01']).to_string(index=False))
