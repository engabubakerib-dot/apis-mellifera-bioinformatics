#!/usr/bin/env python3
"""Reproduces, from LOCAL downloaded files only, the browser-derived secondary analyses:
(1) GSE168740 library-library correlation (C1_1d), (2) paper-gene directions, (3) Vg/Mrjp1 co-expression with
markers, (4) GSE168740 1-day contrast without C1_1d. Output: results/validation/secondary_reproduction.tsv"""
import tarfile,gzip,io,sys,numpy as np,pandas as pd
from scipy import stats
R=sys.argv[1] if len(sys.argv)>1 else '.'; rows=[]
cols={}
with tarfile.open(f'{R}/data/raw/GSE168740/GSE168740_RAW.tar') as t:
    for m in t.getmembers():
        d=pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(m).read())),sep='\t',dtype={'Gene_ID':str}); s=d.columns[2]; cols[s]=d.groupby('Gene_ID')[s].first()
F=pd.DataFrame(cols); keep=F[(F>=1).all(1)]; L=np.log2(keep+1); C=L.corr()
for s in ['C1_1d']+[c for c in C if c.endswith('_1d') and c!='C1_1d']:
    top=C[s].drop(s).sort_values(ascending=False).head(3); rows.append(dict(analysis='GSE168740 library correlation',item=s,value='; '.join(f'{k}={v:.3f}' for k,v in top.items()),n_genes=len(keep)))
for sym,g in [('Hbg3','406131'),('Mrjp1','406090'),('LOC408790','408790'),('Def1','406143'),('GstS4','411045'),('LOC113219351','113219351')]:
    for age in ['1d','5d','10d']:
        c=F.loc[g,[f'C{i}_{age}' for i in (1,2,3)]]; tt=F.loc[g,[f'T{i}_{age}' for i in (1,2,3)]]
        rows.append(dict(analysis='GSE168740 FPKM per library',item=f'{sym} {age}',value='C='+'/'.join(f'{v:.1f}' for v in c)+' T='+'/'.join(f'{v:.1f}' for v in tt)))
c=['C2_1d','C3_1d']; tr=['T1_1d','T2_1d','T3_1d']; k=F[c+tr].mean(1)>=1; X=np.log2(F.loc[k]+1); t,p=stats.ttest_ind(X[tr],X[c],axis=1)
rows.append(dict(analysis='GSE168740 1d without C1_1d (2 v 3, Student t)',item='genes p<0.01 / tested',value=f'{int((p<0.01).sum())}/{int(k.sum())} (expected {k.sum()*0.01:.0f})'))
for sym,g in [('Vg','406088'),('Mrjp1','406090'),('Hymenoptaecin','406142'),('GstS4','411045')]:
    i=list(X.index).index(g) if g in X.index else None
    if i is not None: rows.append(dict(analysis='GSE168740 1d without C1_1d',item=sym,value=f'log2FC={X.loc[g,tr].mean()-X.loc[g,c].mean():.3f} p={p[i]:.3f}'))
def resid(df,grp): return df.sub(df.T.groupby(grp).transform('mean').T)
for exp in ['Apiary_Nurses','Apiary_Foragers']:
    M=pd.read_csv(f'{R}/data/raw/GSE178742/GSE178742_{exp}_Logcpm.txt.gz',sep='\t',index_col=0); M.columns=[x.strip() for x in M.columns]; Rz=resid(M,[x.split('-')[0] for x in M.columns])
    for a,b in [('GB49544','GB51698'),('GB55205','GB53576'),('GB49544','GB55205')]:
        rows.append(dict(analysis=f'GSE178742 {exp} colony-residualised r',item=f'{a}~{b}',value=f'{np.corrcoef(Rz.loc[a],Rz.loc[b])[0,1]:.3f}'))
Cn=pd.read_csv(f'{R}/data/raw/GSE121885/GSE121885_processedData.txt.gz',sep='\t',index_col=0); Cn=Cn[~Cn.index.str.startswith('__')]
pos=Cn[(Cn>0).all(1)]; lg=np.log(pos); sf=np.exp(lg.sub(lg.mean(1),axis=0).median(0)); N=Cn/sf; N=N[N.mean(1)>=10]; Y=np.log2(N+1); Rz=resid(Y,[x.split('.')[0] for x in Y.columns])
for a,b in [('GB49544','GB51698'),('GB55205','GB53576'),('GB49544','GB55205'),('GB49544','GB53755')]:
    rows.append(dict(analysis='GSE121885 group-residualised r',item=f'{a}~{b}',value=f'{np.corrcoef(Rz.loc[a],Rz.loc[b])[0,1]:.3f}'))
out=pd.DataFrame(rows); out.to_csv(f'{R}/results/validation/secondary_reproduction.tsv',sep='\t',index=False); print(out.to_string())
