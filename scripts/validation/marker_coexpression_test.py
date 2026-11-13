#!/usr/bin/env python3
"""Tissue-marker co-expression test (evidence type: CO-EXPRESSION, not proof of contamination).
Gland markers (HPG): apisimin, Hbg3 (alpha-glucosidase), Acph-1; fat-body: Hex70a, Hex110; neuronal: Syt1, brp.
IDs verified via NCBI Gene aliases (OGS3.2 GB / GeneID). Uses local GEO files only.
Prediction under the carry-over hypothesis: MRJP index correlates positively with ALL independent HPG markers and
negatively with neuronal markers (dilution); under intrinsic neural expression: no consistent relation to HPG markers."""
import tarfile,gzip,io,sys,numpy as np,pandas as pd
from scipy import stats
R=sys.argv[1] if len(sys.argv)>1 else '.'; rows=[]
GB=dict(Mrjp1='GB55205',Mrjp2='GB55212',Mrjp3='GB55204',Mrjp4='GB55206',Mrjp5='GB55208',Mrjp7='GB55213',apisimin='GB53576',Hbg3='GB43247',Acph1='GB41326',Hex70a='GB51698',Hex110='GB44996',Vg='GB49544',Syt1='GB54827',brp='GB42866')
GID=dict(Mrjp1='406090',Mrjp2='406091',Mrjp3='406121',Mrjp4='406133',Mrjp5='406116',Mrjp7='409555',apisimin='406093',Hbg3='406131',Acph1='725215',Hex70a='726848',Hex110='551648',Vg='406088',Syt1='408613',brp='409012')
def test(X,lab,unit_note):
    mrjp=X[['Mrjp1','Mrjp2','Mrjp3','Mrjp4','Mrjp7']].mean(1)
    for m in ['apisimin','Hbg3','Acph1','Hex70a','Vg','Syt1','brp']:
        if m in X and X[m].notna().all():
            r,p=stats.spearmanr(mrjp,X[m]); rows.append(dict(dataset=lab,x='MRJP index (Mrjp1,2,3,4,7)',y=m,spearman_rho=round(r,3),p=p,n=len(X),note=unit_note))
    for m in ['Hex70a','Hex110','Syt1','brp']:
        if m in X: r,p=stats.spearmanr(X['Vg'],X[m]); rows.append(dict(dataset=lab,x='Vg',y=m,spearman_rho=round(r,3),p=p,n=len(X),note=unit_note))
for exp in ['Apiary_Nurses','Apiary_Foragers','Field_Nurses','Field_Foragers']:
    M=pd.read_csv(f'{R}/data/raw/GSE178742/GSE178742_{exp}_Logcpm.txt.gz',sep='\t',index_col=0); M.columns=[c.strip() for c in M.columns]
    have={k:v for k,v in GB.items() if v in M.index}; X=M.loc[list(have.values())].T; X.columns=list(have.keys())
    if exp.startswith('Apiary'):
        col=[c.split('-')[0] for c in X.index]; X=X-X.groupby(col).transform('mean'); note='individual bees, colony-mean-centred (within-colony association)'
    else: note='colony-pooled samples (exposure-confounded site not removed)'
    test(X,'GSE178742 '+exp,note); 
    missing=[k for k in GB if k not in have]; 
    if missing: rows.append(dict(dataset='GSE178742 '+exp,x='MISSING IDs (filtered from processed matrix)',y=','.join(missing)))
cols={}
with tarfile.open(f'{R}/data/raw/GSE168740/GSE168740_RAW.tar') as t:
    for m in t.getmembers():
        d=pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(m).read())),sep='\t',dtype={'Gene_ID':str}); s=d.columns[2]; cols[s]=d.groupby('Gene_ID')[s].first()
F=pd.DataFrame(cols); X=np.log2(F.loc[list(GID.values())].T+1); X.columns=list(GID.keys())
test(X,'GSE168740 all 18 libraries','pooled libraries; treatment and age NOT removed')
age=[s.split('_')[1]+s[0] for s in X.index]; Xc=X-X.groupby(age).transform('mean'); test(Xc,'GSE168740 group-centred','centred within age x treatment (n=3 per group)')
out=pd.DataFrame(rows); out['p']=out.p.astype(float).map(lambda v: f'{v:.2e}' if pd.notna(v) else ''); out.to_csv(f'{R}/results/validation/marker_coexpression_test.tsv',sep='\t',index=False); print(out.to_string())
