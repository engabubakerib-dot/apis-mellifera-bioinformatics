#!/usr/bin/env python3
"""Compares genome-wide tables regenerated from LOCAL downloaded files (reanalyse_processed_geo.py) with the
stored target-gene tables produced earlier (browser computation). Writes results/validation/reproduction_check.tsv"""
import os,sys,numpy as np,pandas as pd
R=sys.argv[1] if len(sys.argv)>1 else '.'; D=f'{R}/results/differential_expression'; os.makedirs(f'{R}/results/validation',exist_ok=True); rows=[]
old=pd.read_csv(f'{D}/GSE178742_targets_colony_level_stats.tsv',sep='\t')
for exp,lab in [('Apiary_Foragers','ApiaryForagers'),('Apiary_Nurses','ApiaryNurses'),('Field_Foragers','FieldForagers'),('Field_Nurses','FieldNurses')]:
    f=f'{D}/GSE178742_{exp}_genomewide.tsv.gz'
    if not os.path.exists(f): rows.append(dict(dataset='GSE178742',contrast=exp,status='MISSING regenerated file')); continue
    n=pd.read_csv(f,sep='\t',index_col=0); o=old[(old.exp==lab)&old.logFC.notna()].set_index('gb')
    j=o.join(n,rsuffix='_new',how='inner')
    rows.append(dict(dataset='GSE178742',contrast=exp,genes_compared=len(j),max_abs_diff_logFC=float((j.logFC-j.logFC_new).abs().max()),max_abs_diff_p=float((j.p-j.p_new).abs().max()),
                     n_FDR_lt_005=int((n.fdr<0.05).sum()),n_p_lt_001=int((n.p<0.01).sum()),status='OK' if (j.logFC-j.logFC_new).abs().max()<0.002 else 'DIFFERS'))
old=pd.read_csv(f'{D}/GSE168740_targets_per_age.tsv',sep='\t',dtype=str)
for age in ['1d','5d','10d']:
    f=f'{D}/GSE168740_{age}_genomewide.tsv.gz'
    if not os.path.exists(f): rows.append(dict(dataset='GSE168740',contrast=age,status='MISSING')); continue
    n=pd.read_csv(f,sep='\t',index_col=0); n.index=n.index.astype(str); o=old[(old.age==age)].copy(); o=o[pd.to_numeric(o.log2FC_trt_vs_ctrl,errors='coerce').notna()].set_index('ncbi_gene_id')
    j=o.join(n,how='inner',rsuffix='_new'); d=(pd.to_numeric(j.log2FC_trt_vs_ctrl)-j.log2FC).abs()
    rows.append(dict(dataset='GSE168740',contrast=age,genes_compared=len(j),max_abs_diff_logFC=float(d.max()),max_abs_diff_p=float((pd.to_numeric(j.p)-j.p_new).abs().max()),n_FDR_lt_005=int((n.fdr<0.05).sum()),n_p_lt_001=int((n.p<0.01).sum()),status='OK' if d.max()<0.002 else 'DIFFERS'))
f=f'{D}/GSE121885_factorial_genomewide.tsv.gz'
if os.path.exists(f):
    n=pd.read_csv(f,sep='\t',index_col=0); o=pd.read_csv(f'{D}/GSE121885_targets_factorial.tsv',sep='\t').set_index('gb'); j=o.join(n,how='inner',rsuffix='_new'); d=(j.diet_log2FC_R_vs_C-j.diet_log2FC).abs()
    rows.append(dict(dataset='GSE121885',contrast='diet',genes_compared=len(j),max_abs_diff_logFC=float(d.max()),max_abs_diff_p=float((j.diet_p-j.diet_p_new).abs().max()),n_FDR_lt_005=int((n.diet_fdr<0.05).sum()),status='OK' if d.max()<0.002 else 'DIFFERS'))
out=pd.DataFrame(rows); out.to_csv(f'{R}/results/validation/reproduction_check.tsv',sep='\t',index=False); print(out.to_string())
