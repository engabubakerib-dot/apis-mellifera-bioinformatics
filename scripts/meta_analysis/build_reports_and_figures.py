#!/usr/bin/env python3
"""Builds reports/*.xlsx, results/meta_analysis/*, results/figures/* from the saved result tables."""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
R='.'; os.makedirs('results/meta_analysis',exist_ok=True); os.makedirs('results/figures',exist_ok=True)
g78=pd.read_csv('results/differential_expression/GSE178742_targets_colony_level_stats.tsv',sep='\t')
g68=pd.read_csv('results/differential_expression/GSE168740_targets_per_age.tsv',sep='\t')
g12=pd.read_csv('results/differential_expression/GSE121885_targets_factorial.tsv',sep='\t')
tgt=pd.read_csv('data/metadata/gene_targets_ncbi.tsv',sep='\t')
# ---------------- effect-size table (one row per gene x contrast) ----------------
rows=[]
lab={'ApiaryForagers':'GSE178742 apiary foragers (clothianidin, colony n=3v3)','ApiaryNurses':'GSE178742 apiary nurses (clothianidin, colony n=3v3)',
     'FieldForagers':'GSE178742 field foragers (NNI field, colony n=5v5)','FieldNurses':'GSE178742 field nurses (NNI field, colony n=5v4)'}
for _,r in g78.iterrows():
    if pd.isna(r.logFC): continue
    rows.append(dict(sym=r.sym,contrast=lab[r.exp],study='GSE178742',stressor='neonicotinoid (clothianidin)',tissue='brain',scale='log2 (edgeR logCPM)',log2FC=r.logFC,se=r.se,df=(4 if r.exp.startswith('Apiary') else (8 if r.exp=='FieldForagers' else 7)),p=r.p,fdr=r.fdr_genomewide,independence_group='GSE178742-'+('apiary' if r.exp.startswith('Apiary') else 'field')))
for _,r in g68.iterrows():
    try: l=float(r.log2FC_trt_vs_ctrl); s=float(r.se)
    except: continue
    rows.append(dict(sym=r.sym,contrast=f'GSE168740 dinotefuran {r.age} (pooled libs n=3v3)',study='GSE168740',stressor='neonicotinoid (dinotefuran)',tissue='brain (15-brain pools)',scale='log2(FPKM+1)',log2FC=l,se=s,df=4,p=float(r.p),fdr=float(r.fdr_genomewide),independence_group='GSE168740-'+r.age))
for _,r in g12.iterrows():
    rows.append(dict(sym=r.sym,contrast='GSE121885 diet Rockrose vs Chestnut (main effect, n=12v12)',study='GSE121885',stressor='nutrition (pollen source)',tissue='whole body',scale='log2(norm count+1)',log2FC=r.diet_log2FC_R_vs_C,se=r.diet_se,df=20,p=r.diet_p,fdr=r.diet_fdr,independence_group='GSE121885'))
E=pd.DataFrame(rows); E['ci_low']=E.log2FC-stats.t.ppf(.975,E.df)*E.se; E['ci_high']=E.log2FC+stats.t.ppf(.975,E.df)*E.se
E.to_csv('results/meta_analysis/effect_sizes_all_targets.tsv',sep='\t',index=False)
# ---------------- heterogeneity / random-effects (brain-neonicotinoid stratum only) ----------------
def dl(y,s):
    w=1/s**2; ybar=(w*y).sum()/w.sum(); Q=(w*(y-ybar)**2).sum(); k=len(y); c=w.sum()-(w**2).sum()/w.sum()
    tau2=max(0,(Q-(k-1))/c); ws=1/(s**2+tau2); mu=(ws*y).sum()/ws.sum(); se=np.sqrt(1/ws.sum()); I2=max(0,(Q-(k-1))/Q) if Q>0 else 0
    return mu,se,Q,k,I2,tau2
ma=[]
for sym,g in E[E.tissue.str.startswith('brain')].groupby('sym'):
    if len(g)<3: continue
    mu,se,Q,k,I2,tau2=dl(g.log2FC.values,g.se.values)
    # sensitivity: one contrast per independent cohort group is not possible within GSE178742 caste pairs (same colonies) -> report both
    ma.append(dict(sym=sym,k_contrasts=k,pooled_log2FC_DL=mu,se=se,ci_low=mu-1.96*se,ci_high=mu+1.96*se,p=2*stats.norm.sf(abs(mu/se)),Q=Q,Q_p=stats.chi2.sf(Q,k-1),I2=I2,tau2=tau2,
        n_up=int((g.log2FC>0).sum()),n_down=int((g.log2FC<0).sum())))
M=pd.DataFrame(ma)
M['note']='EXPLORATORY: contrasts not independent (GSE178742 forager/nurse share colonies; GSE168740 ages share source colonies); mixed scales; do not cite as a definitive pooled effect'
M['n_nominal_p_lt_0.05']=[int((E[(E.sym==s)&E.tissue.str.startswith('brain')].p<0.05).sum()) for s in M.sym]
M.sort_values('Q_p').to_csv('results/meta_analysis/brain_neonicotinoid_random_effects_EXPLORATORY.tsv',sep='\t',index=False)
# ---------------- figures ----------------
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
def forest(sym,fn,title):
    g=E[E.sym==sym].reset_index(drop=True); fig,ax=plt.subplots(figsize=(7.2,0.42*len(g)+1.2))
    for i,r in g.iterrows():
        c='#B5522B' if r.study=='GSE121885' else ('#2B6CB0' if r.study=='GSE178742' else '#2F855A')
        ax.plot([r.ci_low,r.ci_high],[i,i],color=c,lw=1.6); ax.plot(r.log2FC,i,'s',color=c,ms=5)
        ax.text(1.02,i,f'p={r.p:.3f}',transform=ax.get_yaxis_transform(),va='center',fontsize=8)
    ax.axvline(0,color='#888',lw=.8,ls='--'); ax.set_yticks(range(len(g))); ax.set_yticklabels(g.contrast,fontsize=8); ax.invert_yaxis()
    ax.set_xlabel('log2 fold change (stressed vs control; 95% CI, t-based)'); fig.suptitle(title,fontsize=10,x=0.01,ha='left')
    fig.tight_layout(); fig.savefig(fn,dpi=200); plt.close(fig)
forest('Vg','results/figures/Fig1_Vg_forest_by_contrast.png','Vitellogenin (Vg; GeneID 406088 / GB49544) across stressor contrasts')
forest('Mrjp1','results/figures/Fig2_Mrjp1_forest_by_contrast.png','Mrjp1 (GeneID 406090 / GB55205) across stressor contrasts')
# carryover figure
P=pd.read_csv('results/gene_expression/GSE178742_key_genes_logCPM_per_sample.csv')
fig,axs=plt.subplots(1,2,figsize=(8.4,3.6))
for e,c in zip(['Apiary_Foragers','Apiary_Nurses','Field_Foragers','Field_Nurses'],['#2B6CB0','#B5522B','#2F855A','#6B46C1']):
    g=P[P.exp==e]; axs[0].scatter(g.GB53576,g.GB55205,s=14,color=c,label=e.replace('_',' ')); axs[1].scatter(g.GB51698,g.GB49544,s=14,color=c)
axs[0].set_xlabel('apisimin logCPM (HPG marker)'); axs[0].set_ylabel('Mrjp1 logCPM'); axs[0].legend(fontsize=7,frameon=False)
axs[1].set_xlabel('Hex70a logCPM (fat-body storage protein)'); axs[1].set_ylabel('Vg logCPM')
fig.suptitle('GSE178742 brain libraries: MRJP and Vg signals track gland / fat-body marker genes',fontsize=10,x=0.02,ha='left')
fig.tight_layout(); fig.savefig('results/figures/Fig3_GSE178742_tissue_carryover_markers.png',dpi=200); plt.close(fig)
# GSE168740 per-sample
F=pd.read_csv('results/gene_expression/GSE168740_key_genes_FPKM_per_sample.tsv',sep='\t'); F['age']=F['sample'].str.split('_').str[1]; F['grp']=F['sample'].str[0]
fig,axs=plt.subplots(1,3,figsize=(9,3.2))
for ax,gn in zip(axs,['Mrjp1','Vg','Apisimin']):
    for j,age in enumerate(['1d','5d','10d']):
        for k,(gr,c) in enumerate([('C','#4A5568'),('T','#C05621')]):
            v=np.log2(F[(F.age==age)&(F.grp==gr)][gn]+1); ax.scatter(np.full(len(v),j*3+k)+np.linspace(-.15,.15,len(v)),v,color=c,s=16)
    ax.set_xticks([0.5,3.5,6.5]); ax.set_xticklabels(['1 d','5 d','10 d']); ax.set_title(gn,fontsize=9); ax.set_ylabel('log2(FPKM+1)')
axs[0].scatter([],[],color='#4A5568',label='control'); axs[0].scatter([],[],color='#C05621',label='dinotefuran'); axs[0].legend(fontsize=7,frameon=False)
fig.suptitle('GSE168740 (all 18 libraries): per-library values; note C1_1d outlier',fontsize=10,x=0.02,ha='left'); fig.tight_layout()
fig.savefig('results/figures/Fig4_GSE168740_Mrjp1_Vg_apisimin_per_library.png',dpi=200); plt.close(fig)
# gene-ID mapping workbook: now built by scripts/gene_mapping/build_mapping_workbook (see v2 notes in PROJECT_STATUS)
# ---------------- dataset audit workbook ----------------
A=pd.DataFrame([
 dict(accession='GSE168740 / PRJNA713798 / SRP310279',title='NGS of honey bee brain transcriptomes treated by dinotefuran',publication='BMC Genomics 2021 (PMID 34217210; PMC8254963)',subspecies='not stated (Hangzhou, China)',tissue='brain (15 pooled brains/library)',stressor='dinotefuran 0.01 mg/L oral, chronic 1/5/10 d (lab cage)',control='50% sucrose',n='18 libraries = 3 ages x 2 trt x 3 pools (VERIFIED from GEO SOFT)',experimental_unit='pooled library (source: brood frames from 3 colonies; colony not tracked per library)',age='newly emerged, exposed 1/5/10 d',caste='young workers (caged)',platform='HiSeq 2500 PE150',raw_reads='yes (SRA)',processed='FPKM per sample (Cufflinks); duplicate Gene_ID rows (14,480/file) all identical',annotation='GCF_003254395.2 (Amel_HAv3.1), NCBI GeneIDs',confounders='pooling; C1_1d outlier profile; age = exposure duration; massive MRJP/apisimin FPKM in "brain" suggests gland carryover',DE_suitability='exploratory only from FPKM (n=3v3); raw-count reanalysis recommended',meta_suitability='brain/neonic stratum, with caveats',status='ANALYSED (all 18 libraries)'),
 dict(accession='GSE178742 / PRJNA740269 / SRP325318',title='Searching beyond the streetlight: neonicotinoid exposure alters the neurogenomic state of worker honey bees',publication='Ecology & Evolution 2021/2022 (PMID 35003705; PMC8717355)',subspecies='not stated (Ontario, Canada)',tissue='brain',stressor='Field: NNI-treated corn proximity (2014); Apiary: clothianidin pollen patties 4.9->2.0 ppb, 12 wk (collected day 30)',control='Field: site >=3 km from agriculture; Apiary: unspiked patties',n='79 GEO samples: apiary 30 foragers + 30 nurses (3 exposed + 3 control colonies x 5 bees); field 10 foragers + 9 nurses (5 pooled samples/colony per caste; 1 nurse sample lost) (VERIFIED)',experimental_unit='colony (apiary: bees nested in colony -> n=3 v 3; field: 5 v 5 / 5 v 4)',age='not controlled (behavioural castes)',caste='nurse, forager',platform='HiSeq 4000 PE100',raw_reads='yes (SRA)',processed='edgeR TMM logCPM matrices (4)',annotation='Amel_4.5 + OGS v3.2 (GB IDs)',confounders='field: exposure fully confounded with site (2 sites); apiary: series text says 10 colonies but only 6 colonies sequenced; matrix aliases C25-F1.1->C25-F10, C27-F1->C27-F10',DE_suitability='yes for apiary at colony level (low power, n=3v3); field confounded',meta_suitability='yes (brain/neonic), forager & nurse contrasts non-independent',status='ANALYSED (processed logCPM)'),
 dict(accession='GSE216021 / PRJNA891719',title='Genetics of tolerance in honey bees to the neonicotinoid clothianidin',publication='PMID 36843853 (2023)',subspecies='not stated',tissue='brain, Malpighian tubules, ventriculus',stressor='clothianidin 4.27 ppb, all bees exposed (24 h)',control='NONE unexposed; contrast = tolerant vs susceptible patrilines',n='18 = 6 patrilines x 3 tissues (VERIFIED)',experimental_unit='patriline',age='8 d',caste='worker',platform='NovaSeq 6000',raw_reads='yes',processed='per-sample counts (HTSeq) in RAW tar',annotation='Amel_4.5',confounders='genotype only; no exposure contrast',DE_suitability='NOT for exposure effects; usable for tolerance/tissue baseline only (n=1 per patriline-tissue)',meta_suitability='no',status='AUDITED - excluded from exposure meta-analysis'),
 dict(accession='GSE121885 / PRJNA498919',title='Transcriptomic responses to diet quality and viral infection in Apis mellifera',publication='PMIDs 31117959, 32542031 (Dolezal et al.)',subspecies='not stated',tissue='whole body',stressor='pollen diet (Rockrose vs Chestnut) x IAPV infection',control='factorial',n='24 = 4 groups x 6 (VERIFIED)',experimental_unit='individual bee per GEO; cage structure not in GEO metadata - verify in paper',age='not in GEO metadata',caste='worker (caged)',platform='HiSeq 2500',raw_reads='yes',processed='raw HTSeq counts (OGS3.2 GB IDs)',annotation='BeeBase OGS v3.2',confounders='possible cage effects; diet-quality ranking of the two pollens not verified here',DE_suitability='yes (raw counts)',meta_suitability='nutrition stratum (whole body; not poolable with brain)',status='ANALYSED (2x2 factorial)'),
 dict(accession='GSE71750 / PRJNA292006',title='Regulation of hypopharyngeal gland activity under nutrition and brood pheromone',publication='not linked in GEO',subspecies='not stated',tissue='hypopharyngeal gland',stressor='pollen vs honey diet; brood pheromone',control='Pollen (per GEO titles)',n='6 = 3 groups x 2 (VERIFIED)',experimental_unit='unclear',age='GEO lists day 9-14 per sample (likely mis-entered)',caste='worker',platform='HiSeq 2500',raw_reads='yes',processed='FPKM (DFCI gene index - legacy)',annotation='DFCI Apis mellifera gene index (legacy)',confounders='n=2; legacy annotation; inconsistent age metadata',DE_suitability='descriptive only',meta_suitability='no (n=2)',status='AUDITED'),
 dict(accession='GSE300082 / PRJNA1278995',title='Imidacloprid susceptibility is linked to life history regulation in honey bees',publication='PMID 41814586 (2026)',subspecies='not stated',tissue='whole bee',stressor='none at sequencing (all samples "untreated"); susceptibility measured separately per colony',control='-',n='96 (VERIFIED)',experimental_unit='colony of origin',age='Day 7 etc.',caste='nurse, forager',platform='NovaSeq X',raw_reads='yes',processed='DESeq2 count table',annotation='GCF_003254395.2',confounders='-',DE_suitability='not an exposure contrast; could test colony susceptibility-associated baseline expression of Vg/MRJP',meta_suitability='no (different question)',status='AUDITED'),
 dict(accession='PRJNA1221518',title='Transcriptomic data of pupal A. m. ligustica worker post cold-stress',publication='not identified',subspecies='A. m. ligustica',tissue='pupa (whole?)',stressor='cold stress',control='unknown - sample titles empty',n='9 runs (ENA)',experimental_unit='unknown',age='pupa',caste='worker',platform='NovaSeq 6000 PE',raw_reads='yes (~2.4-3.9 GB/run)',processed='none',annotation='-',confounders='design labels missing in ENA',DE_suitability='only after raw-read processing and design recovery from BioSample/paper',meta_suitability='cold stratum candidate (pupal - not comparable to adults)',status='AUDITED - deferred (raw reads only)'),
 dict(accession='PRJNA495845',title='Apis mellifera raw sequence reads (whole adult workers)',publication='not identified',subspecies='not stated',tissue='whole adult',stressor='presumed transport/temperature (unverified)',control='unknown',n='>=26 runs (ENA, first page)',experimental_unit='hive/day labels in some titles',age='-',caste='worker',platform='NextSeq 500 PE',raw_reads='yes',processed='none',annotation='-',confounders='design labels mostly generic',DE_suitability='needs BioSample attributes + raw processing',meta_suitability='unknown',status='AUDITED - deferred'),
])
with pd.ExcelWriter('reports/Dataset_Audit.xlsx') as w:
    A.to_excel(w,'audit',index=False)
    pd.read_csv('data/metadata/GSE178742_sample_metadata_verified.tsv',sep='\t').to_excel(w,'GSE178742_samples',index=False)
    pd.read_csv('data/metadata/GSE168740_sample_metadata_verified.tsv',sep='\t').to_excel(w,'GSE168740_samples',index=False)
    pd.read_csv('data/metadata/GSE121885_sample_metadata.tsv',sep='\t').to_excel(w,'GSE121885_samples',index=False)
    pd.read_csv('results/QC/genomewide_and_QC_summary.tsv',sep='\t').to_excel(w,'genomewide_summary',index=False)
print(M[M.sym.isin(['Vg','Mrjp1','Mrjp3','Cat','Sod1','CYP9Q1','Hsp90','GstS4'])][['sym','k_contrasts','pooled_log2FC_DL','ci_low','ci_high','p','I2','Q_p','n_up','n_down','n_nominal_p_lt_0.05']].round(3).to_string(index=False))
