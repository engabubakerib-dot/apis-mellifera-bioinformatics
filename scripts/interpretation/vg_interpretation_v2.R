## Vg-centred interpretation tables and publication figures - REVISION 2 (2026-11-13).
## Implements the nine revisions of reports/VG_FIGURES_POST_RUN_AUDIT_2026-11-13.md. Version 1 (vg_interpretation.R) and its outputs are kept unchanged.
## Reads ONLY the validated DE run results/rawreads/DE/run_20261113_044101/ and data/metadata/gene_targets_ncbi.tsv. Nothing is simulated, re-fitted or
## re-tested: every statistic is copied from an edgeR/DESeq2 output table; Fig. 1 per-library values use the run's own count matrices and TMM factors.
## Writes a NEW folder results/rawreads/Vg_interpretation/run_v2_<time>/ (tables/, figures/ as PDF + 300-dpi LZW TIFF, logs). Installs nothing.
## Run from the project folder:  Rscript --vanilla scripts/interpretation/vg_interpretation_v2.R   (or double-click RUN_VG_FIGURES_V2.bat)
rl=file.path('envs','R-library',paste(R.version$major,sub('\\..*','',R.version$minor),sep='.')); if(dir.exists(rl)) .libPaths(c(normalizePath(rl),.libPaths()))
for(p in c('ggplot2','edgeR')) if(!requireNamespace(p,quietly=TRUE)) stop('required package not installed: ',p,' (nothing is installed by this script)')
suppressMessages({library(ggplot2);library(edgeR)})
HAS_REPEL=requireNamespace('ggrepel',quietly=TRUE)
DE='results/rawreads/DE/run_20261113_044101'; if(!file.exists(file.path(DE,'_DE_FINISHED'))) stop('validated DE run not found: ',DE)
OUT=file.path('results/rawreads/Vg_interpretation',format(Sys.time(),'run_v2_%Y%m%d_%H%M%S')); if(dir.exists(OUT)) stop('output folder exists: ',OUT)
TD=file.path(OUT,'tables'); FD=file.path(OUT,'figures'); dir.create(TD,recursive=TRUE); dir.create(FD,recursive=TRUE)
LOG=file.path(OUT,'RUN_LOG.txt'); say=function(...){m=paste0(...); cat(m,'\n'); cat(m,'\n',file=LOG,append=TRUE)}
STATUS=character(0); step=function(name,expr){r=tryCatch({force(expr);'OK'},error=function(e) paste('ERROR:',conditionMessage(e))); STATUS[name]<<-r; say(sprintf('%-48s %s',name,r))}
say('vg_interpretation_v2.R  ',format(Sys.time()),'  ',R.version.string,'  ggplot2 ',as.character(packageVersion('ggplot2')),'  edgeR ',as.character(packageVersion('edgeR')),
    '  ggrepel ',if(HAS_REPEL) as.character(packageVersion('ggrepel')) else 'not installed (built-in label spacing used)')
say('DE input folder: ',DE)
invisible(file.copy('scripts/interpretation/vg_interpretation_v2.R',file.path(OUT,'vg_interpretation_v2.R.used')))
VG='GB49544'
## ---------------- consistent styling ----------------
COL_TIER=c('Primary'='#2a78d6','Sensitivity'='#eb6834','Secondary method'='#e87ba4','Exploratory'='#1baf7a','Diagnostic/QC'='#7d7c78')
TIER_ORDER=names(COL_TIER)
COL_DIR=c('Up'='#e34948','Down'='#2a78d6','Not significant'='#c9c8c2')
TH=theme_bw(base_size=8,base_family='sans')+theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(linewidth=0.25,colour='#e6e5e1'),
   strip.background=element_rect(fill='#f0efec',colour='#b9b8b3',linewidth=0.3),strip.text=element_text(size=7,face='bold',margin=margin(2,2,2,2)),
   panel.border=element_rect(colour='#8f8e89',linewidth=0.4),axis.ticks=element_line(linewidth=0.3),axis.text=element_text(size=6.5,colour='#2b2b29'),
   axis.title=element_text(size=7.5),legend.key.size=unit(3.5,'mm'),legend.title=element_text(size=7,face='bold'),legend.text=element_text(size=6.8),
   plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7,colour='#52514e',lineheight=1.05),
   plot.caption=element_text(size=6.3,colour='#52514e',hjust=0,lineheight=1.05),plot.title.position='plot',plot.caption.position='plot',
   plot.margin=margin(6,8,6,6),legend.box.just='left')
WR=function(s,n) paste(strwrap(s,n),collapse='\n')          # wrap long subtitles / captions so nothing is clipped
W_SUB=120; W_CAP=140
TIFTYPE=if(isTRUE(capabilities('cairo'))) 'cairo' else NULL
save_fig=function(p,name,w,h){  # w, h in mm; PDF (editable text, no dingbats) + 300-dpi LZW TIFF
  ggsave(file.path(FD,paste0(name,'.pdf')),p,width=w,height=h,units='mm',device='pdf',useDingbats=FALSE)
  if(is.null(TIFTYPE)) ggsave(file.path(FD,paste0(name,'.tiff')),p,width=w,height=h,units='mm',dpi=300,device='tiff',compression='lzw')
  else ggsave(file.path(FD,paste0(name,'.tiff')),p,width=w,height=h,units='mm',dpi=300,device='tiff',compression='lzw',type=TIFTYPE)}
fmtp=function(x) ifelse(is.na(x),'NA',ifelse(x<0.001,formatC(x,format='e',digits=1),formatC(x,format='f',digits=3)))
TIER_NOTE=paste('Tiers: Primary = pre-specified main model; Sensitivity = same question with a different sample set (PRJNA1017469: all libraries of the route, 24 ingestion and 23 contact,',
  'QC exclusions reversed; GSE168740: 1 d without library C1_1d); Secondary method = DESeq2 on the same GSE168740 data;',
  'Exploratory = PRJNA1017469 route difference (ingestion effect minus contact effect); Diagnostic/QC = control 4 h vs 1 h, not a treatment effect.')
## ---------------- inputs ----------------
ss=read.delim(file.path(DE,'PRJNA1017469_sample_sheet_and_QC.tsv'),colClasses=c(run='character'))
tg=read.delim('data/metadata/gene_targets_ncbi.tsv',colClasses='character'); tg=tg[nzchar(tg$ogs32_id),]
files=sort(list.files(DE,pattern='_(edgeR_QLF|DESeq2)\\.tsv$')); AN=sub('\\.tsv$','',files)
vgt=read.delim(file.path(DE,'Vg_GB49544_all_analyses.tsv'),colClasses=c(analysis='character'))
if(length(AN)!=40 || !setequal(AN,vgt$analysis)) stop('expected the 40 analyses listed in Vg_GB49544_all_analyses.tsv; found ',length(AN))
say('analyses found: ',length(AN))
inp=c(files,'Vg_GB49544_all_analyses.tsv','PRJNA1017469_sample_sheet_and_QC.tsv','GSE168740_subsetA_lengthScaledTPM_counts.tsv','PRJNA1017469_subsetA_lengthScaledTPM_counts.tsv',
      'PRJNA1017469_sensitivity_allLibraries_ingestion_TMM_factors.tsv','PRJNA1017469_sensitivity_allLibraries_contact_TMM_factors.tsv')
write.table(data.frame(file=file.path(DE,inp),bytes=file.size(file.path(DE,inp)),md5=unname(tools::md5sum(file.path(DE,inp)))),file.path(OUT,'INPUT_MANIFEST.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
std=function(a){t=read.delim(file.path(DE,paste0(a,'.tsv')),row.names=1,check.names=FALSE)
  if('logFC' %in% names(t)) data.frame(gene=rownames(t),log2FC=t$logFC,expr=t$logCPM,expr_measure='logCPM',SE=NA_real_,P=t$PValue,FDR=t$FDR)
  else data.frame(gene=rownames(t),log2FC=t$log2FoldChange,expr=t$baseMean,expr_measure='baseMean',SE=t$lfcSE,P=t$pvalue,FDR=t$padj)}
RES=setNames(lapply(AN,std),AN)
## ---------------- analysis metadata (parsed from the file names written by de_from_kallisto.R) ----------------
CK=c('SUB_vs_CTL_1h','D50_vs_CTL_1h','SUB_vs_CTL_4h','D50_vs_CTL_4h','SUBxTime_interaction','D50xTime_interaction','CTL_4h_vs_1h_drift')
CLAB=c(SUB_vs_CTL_1h='SUB vs CTL, 1 h',D50_vs_CTL_1h='D50 vs CTL, 1 h',SUB_vs_CTL_4h='SUB vs CTL, 4 h',D50_vs_CTL_4h='D50 vs CTL, 4 h',
       SUBxTime_interaction='SUB x time interaction',D50xTime_interaction='D50 x time interaction',CTL_4h_vs_1h_drift='CTL 4 h vs 1 h')
CGRP=list(SUB_vs_CTL_1h=c('SUB_1h','CTL_1h'),D50_vs_CTL_1h=c('D50_1h','CTL_1h'),SUB_vs_CTL_4h=c('SUB_4h','CTL_4h'),D50_vs_CTL_4h=c('D50_4h','CTL_4h'),
       SUBxTime_interaction=c('SUB_4h','CTL_4h','SUB_1h','CTL_1h'),D50xTime_interaction=c('D50_4h','CTL_4h','D50_1h','CTL_1h'),CTL_4h_vs_1h_drift=c('CTL_4h','CTL_1h'))
EXK=c('SUB_1h','D50_1h','SUB_4h','D50_4h')
ngrp=function(route,set,groups){s=ss[ss$route %in% route & (set=='all' | ss$decision=='RETAIN'),]; paste(sapply(groups,function(g) paste(g,sum(s$group==g))),collapse=' v ')}
PANELS=c('GSE168740 (brain, GEO annotation)','PRJNA1017469 ingestion (head)','PRJNA1017469 contact (head)','PRJNA1017469 route difference (head)')
meta=function(a){
  if(startsWith(a,'GSE168740')){age=regmatches(a,regexpr('(10d|5d|1d)',a)); noC1=grepl('noC1_1d',a); d2=grepl('DESeq2$',a); A=sub('d',' d',age)
    tier=ifelse(noC1,'Sensitivity',ifelse(d2,'Secondary method','Primary')); meth=ifelse(d2,'DESeq2','edgeR QL')
    return(data.frame(analysis=a,dataset='GSE168740',tissue='brain per GEO annotation (pools of 15)',panel=PANELS[1],route=NA,model=ifelse(noC1,'sensitivity (without C1_1d)','per-age'),
      method=meth,tier=tier,comparison=paste0('Dinotefuran vs control, ',A,ifelse(noC1,' (without C1_1d)','')),key=age,n_groups=ifelse(noC1,'control 2 v treated 3','control 3 v treated 3'),
      label=paste0('Dinotefuran v control, ',A,ifelse(noC1,', no C1_1d',''),' | ',ifelse(d2,'DESeq2','edgeR'),', ',tolower(tier)),
      colname=paste0('GSE168740 ',A,ifelse(noC1,' no C1_1d',''),' ',ifelse(d2,'DESeq2','edgeR')),
      ord=match(age,c('1d','5d','10d'))+ifelse(noC1,0.5,0)))}
  ex=grepl('_EXPLORATORY_',a); mdl=ifelse(ex,'exploratory',ifelse(grepl('_primary_',a),'primary','sensitivity (all libraries of the route)'))
  if(ex){k=sub('.*routeDifference_(SUB|D50)_(1h|4h)_edgeR_QLF$','\\1_\\2',a); d=sub('_.*','',k); t=sub('h$',' h',sub('.*_','',k))
    return(data.frame(analysis=a,dataset='PRJNA1017469',tissue='head (pools of 3)',panel=PANELS[4],route='ingestion - contact',model=mdl,method='edgeR QL',tier='Exploratory',
      comparison=paste0('Route difference of the ',d,' effect at ',t,' (ingestion - contact)'),key=k,
      n_groups=paste0('ingestion ',ngrp('ingestion','primary',c(paste0(d,'_',sub(' ','',t)),paste0('CTL_',sub(' ','',t)))),'; contact ',ngrp('contact','primary',c(paste0(d,'_',sub(' ','',t)),paste0('CTL_',sub(' ','',t))))),
      label=paste0(d,' effect, ',t,' (ingestion - contact) | exploratory'),colname=paste0('Route diff. ',d,' ',t),ord=100+match(k,EXK)))}
  rt=ifelse(grepl('_ingestion_',a),'ingestion','contact'); k=sub(paste0('.*_',rt,'_(.*)_edgeR_QLF$'),'\\1',a); prim=mdl=='primary'
  tier=ifelse(k=='CTL_4h_vs_1h_drift','Diagnostic/QC',ifelse(prim,'Primary','Sensitivity'))
  lab=ifelse(k=='CTL_4h_vs_1h_drift',paste0(CLAB[[k]],' | diagnostic, ',ifelse(prim,'primary model','sensitivity model')),paste0(CLAB[[k]],' | ',ifelse(prim,'primary','sensitivity (all libraries)')))
  data.frame(analysis=a,dataset='PRJNA1017469',tissue='head (pools of 3)',panel=ifelse(rt=='ingestion',PANELS[2],PANELS[3]),route=rt,model=mdl,method='edgeR QL',tier=tier,
    comparison=paste0(CLAB[[k]],ifelse(k=='CTL_4h_vs_1h_drift',' (diagnostic)','')),key=k,n_groups=ngrp(rt,ifelse(prim,'primary','all'),CGRP[[k]]),label=lab,
    colname=paste0(ifelse(rt=='ingestion','Ing. ','Con. '),CLAB[[k]],ifelse(prim,'',' [sens.]')),ord=ifelse(rt=='ingestion',10,20)+match(k,CK)+ifelse(prim,0,0.5))}
M=do.call(rbind,lapply(AN,meta))
if(anyDuplicated(M$colname)||anyDuplicated(paste(M$panel,M$label))) stop('analysis labels are not unique')
M$tier_rank=match(M$tier,TIER_ORDER)
## =============== TABLE 1: Vg across all 40 analyses (revision 7: gene_id and gene_symbol columns) ===============
step('T1 Vg consolidated table',{
  V=do.call(rbind,lapply(AN,function(a){r=RES[[a]]; o=order(r$P); i=match(VG,r$gene)
    data.frame(analysis=a,tested=!is.na(i),log2FC=r$log2FC[i],expression=r$expr[i],expression_measure=r$expr_measure[1],lfcSE_DESeq2=r$SE[i],P=r$P[i],FDR=r$FDR[i],
      rank_by_P=if(is.na(i)) NA else match(VG,r$gene[o]),genes_tested=nrow(r),DEG_FDR05=sum(r$FDR<0.05,na.rm=TRUE))}))
  V=merge(M,V,by='analysis',sort=FALSE)
  V$direction=ifelse(is.na(V$log2FC),'not available',ifelse(V$log2FC>0,'higher in first-named group (up)','lower in first-named group (down)'))
  V$significance=ifelse(is.na(V$FDR),'FDR not available (DESeq2 independent filtering)',ifelse(V$FDR<0.05,'significant (FDR<0.05)',ifelse(V$P<0.05,'not significant (nominal P<0.05 only; FDR>=0.05)','not significant')))
  chk=merge(V[,c('analysis','log2FC','P')],vgt[,c('analysis','logFC','PValue')],by='analysis')
  if(nrow(chk)!=40 || any(abs(chk$log2FC-chk$logFC)>1e-9) || any(abs(chk$P-chk$PValue)>1e-12)) stop('Vg values disagree with Vg_GB49544_all_analyses.tsv')
  V=V[order(match(V$panel,PANELS),V$tier_rank,V$ord),]
  V=data.frame(gene_id=VG,gene_symbol='Vg',NCBI_gene_id=tg$ncbi_gene_id[match(VG,tg$ogs32_id)],V[,setdiff(names(V),c('label','colname','ord','tier_rank'))],check.names=FALSE)
  V<<-V; write.table(V,file.path(TD,'T1_Vg_GB49544_all_40_analyses.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  say('  Vg significant (FDR<0.05) in: ',paste(V$analysis[which(V$FDR<0.05)],collapse=', '))})
VM=merge(V,M[,c('analysis','label','colname','ord','tier_rank')],by='analysis')
## =============== TABLE 2/3: edgeR vs DESeq2 (revision 8: totals + counts among common genes; r to 4 decimals) ===============
step('T2-T3 edgeR vs DESeq2 concordance',{
  K=c('1d','5d','10d','1d_noC1_1d'); CMP<<-do.call(rbind,lapply(K,function(k){e=RES[[paste0('GSE168740_subsetA_',k,'_edgeR_QLF')]]; d=RES[[paste0('GSE168740_subsetA_',k,'_DESeq2')]]
    m=merge(e,d,by='gene',suffixes=c('.edgeR','.DESeq2')); m$age=k; m}))
  S=do.call(rbind,lapply(K,function(k){m=CMP[CMP$age==k,]; e=RES[[paste0('GSE168740_subsetA_',k,'_edgeR_QLF')]]; d=RES[[paste0('GSE168740_subsetA_',k,'_DESeq2')]]
    se=!is.na(m$FDR.edgeR)&m$FDR.edgeR<0.05; sd=!is.na(m$FDR.DESeq2)&m$FDR.DESeq2<0.05; b=se&sd
    data.frame(comparison=k,genes_tested_edgeR=nrow(e),genes_tested_DESeq2=nrow(d),genes_tested_by_both=nrow(m),
      pearson_r_log2FC_common=round(cor(m$log2FC.edgeR,m$log2FC.DESeq2),4),spearman_rho_log2FC_common=round(cor(m$log2FC.edgeR,m$log2FC.DESeq2,method='spearman'),4),
      DEG_edgeR_total=sum(e$FDR<0.05,na.rm=TRUE),DEG_DESeq2_total=sum(d$FDR<0.05,na.rm=TRUE),
      DEG_edgeR_among_common=sum(se),DEG_DESeq2_among_common=sum(sd),DEG_both_among_common=sum(b),DEG_edgeR_only_among_common=sum(se&!sd),DEG_DESeq2_only_among_common=sum(sd&!se),
      same_sign_among_both=sum(sign(m$log2FC.edgeR[b])==sign(m$log2FC.DESeq2[b])),DESeq2_FDR_NA_among_edgeR_DEG_common=sum(se&is.na(m$FDR.DESeq2)))}))
  write.table(S,file.path(TD,'T2_edgeR_vs_DESeq2_genomewide_GSE168740.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  v=CMP[CMP$gene==VG,c('age','log2FC.edgeR','log2FC.DESeq2','P.edgeR','P.DESeq2','FDR.edgeR','FDR.DESeq2')]
  v$agreement=ifelse(is.na(v$FDR.DESeq2),'DESeq2 FDR not available',ifelse((v$FDR.edgeR<0.05)==(v$FDR.DESeq2<0.05),
     ifelse(v$FDR.edgeR<0.05,'both significant','both not significant'),'DISAGREE on significance'))
  v$same_direction=sign(v$log2FC.edgeR)==sign(v$log2FC.DESeq2)
  write.table(data.frame(gene_id=VG,gene_symbol='Vg',v),file.path(TD,'T3_Vg_edgeR_vs_DESeq2_GSE168740.tsv'),sep='\t',quote=FALSE,row.names=FALSE)})
## =============== TABLE 4/5 ===============
step('T4-T5 DEG counts and lists',{
  D=do.call(rbind,lapply(AN,function(a){r=RES[[a]]; s=!is.na(r$FDR)&r$FDR<0.05; data.frame(analysis=a,genes_tested=nrow(r),DEG_FDR05=sum(s),up=sum(s&r$log2FC>0),down=sum(s&r$log2FC<0),FDR_NA=sum(is.na(r$FDR)))}))
  DEGN<<-merge(M,D,by='analysis'); write.table(DEGN[order(DEGN$tier_rank,DEGN$ord),setdiff(names(DEGN),c('label','colname','ord','tier_rank'))],file.path(TD,'T4_DEG_counts_all_40_analyses.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  L=do.call(rbind,lapply(AN,function(a){r=RES[[a]]; r=r[!is.na(r$FDR)&r$FDR<0.05,]; if(!nrow(r)) return(NULL); r$analysis=a; r}))
  i=match(L$gene,tg$ogs32_id); L$documented_symbol=tg$symbol[i]; L$NCBI_gene_id=tg$ncbi_gene_id[i]; L$NCBI_description=tg$ncbi_description[i]; L$project_category=tg$category[i]
  L$annotation_status=ifelse(is.na(i),'no annotation in project resources (OGS v3.2 ID only)','documented: NCBI Gene record (curated project target list)')
  L=merge(M[,c('analysis','dataset','panel','tier','method','comparison')],L,by='analysis'); L=L[order(L$analysis,L$P),]
  write.table(L,file.path(TD,'T5_DEG_FDR05_all_analyses_with_documented_annotation.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  say('  DEG rows: ',nrow(L),'; unique genes: ',length(unique(L$gene)),'; with documented annotation: ',length(unique(L$gene[!is.na(L$documented_symbol)])))})
## =============== TABLE 6/7: curated candidate panel ===============
step('T6-T7 candidate panel',{
  C=do.call(rbind,lapply(AN,function(a){r=RES[[a]]; i=match(tg$ogs32_id,r$gene)
    data.frame(analysis=a,symbol=tg$symbol,ogs32_id=tg$ogs32_id,NCBI_gene_id=tg$ncbi_gene_id,NCBI_description=tg$ncbi_description,category=tg$category,
      tested=!is.na(i),log2FC=r$log2FC[i],expression=r$expr[i],P=r$P[i],FDR=r$FDR[i])}))
  CAND<<-merge(M[,c('analysis','dataset','panel','tier','method','comparison','colname','ord','tier_rank')],C,by='analysis')
  write.table(CAND[,setdiff(names(CAND),c('colname','ord','tier_rank'))],file.path(TD,'T6_candidate_panel_all_analyses.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  h=CAND[!is.na(CAND$FDR)&CAND$FDR<0.05,]; write.table(h[order(h$category,h$symbol),setdiff(names(h),c('colname','ord','tier_rank'))],file.path(TD,'T7_candidate_panel_FDR05_hits.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  say('  candidate rows with FDR<0.05: ',nrow(h))})
## =============== TABLE 8 + FIGURE 1 ===============
step('T8 + Fig1 Vg expression per library',{
  g=read.delim(file.path(DE,'GSE168740_subsetA_lengthScaledTPM_counts.tsv'),row.names=1,check.names=FALSE)
  age=sub('.*_','',colnames(g)); trt=ifelse(substr(colnames(g),1,1)=='C','Control','Dinotefuran'); y=DGEList(as.matrix(g))
  k=filterByExpr(y,group=paste(age,trt)); y2=normLibSizes(y[k,,keep.lib.sizes=FALSE]); eff=y2$samples$lib.size*y2$samples$norm.factors
  G=data.frame(dataset='GSE168740',library=colnames(g),route=NA,group=trt,time_or_age=age,vg_counts=as.numeric(unlist(g[VG,])),eff_lib_size=eff,
    vg_log2CPM=as.numeric(cpm(as.matrix(g[VG,,drop=FALSE]),lib.size=eff,log=TRUE,prior.count=2)),in_primary_model=TRUE,
    normalisation='TMM across all 18 libraries (display only; DE normalised within each age)')
  p=read.delim(file.path(DE,'PRJNA1017469_subsetA_lengthScaledTPM_counts.tsv'),row.names=1,check.names=FALSE)
  P=do.call(rbind,lapply(c('ingestion','contact'),function(rt){f=read.delim(file.path(DE,paste0('PRJNA1017469_sensitivity_allLibraries_',rt,'_TMM_factors.tsv')),colClasses=c(run='character'))
    e=f$lib.size*f$norm.factors; s=ss[match(f$run,ss$run),]
    data.frame(dataset='PRJNA1017469',library=f$isolate,route=rt,group=s$dose,time_or_age=s$time,vg_counts=as.numeric(unlist(p[VG,f$run])),eff_lib_size=e,
      vg_log2CPM=as.numeric(cpm(as.matrix(p[VG,f$run,drop=FALSE]),lib.size=e,log=TRUE,prior.count=2)),in_primary_model=s$decision=='RETAIN',
      normalisation='TMM factors of the all-libraries (sensitivity) model of this route, as written by the DE run')}))
  E=rbind(G,P); write.table(data.frame(gene_id=VG,gene_symbol='Vg',E),file.path(TD,'T8_Vg_per_library_expression.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  ## Fig 1a
  G$age=factor(sub('d',' d',G$time_or_age),levels=c('1 d','5 d','10 d')); G$group=factor(G$group,levels=c('Control','Dinotefuran'))
  a=VM[VM$dataset=='GSE168740'&VM$model=='per-age',]; a$age=factor(sub('d',' d',a$key),levels=levels(G$age))
  lab=do.call(rbind,lapply(levels(a$age),function(z){e=a[a$age==z&a$method=='edgeR QL',]; d=a[a$age==z&a$method=='DESeq2',]
    data.frame(age=factor(z,levels=levels(G$age)),label=sprintf('edgeR (primary): log2FC %.2f, FDR %s\nDESeq2 (secondary): log2FC %.2f, FDR %s',e$log2FC,fmtp(e$FDR),d$log2FC,fmtp(d$FDR)))}))
  yt=max(G$vg_log2CPM)+1.0
  f1=ggplot(G,aes(group,vg_log2CPM))+stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom='crossbar',width=0.45,linewidth=0.3,colour='#52514e')+
    geom_point(aes(fill=group),shape=21,size=2.4,stroke=0.3,colour='#2b2b29',position=position_jitter(width=0.12,height=0,seed=1))+
    geom_text(data=G[G$library=='C1_1d',],aes(label='C1_1d'),nudge_x=-0.3,size=2.2,colour='#52514e')+
    geom_text(data=lab,aes(x=1.5,y=yt,label=label),inherit.aes=FALSE,size=2.1,lineheight=0.95,colour='#0b0b0b')+
    facet_wrap(~age,nrow=1)+scale_fill_manual(values=c(Control='#c9c8c2',Dinotefuran='#2a78d6'),name=NULL)+
    scale_y_continuous(expand=expansion(mult=c(0.05,0.25)))+labs(x=NULL,y='Vg (GB49544) log2 CPM',
    title='Vitellogenin (GB49544) in GSE168740 libraries (tissue: brain, per GEO annotation)',
    subtitle=WR('Points = pooled libraries (15 brains per library, per the GEO record); bar = group mean. Labels: per-age DE results, primary edgeR QL and secondary-method DESeq2. C1_1d is the library dropped in the 1 d sensitivity analysis.',W_SUB),
    caption=WR('log2 CPM from length-scaled counts with TMM normalisation across all 18 libraries (prior count 2); display only, because each DE test normalised within its own age. FDR: Benjamini-Hochberg within each contrast. Preliminary 500,000-read subsets. RNA abundance only.',W_CAP))+
    TH+theme(legend.position='bottom')
  save_fig(f1,'Fig1a_Vg_expression_GSE168740',180,100)
  ## Fig 1b (revision 5: floor note for zero-count excluded libraries)
  P$dose=factor(P$group,levels=c('CTL','SUB','D50')); P$time=factor(sub('h',' h',P$time_or_age),levels=c('1 h','4 h')); P$route=factor(P$route,levels=c('ingestion','contact'))
  P$inc=factor(ifelse(P$in_primary_model,'included in primary model','excluded from primary model by Final QC'),levels=c('included in primary model','excluded from primary model by Final QC'))
  z0=P[!P$in_primary_model & P$vg_counts==0,]; zp=P[!P$in_primary_model & P$vg_counts>0,]
  b=VM[VM$dataset=='PRJNA1017469'&VM$tier=='Primary'&grepl('_vs_CTL_',VM$key),]
  b$dose=factor(sub('_.*','',b$key),levels=levels(P$dose)); b$time=factor(sub('h',' h',sub('.*_','',b$key)),levels=levels(P$time)); b$route=factor(b$route,levels=levels(P$route))
  b$label=sprintf('log2FC %.2f\nFDR %s',b$log2FC,fmtp(b$FDR)); yb=max(P$vg_log2CPM[P$in_primary_model])+1.3
  f2=ggplot(P,aes(dose,vg_log2CPM))+stat_summary(data=P[P$in_primary_model,],fun=mean,fun.min=mean,fun.max=mean,geom='crossbar',width=0.45,linewidth=0.3,colour='#52514e')+
    geom_point(aes(shape=inc,fill=dose),size=2.2,stroke=0.4,colour='#2b2b29',position=position_jitter(width=0.12,height=0,seed=2))+
    geom_text(data=b,aes(x=dose,y=yb,label=label),inherit.aes=FALSE,size=2.1,lineheight=0.95)+facet_grid(route~time)+
    scale_shape_manual(values=c(21,4),name=NULL,drop=FALSE)+scale_fill_manual(values=c(CTL='#c9c8c2',SUB='#86b6ef',D50='#1c5cab'),name='Dose group')+
    scale_y_continuous(expand=expansion(mult=c(0.05,0.2)))+labs(x=NULL,y='Vg (GB49544) log2 CPM',title='Vitellogenin (GB49544) in PRJNA1017469 head libraries',
    subtitle=WR('Points = pooled libraries (3 heads per library); bar = mean of the libraries in the primary model. Labels: primary edgeR QL contrast against CTL at the same time point; none reaches FDR < 0.05.',W_SUB),
    caption=WR(paste0('log2 CPM from length-scaled counts with the TMM factors of the all-libraries (sensitivity) model of each route (prior count 2). Crosses = libraries excluded from the primary model by Final QC; they enter only the sensitivity model. ',
      if(nrow(z0)) paste0('Crosses at the bottom of a panel (',paste(z0$library,collapse=', '),') have 0 Vg counts and are drawn at the prior-count floor of a very small library, not at a measured expression level. ') else '',
      if(nrow(zp)) paste0('Other excluded libraries: ',paste0(zp$library,' (',zp$vg_counts,' Vg counts)',collapse=', '),'. ') else '','Preliminary 500,000-read subsets. RNA abundance only.'),W_CAP))+
    TH+theme(legend.position='bottom',legend.box='vertical')+guides(fill=guide_legend(override.aes=list(shape=21)))
  save_fig(f2,'Fig1b_Vg_expression_PRJNA1017469',180,140)})
## =============== FIGURE 2 (revision 3: unique labels; revision 1: no clipping) ===============
step('Fig2 Vg forest (40 analyses)',{
  F=VM; F$sig=ifelse(is.na(F$FDR),'FDR not available',ifelse(F$FDR<0.05,'FDR < 0.05','FDR >= 0.05'))
  F=F[order(match(F$panel,PANELS),F$tier_rank,F$ord),]; F$lab=factor(F$label,levels=rev(unique(F$label))); F$panel=factor(F$panel,levels=PANELS)
  L=rbind(data.frame(F,measure='log2 fold change',value=F$log2FC),data.frame(F,measure='-log10 FDR',value=-log10(F$FDR)))
  L$measure=factor(L$measure,levels=c('log2 fold change','-log10 FDR')); L$fillk=ifelse(L$sig=='FDR < 0.05',L$tier,'open')
  ref=data.frame(measure=factor(c('log2 fold change','-log10 FDR'),levels=levels(L$measure)),x=c(0,-log10(0.05)),lt=c('solid','dashed'))
  na=L[L$measure=='-log10 FDR'&is.na(L$value),]
  f=ggplot(L[!is.na(L$value),],aes(value,lab))+geom_vline(data=ref,aes(xintercept=x,linetype=lt),linewidth=0.35,colour='#7d7c78',show.legend=FALSE)+
    geom_point(aes(colour=tier,fill=fillk,shape=method),size=2,stroke=0.6)+
    geom_text(data=na,aes(x=0.03,y=lab,label='NA (filtered)'),hjust=0,size=2,colour='#52514e',inherit.aes=FALSE)+
    facet_grid(panel~measure,scales='free',space='free_y',labeller=labeller(panel=label_wrap_gen(16)))+
    scale_colour_manual(values=COL_TIER,breaks=TIER_ORDER,name='Analysis tier')+scale_fill_manual(values=c(COL_TIER,open='white'),guide='none')+
    scale_shape_manual(values=c('edgeR QL'=21,'DESeq2'=22),name='Method')+scale_linetype_identity()+
    labs(x=NULL,y=NULL,title='Vitellogenin (GB49544): log2 fold change and FDR in all 40 analyses',
      subtitle=WR('Filled symbol = FDR < 0.05 (only GSE168740 10 d, edgeR QL, primary); open symbol = not significant. Dashed line = FDR 0.05. Route and tissue are given in the row panels.',W_SUB),
      caption=WR(paste('Values copied from DE run run_20261113_044101. Positive log2FC = higher in the first-named group (treated; later time; or ingestion for route differences).',
        'NA (filtered) = DESeq2 independent filtering removed the gene from multiple-testing adjustment (GSE168740 5 d). No confidence intervals are drawn because the edgeR QL tables do not report standard errors.',
        TIER_NOTE,'RNA level only; preliminary 500,000-read subsets.'),W_CAP))+
    TH+theme(legend.position='bottom',axis.text.y=element_text(size=6.2),strip.text.y=element_text(angle=0,hjust=0,size=6.6))
  save_fig(f,'Fig2_Vg_log2FC_all_analyses',180,240)})
## =============== FIGURE 3 (revision 6: r to 4 decimals) ===============
step('Fig3 edgeR vs DESeq2 (GSE168740)',{
  X=CMP; se=!is.na(X$FDR.edgeR)&X$FDR.edgeR<0.05; sd=!is.na(X$FDR.DESeq2)&X$FDR.DESeq2<0.05
  X$cat=factor(ifelse(se&sd,'FDR < 0.05 in both',ifelse(se,'edgeR only',ifelse(sd,'DESeq2 only','neither'))),levels=c('neither','edgeR only','DESeq2 only','FDR < 0.05 in both'))
  AL=c('1d'='1 d (primary vs secondary method)','5d'='5 d (primary vs secondary method)','10d'='10 d (primary vs secondary method)','1d_noC1_1d'='1 d without C1_1d (sensitivity, both methods)')
  X$age=factor(AL[X$age],levels=AL); X=X[order(X$cat),]
  r=do.call(rbind,lapply(levels(X$age),function(z){m=X[X$age==z,]; data.frame(age=factor(z,levels=AL),lab=sprintf('r = %.4f\nn = %d genes',cor(m$log2FC.edgeR,m$log2FC.DESeq2),nrow(m)))}))
  v=X[X$gene==VG,]; v$lab=sprintf('Vg: edgeR FDR %s\nDESeq2 FDR %s',fmtp(v$FDR.edgeR),fmtp(v$FDR.DESeq2))
  f=ggplot(X,aes(log2FC.edgeR,log2FC.DESeq2))+geom_abline(slope=1,intercept=0,linewidth=0.3,colour='#7d7c78',linetype='dashed')+
    geom_point(aes(colour=cat),size=0.7,alpha=0.85)+geom_point(data=v,shape=21,size=3,stroke=0.8,colour='#0b0b0b',fill=NA)+
    geom_text(data=v,aes(label=lab),hjust=-0.1,vjust=1.2,size=2.1,lineheight=0.95)+geom_text(data=r,aes(x=-Inf,y=Inf,label=lab),hjust=-0.08,vjust=1.2,size=2.2,inherit.aes=FALSE)+
    facet_wrap(~age,nrow=2)+scale_colour_manual(values=c(neither='#c9c8c2','edgeR only'='#2a78d6','DESeq2 only'='#eb6834','FDR < 0.05 in both'='#0b0b0b'),name=NULL)+
    coord_equal()+labs(x='log2 fold change, edgeR QL',y='log2 fold change, DESeq2',title='GSE168740: agreement between edgeR QL and DESeq2',
      subtitle=WR('Genes tested by both methods; circle = vitellogenin (GB49544). Colour = FDR < 0.05 status in each method. Dinotefuran vs control within each age.',W_SUB),
      caption=WR('edgeR QL is the primary method and DESeq2 the secondary method; the 1 d without C1_1d panel is a sensitivity analysis. Values copied from the DE run tables; r = Pearson correlation of log2FC over the genes shown. Effect sizes agree closely; the methods differ in which genes pass FDR < 0.05.',W_CAP))+
    TH+theme(legend.position='bottom')+guides(colour=guide_legend(override.aes=list(size=2)))
  save_fig(f,'Fig3_edgeR_vs_DESeq2_GSE168740',180,195)})
## =============== FIGURE 4 ===============
step('Fig4 DEG overview (40 analyses)',{
  D=DEGN; D=D[order(match(D$panel,PANELS),D$tier_rank,D$ord),]; D$lab=paste0(D$label,' (',D$genes_tested,' tested)')
  D$lab=factor(D$lab,levels=rev(unique(D$lab))); D$panel=factor(D$panel,levels=PANELS)
  B=rbind(data.frame(D,dir='Up',n=D$up),data.frame(D,dir='Down',n=-D$down)); B$tl=ifelse(B$n==0,'',abs(B$n))
  f=ggplot(B,aes(n,lab))+geom_vline(xintercept=0,linewidth=0.3,colour='#7d7c78')+geom_col(aes(fill=dir),width=0.7)+
    geom_text(aes(label=tl,hjust=ifelse(n>0,-0.25,1.25)),size=2)+geom_point(data=D,aes(x=0,y=lab,colour=tier),shape=15,size=1.5,inherit.aes=FALSE)+
    facet_grid(panel~.,scales='free_y',space='free_y',labeller=labeller(panel=label_wrap_gen(16)))+scale_fill_manual(values=COL_DIR[c('Up','Down')],name='Direction (FDR < 0.05)')+
    scale_colour_manual(values=COL_TIER,breaks=TIER_ORDER,name='Analysis tier')+scale_x_continuous(labels=abs,expand=expansion(mult=0.15))+
    labs(x='Number of genes with FDR < 0.05 (down | up)',y=NULL,title='Differentially expressed genes in all 40 analyses',
      subtitle=WR('Square at zero = analysis tier. Counts recomputed from the DE run tables; FDR is Benjamini-Hochberg within each contrast.',W_SUB),
      caption=WR(paste('Counts are not comparable across analyses with different numbers of tested genes, sample sizes or methods.',TIER_NOTE,'Preliminary 500,000-read subsets.'),W_CAP))+
    TH+theme(legend.position='bottom',legend.box='vertical',axis.text.y=element_text(size=6.1),strip.text.y=element_text(angle=0,hjust=0,size=6.6))
  save_fig(f,'Fig4_DEG_overview_all_analyses',180,235)})
## =============== FIGURE 5 / S1 (revision 2: readable headings, non-overlapping labels) ===============
place_labels=function(lab,Z){  # fallback when ggrepel is absent: per panel, stack labels so they never overlap
  out=NULL; for(fc in unique(as.character(lab$facet))){l=lab[as.character(lab$facet)==fc,]; z=Z[as.character(Z$facet)==fc,]
    l$ly=l$y; l$lx=l$log2FC; l$hj=0; yr=diff(range(z$y)); xr=diff(range(z$log2FC)); gap=0.075*yr; side=ifelse(l$log2FC>mean(range(z$log2FC)),-1,1)
    for(s in c(-1,1)){i=which(side==s); if(!length(i)) next; i=i[order(-l$y[i])]; ly=l$y[i]
      if(length(i)>1) for(j in 2:length(i)) ly[j]=min(ly[j],ly[j-1]-gap)
      l$ly[i]=ly; l$lx[i]=l$log2FC[i]+s*0.06*xr; l$hj[i]=ifelse(s>0,0,1)}
    out=rbind(out,l)}; out}
volc=function(sel,name,title,w,h,ncol){
  Z=do.call(rbind,lapply(sel,function(a){r=RES[[a]]; r$analysis=a; r}))
  Z=merge(Z,M[,c('analysis','comparison','panel','tier','dataset','route','key','ord')],by='analysis'); Z$sig=!is.na(Z$FDR)&Z$FDR<0.05
  Z$dir=factor(ifelse(Z$sig,ifelse(Z$log2FC>0,'Up','Down'),'Not significant'),levels=names(COL_DIR)); Z$cand=Z$gene %in% tg$ogs32_id
  hd=function(a){m=M[M$analysis==a,]; r=RES[[a]]
    l1=if(m$dataset=='GSE168740') paste0('GSE168740 (brain, GEO) | ',sub('d',' d',m$key)) else paste0('PRJNA1017469 head | ',m$route)
    l2=if(m$dataset=='GSE168740') paste0('Dinotefuran v control | ',m$tier) else paste0(sub(' \\(diagnostic\\)','',m$comparison),' | ',m$tier)
    paste0(l1,'\n',l2,'\n',sum(r$FDR<0.05,na.rm=TRUE),' of ',nrow(r),' genes FDR < 0.05')}
  sel=sel[order(M$ord[match(sel,M$analysis)])]; ttl=data.frame(analysis=sel,facet=sapply(sel,hd))
  Z$facet=factor(ttl$facet[match(Z$analysis,ttl$analysis)],levels=ttl$facet); Z$y=-log10(Z$P); Z=Z[order(Z$sig),]
  thr=do.call(rbind,lapply(split(Z,Z$facet,drop=TRUE),function(z) if(any(z$sig)) data.frame(facet=z$facet[1],y=-log10(max(z$P[z$sig]))) else NULL))
  v=Z[Z$gene==VG,]; v$lab=sprintf('Vg  log2FC %.2f\nFDR %s',v$log2FC,fmtp(v$FDR)); cs=Z[Z$cand & Z$sig & Z$gene!=VG,]
  cs$lab=tg$symbol[match(cs$gene,tg$ogs32_id)]; lab=rbind(v[,c('facet','log2FC','y','lab')],cs[,c('facet','log2FC','y','lab')])
  f=ggplot(Z,aes(log2FC,y))+geom_vline(xintercept=0,linewidth=0.25,colour='#b9b8b3')+geom_point(aes(colour=dir),size=0.6,alpha=0.9)+
    {if(!is.null(thr)) geom_hline(data=thr,aes(yintercept=y),linetype='dashed',linewidth=0.3,colour='#52514e')}+
    geom_point(data=Z[Z$cand,],shape=21,size=1.6,stroke=0.35,colour='#52514e',fill=NA)+
    geom_point(data=v,shape=21,size=2.6,stroke=0.8,colour='#0b0b0b',fill='#eda100')
  if(HAS_REPEL) f=f+ggrepel::geom_text_repel(data=lab,aes(label=lab),size=2,lineheight=0.9,min.segment.length=0,segment.size=0.2,segment.colour='#52514e',
       box.padding=0.35,point.padding=0.15,max.overlaps=Inf,seed=1,force=2)
  else {pl=place_labels(lab,Z); f=f+geom_segment(data=pl,aes(x=log2FC,y=y,xend=lx,yend=ly),linewidth=0.2,colour='#52514e')+
       geom_text(data=pl,aes(x=lx,y=ly,label=lab,hjust=hj),size=2,lineheight=0.9)}
  f=f+facet_wrap(~facet,ncol=ncol,scales='free')+scale_colour_manual(values=COL_DIR,name='FDR < 0.05',drop=FALSE)+
    scale_x_continuous(expand=expansion(mult=0.12))+scale_y_continuous(expand=expansion(mult=c(0.03,0.12)))+
    labs(x='log2 fold change',y='-log10 P',title=title,
      subtitle=WR('Yellow = vitellogenin (GB49544); ringed = curated project target genes with documented NCBI Gene annotation (labelled when FDR < 0.05). Dashed line = P value at the FDR 0.05 boundary (drawn only when at least one gene passes).',W_SUB),
      caption=WR(paste('Values copied from the DE run tables (edgeR quasi-likelihood F-tests). FDR is Benjamini-Hochberg within each contrast. Positive log2FC = higher in the first-named group.',TIER_NOTE,'Preliminary 500,000-read subsets.'),W_CAP))+
    TH+theme(legend.position='bottom',strip.text=element_text(size=6.3,lineheight=1))+guides(colour=guide_legend(override.aes=list(size=2)))
  save_fig(f,name,w,h)}
step('Fig5 volcano, selected primary contrasts',volc(c('GSE168740_subsetA_1d_edgeR_QLF','GSE168740_subsetA_5d_edgeR_QLF','GSE168740_subsetA_10d_edgeR_QLF',
  'PRJNA1017469_primary_ingestion_SUB_vs_CTL_1h_edgeR_QLF','PRJNA1017469_primary_ingestion_D50_vs_CTL_4h_edgeR_QLF','PRJNA1017469_primary_contact_D50_vs_CTL_4h_edgeR_QLF'),
  'Fig5_volcano_selected_primary_contrasts','Volcano plots of selected primary contrasts',180,160,3))
step('FigS1 volcano, all PRJNA1017469 primary',volc(M$analysis[M$dataset=='PRJNA1017469'&M$model=='primary'],'FigS1_volcano_PRJNA1017469_all_primary_contrasts',
  'PRJNA1017469 primary models: all 14 contrasts (CTL 4 h vs 1 h is diagnostic/QC)',180,265,3))
## =============== FIGURE 6 / S2 (revision 4: distinct not-tested marking, chronological columns, readable categories; revision 9: omitted genes) ===============
CAT_SHORT=c('storage (fat-body index)'='fat-body index','HPG marker (contamination index)'='HPG marker','NOT Vg (mapping correction)'='mapping control',
            'endocrine/behaviour'='endocrine/behav.','immune signalling'='immune signal.','yellow/MRJP family'='yellow/MRJP')
heat=function(an,name,title,w,h,colfacet){
  H=CAND[CAND$analysis %in% an,]; keep=unique(H$ogs32_id[H$tested]); omitted=setdiff(unique(CAND$symbol),unique(H$symbol[H$tested]))
  never=setdiff(unique(CAND$symbol),unique(CAND$symbol[CAND$tested])); only_else=setdiff(omitted,never)
  H=H[H$ogs32_id %in% keep,]; H$row=paste0(H$symbol,' (',H$ogs32_id,')')
  H$cat=ifelse(H$category %in% names(CAT_SHORT),CAT_SHORT[H$category],H$category)
  co=unique(H[,c('colname','tier_rank','dataset','ord')]); co=co[order(co$tier_rank,co$dataset,co$ord),]; H$col=factor(H$colname,levels=co$colname)
  H$cf=if(colfacet=='dataset') factor(H$dataset,levels=c('GSE168740','PRJNA1017469')) else factor(sub('Secondary method','Secondary\nmethod',sub('Diagnostic/QC','Diagnostic\n/QC',H$tier)),
       levels=sub('Secondary method','Secondary\nmethod',sub('Diagnostic/QC','Diagnostic\n/QC',TIER_ORDER)))
  H$val=ifelse(H$tested,pmax(pmin(H$log2FC,3),-3),NA); H$star=ifelse(!is.na(H$FDR)&H$FDR<0.05,'*','')
  f=ggplot(H,aes(col,row))+geom_tile(aes(fill=val),colour='white',linewidth=0.3)+
    geom_point(data=H[!H$tested,],aes(shape='not tested (removed by the expression filter)'),size=1.1,stroke=0.45,colour='#6f6e6a')+
    geom_text(aes(label=star),size=3,vjust=0.75)+
    facet_grid(cat~cf,scales='free',space='free')+
    scale_fill_gradient2(low='#2a78d6',mid='#f0efec',high='#e34948',midpoint=0,limits=c(-3,3),na.value='white',name='log2FC\n(capped at +/-3)')+
    scale_shape_manual(values=c('not tested (removed by the expression filter)'=4),name=NULL)+
    labs(x=NULL,y=NULL,title=title,
      subtitle=WR('Curated project target genes with documented NCBI Gene annotation; * = FDR < 0.05; white tile with a cross = not tested in that analysis (removed by the expression filter). Columns in chronological order within each group.',W_SUB),
      caption=WR(paste0('Row groups are the project curation labels (data/metadata/gene_targets_ncbi.tsv), not enrichment results; abbreviated here: fat-body index = storage (fat-body index); HPG marker = HPG marker (contamination index); mapping control = NOT Vg (mapping correction); endocrine/behav. = endocrine/behaviour; immune signal. = immune signalling. ',
        'Of the 64 curated genes, ',length(never),' were never tested in any analysis and are not shown: ',paste(sort(never),collapse=', '),'. ',
        if(length(only_else)) paste0('Tested only outside the analyses shown and therefore not shown: ',paste(sort(only_else),collapse=', '),'. ') else '',
        'Values copied from the DE run tables. ',TIER_NOTE),W_CAP))+
    TH+theme(axis.text.x=element_text(angle=60,hjust=1,size=6),axis.text.y=element_text(size=6),strip.text.y=element_text(angle=0,hjust=0,size=6,margin=margin(0.5,2,0.5,2)),
      strip.text.x=element_text(size=6.5),panel.grid=element_blank(),panel.spacing=unit(1.2,'mm'),legend.position='right')+
    guides(shape=guide_legend(order=2),fill=guide_colourbar(order=1))
  save_fig(f,name,w,h)
  say('  ',name,': ',length(keep),' genes shown; never tested ',length(never),'; tested only elsewhere ',length(only_else))}
step('Fig6 candidate heatmap, primary',heat(M$analysis[M$tier=='Primary'],'Fig6_candidate_genes_primary_analyses','Curated candidate genes in the primary analyses',180,255,'dataset'))
step('FigS2 candidate heatmap, all 40',heat(AN,'FigS2_candidate_genes_all_40_analyses','Curated candidate genes in all 40 analyses',270,265,'tier'))
## ---------------- finish ----------------
writeLines(capture.output(sessionInfo()),file.path(OUT,'R_sessionInfo.txt'))
write.table(data.frame(step=names(STATUS),status=STATUS),file.path(OUT,'STEP_STATUS.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
nerr=sum(STATUS!='OK'); say(sprintf('\nSUMMARY: %d of %d steps OK. Output: %s',sum(STATUS=='OK'),length(STATUS),normalizePath(OUT)))
if(nerr) {say('STEPS WITH ERRORS: ',paste(names(STATUS)[STATUS!='OK'],collapse=', ')); quit(status=1)}
writeLines(format(Sys.time(),'%Y-%m-%dT%H:%M:%S'),file.path(OUT,'_FINISHED'))
