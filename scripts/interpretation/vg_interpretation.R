## Vg-centred interpretation tables and publication figures (2026-11-13).
## Reads ONLY the validated DE outputs in results/rawreads/DE/run_20261113_044101/ plus the curated target-gene annotation
## data/metadata/gene_targets_ncbi.tsv (NCBI Gene descriptions). Nothing is simulated, re-fitted or re-tested:
## every number plotted is copied from an edgeR/DESeq2 output table, or (Fig. 1 only) is a per-library CPM computed from the
## run's own length-scaled count matrix and the run's own TMM factors.
## Writes a NEW folder results/rawreads/Vg_interpretation/run_<time>/ (tables/, figures/ as PDF + 300-dpi LZW TIFF, logs).
## Run from the project folder:  Rscript --vanilla scripts/interpretation/vg_interpretation.R   (or double-click RUN_VG_FIGURES.bat)
rl=file.path('envs','R-library',paste(R.version$major,sub('\\..*','',R.version$minor),sep='.')); if(dir.exists(rl)) .libPaths(c(normalizePath(rl),.libPaths()))
for(p in c('ggplot2','edgeR')) if(!requireNamespace(p,quietly=TRUE)) stop('required package not installed: ',p,' (nothing is installed by this script)')
suppressMessages({library(ggplot2);library(edgeR)})
DE='results/rawreads/DE/run_20261113_044101'; if(!file.exists(file.path(DE,'_DE_FINISHED'))) stop('validated DE run not found: ',DE)
OUT=file.path('results/rawreads/Vg_interpretation',format(Sys.time(),'run_%Y%m%d_%H%M%S')); TD=file.path(OUT,'tables'); FD=file.path(OUT,'figures')
dir.create(TD,recursive=TRUE); dir.create(FD,recursive=TRUE)
LOG=file.path(OUT,'RUN_LOG.txt'); say=function(...){m=paste0(...); cat(m,'\n'); cat(m,'\n',file=LOG,append=TRUE)}
STATUS=character(0); step=function(name,expr){r=tryCatch({force(expr);'OK'},error=function(e) paste('ERROR:',conditionMessage(e))); STATUS[name]<<-r; say(sprintf('%-45s %s',name,r))}
say('vg_interpretation.R  ',format(Sys.time()),'  ',R.version.string,'  ggplot2 ',as.character(packageVersion('ggplot2')),'  edgeR ',as.character(packageVersion('edgeR')))
say('DE input folder: ',DE)
file.copy('scripts/interpretation/vg_interpretation.R',file.path(OUT,'vg_interpretation.R.used'))
VG='GB49544'
## ---------------- consistent styling ----------------
COL_TIER=c('Primary'='#2a78d6','Sensitivity'='#eb6834','Secondary method'='#e87ba4','Exploratory'='#1baf7a','Diagnostic/QC'='#7d7c78')
COL_DIR=c('Up'='#e34948','Down'='#2a78d6','Not significant'='#c9c8c2')
TH=theme_bw(base_size=8,base_family='sans')+theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(linewidth=0.25,colour='#e6e5e1'),
   strip.background=element_rect(fill='#f0efec',colour='#b9b8b3',linewidth=0.3),strip.text=element_text(size=7.5,face='bold'),
   panel.border=element_rect(colour='#8f8e89',linewidth=0.4),axis.ticks=element_line(linewidth=0.3),legend.key.size=unit(3.5,'mm'),
   legend.title=element_text(size=7.5,face='bold'),legend.text=element_text(size=7),plot.title=element_text(size=9,face='bold'),
   plot.subtitle=element_text(size=7.5,colour='#52514e'),plot.caption=element_text(size=6.5,colour='#52514e',hjust=0))
TIFTYPE=if(isTRUE(capabilities('cairo'))) 'cairo' else NULL
save_fig=function(p,name,w,h){  # w, h in mm; PDF (editable text, no dingbats) + 300-dpi LZW TIFF
  ggsave(file.path(FD,paste0(name,'.pdf')),p,width=w,height=h,units='mm',device='pdf',useDingbats=FALSE)
  if(is.null(TIFTYPE)) ggsave(file.path(FD,paste0(name,'.tiff')),p,width=w,height=h,units='mm',dpi=300,device='tiff',compression='lzw')
  else ggsave(file.path(FD,paste0(name,'.tiff')),p,width=w,height=h,units='mm',dpi=300,device='tiff',compression='lzw',type=TIFTYPE)}
fmtp=function(x) ifelse(is.na(x),'NA',ifelse(x<0.001,formatC(x,format='e',digits=1),formatC(x,format='f',digits=3)))
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
## standardise one result table
std=function(a){t=read.delim(file.path(DE,paste0(a,'.tsv')),row.names=1,check.names=FALSE)
  if('logFC' %in% names(t)) data.frame(gene=rownames(t),log2FC=t$logFC,expr=t$logCPM,expr_measure='logCPM',SE=NA_real_,P=t$PValue,FDR=t$FDR)
  else data.frame(gene=rownames(t),log2FC=t$log2FoldChange,expr=t$baseMean,expr_measure='baseMean',SE=t$lfcSE,P=t$pvalue,FDR=t$padj)}
RES=setNames(lapply(AN,std),AN)
## analysis metadata (parsed from the file names written by de_from_kallisto.R)
CLAB=c(SUB_vs_CTL_1h='SUB vs CTL, 1 h',D50_vs_CTL_1h='D50 vs CTL, 1 h',SUB_vs_CTL_4h='SUB vs CTL, 4 h',D50_vs_CTL_4h='D50 vs CTL, 4 h',
       SUBxTime_interaction='SUB x time interaction',D50xTime_interaction='D50 x time interaction',CTL_4h_vs_1h_drift='CTL 4 h vs 1 h (diagnostic)')
CGRP=list(SUB_vs_CTL_1h=c('SUB_1h','CTL_1h'),D50_vs_CTL_1h=c('D50_1h','CTL_1h'),SUB_vs_CTL_4h=c('SUB_4h','CTL_4h'),D50_vs_CTL_4h=c('D50_4h','CTL_4h'),
       SUBxTime_interaction=c('SUB_4h','CTL_4h','SUB_1h','CTL_1h'),D50xTime_interaction=c('D50_4h','CTL_4h','D50_1h','CTL_1h'),CTL_4h_vs_1h_drift=c('CTL_4h','CTL_1h'))
ngrp=function(route,set,groups){s=ss[ss$route %in% route & (set=='all' | ss$decision=='RETAIN'),]; paste(sapply(groups,function(g) paste(g,sum(s$group==g))),collapse=' v ')}
meta=function(a){
  if(startsWith(a,'GSE168740')){age=regmatches(a,regexpr('(10d|5d|1d)',a)); noC1=grepl('noC1_1d',a); d2=grepl('DESeq2$',a)
    return(data.frame(analysis=a,dataset='GSE168740',tissue='brain (pools of 15)',panel='GSE168740 (brain)',route=NA,model=ifelse(noC1,'sensitivity (without C1_1d)','per-age'),
      method=ifelse(d2,'DESeq2','edgeR QL'),tier=ifelse(noC1,'Sensitivity',ifelse(d2,'Secondary method','Primary')),
      comparison=paste0('Dinotefuran vs control, ',sub('d',' d',age),ifelse(noC1,' (without C1_1d)','')),key=age,n_groups=ifelse(noC1,'control 2 v treated 3','control 3 v treated 3'),
      short=paste0(sub('d',' d',age),ifelse(noC1,' -C1',''),' ',ifelse(d2,'DESeq2','edgeR'))))}
  ex=grepl('_EXPLORATORY_',a); mdl=ifelse(ex,'exploratory',ifelse(grepl('_primary_',a),'primary','sensitivity (all 47 libraries)'))
  if(ex){k=sub('.*routeDifference_(SUB|D50)_(1h|4h)_edgeR_QLF$','\\1_\\2',a); d=sub('_.*','',k); t=sub('.*_','',k)
    return(data.frame(analysis=a,dataset='PRJNA1017469',tissue='head (pools of 3)',panel='PRJNA1017469 route difference',route='ingestion - contact',model=mdl,method='edgeR QL',tier='Exploratory',
      comparison=paste0('Route difference of the ',d,' effect at ',sub('h',' h',t),' (ingestion - contact)'),key=k,
      n_groups=paste0('ingestion ',ngrp('ingestion','primary',c(paste0(d,'_',t),paste0('CTL_',t))),'; contact ',ngrp('contact','primary',c(paste0(d,'_',t),paste0('CTL_',t)))),
      short=paste0('Route diff ',d,' ',t)))}
  rt=ifelse(grepl('_ingestion_',a),'ingestion','contact'); k=sub(paste0('.*_',rt,'_(.*)_edgeR_QLF$'),'\\1',a)
  data.frame(analysis=a,dataset='PRJNA1017469',tissue='head (pools of 3)',panel=paste('PRJNA1017469',rt),route=rt,model=mdl,method='edgeR QL',
    tier=ifelse(k=='CTL_4h_vs_1h_drift','Diagnostic/QC',ifelse(mdl=='primary','Primary','Sensitivity')),comparison=CLAB[[k]],key=k,
    n_groups=ngrp(rt,ifelse(mdl=='primary','primary','all'),CGRP[[k]]),short=paste0(ifelse(rt=='ingestion','Ing ','Con '),sub(' \\(diagnostic\\)','',CLAB[[k]]),ifelse(mdl=='primary','',' [sens]')))}
M=do.call(rbind,lapply(AN,meta))
## =============== TABLE 1: Vg across all 40 analyses ===============
step('T1 Vg consolidated table',{
  V=do.call(rbind,lapply(AN,function(a){r=RES[[a]]; o=order(r$P); i=match(VG,r$gene)
    data.frame(analysis=a,tested=!is.na(i),log2FC=r$log2FC[i],expression=r$expr[i],expression_measure=r$expr_measure[1],lfcSE_DESeq2=r$SE[i],P=r$P[i],FDR=r$FDR[i],
      rank_by_P=if(is.na(i)) NA else match(VG,r$gene[o]),genes_tested=nrow(r),DEG_FDR05=sum(r$FDR<0.05,na.rm=TRUE))}))
  V=merge(M,V,by='analysis',sort=FALSE)
  V$direction=ifelse(is.na(V$log2FC),'not available',ifelse(V$log2FC>0,'higher in first group (up)','lower in first group (down)'))
  V$significance=ifelse(is.na(V$FDR),'FDR not available (DESeq2 independent filtering)',ifelse(V$FDR<0.05,'significant (FDR<0.05)',ifelse(V$P<0.05,'not significant (nominal P<0.05 only; FDR>=0.05)','not significant')))
  chk=merge(V[,c('analysis','log2FC','P','FDR')],vgt[,c('analysis','logFC','PValue','FDR')],by='analysis',suffixes=c('','.vgtable'))
  if(any(abs(chk$log2FC-chk$logFC)>1e-9) || any(abs(chk$P-chk$PValue)>1e-12)) stop('Vg values disagree with Vg_GB49544_all_analyses.tsv')
  ord=order(match(V$panel,c('GSE168740 (brain)','PRJNA1017469 ingestion','PRJNA1017469 contact','PRJNA1017469 route difference')),match(V$tier,names(COL_TIER)),V$analysis)
  V<<-V[ord,]; write.table(V,file.path(TD,'T1_Vg_GB49544_all_40_analyses.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  say('  Vg significant (FDR<0.05) in: ',paste(V$analysis[which(V$FDR<0.05)],collapse=', '))})
## =============== TABLE 2/3: edgeR vs DESeq2 (GSE168740) ===============
step('T2 edgeR vs DESeq2 concordance',{
  K=c('1d','5d','10d','1d_noC1_1d'); CMP<<-do.call(rbind,lapply(K,function(k){e=RES[[paste0('GSE168740_subsetA_',k,'_edgeR_QLF')]]; d=RES[[paste0('GSE168740_subsetA_',k,'_DESeq2')]]
    m=merge(e,d,by='gene',suffixes=c('.edgeR','.DESeq2')); m$age=k; m}))
  S=do.call(rbind,lapply(split(CMP,CMP$age),function(m){se=!is.na(m$FDR.edgeR)&m$FDR.edgeR<0.05; sd=!is.na(m$FDR.DESeq2)&m$FDR.DESeq2<0.05; b=se&sd
    data.frame(comparison=m$age[1],genes_edgeR=nrow(RES[[paste0('GSE168740_subsetA_',m$age[1],'_edgeR_QLF')]]),genes_DESeq2=nrow(RES[[paste0('GSE168740_subsetA_',m$age[1],'_DESeq2')]]),
      genes_in_both=nrow(m),pearson_r_log2FC=round(cor(m$log2FC.edgeR,m$log2FC.DESeq2),3),spearman_rho_log2FC=round(cor(m$log2FC.edgeR,m$log2FC.DESeq2,method='spearman'),3),
      DEG_edgeR=sum(se),DEG_DESeq2=sum(sd),DEG_both=sum(b),DEG_edgeR_only=sum(se&!sd),DEG_DESeq2_only=sum(sd&!se),
      same_sign_among_both=sum(sign(m$log2FC.edgeR[b])==sign(m$log2FC.DESeq2[b])),DESeq2_FDR_NA_among_edgeR_DEG=sum(se&is.na(m$FDR.DESeq2)))}))
  write.table(S,file.path(TD,'T2_edgeR_vs_DESeq2_genomewide_GSE168740.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  v=CMP[CMP$gene==VG,c('age','log2FC.edgeR','log2FC.DESeq2','P.edgeR','P.DESeq2','FDR.edgeR','FDR.DESeq2')]
  v$agreement=ifelse(is.na(v$FDR.DESeq2),'DESeq2 FDR not available',ifelse((v$FDR.edgeR<0.05)==(v$FDR.DESeq2<0.05),
     ifelse(v$FDR.edgeR<0.05,'both significant','both not significant'),'DISAGREE on significance'))
  v$same_direction=sign(v$log2FC.edgeR)==sign(v$log2FC.DESeq2)
  write.table(v,file.path(TD,'T3_Vg_edgeR_vs_DESeq2_GSE168740.tsv'),sep='\t',quote=FALSE,row.names=FALSE)})
## =============== TABLE 4/5: DEG counts and DEG lists with documented annotation ===============
step('T4-T5 DEG counts and lists',{
  D=do.call(rbind,lapply(AN,function(a){r=RES[[a]]; s=!is.na(r$FDR)&r$FDR<0.05; data.frame(analysis=a,genes_tested=nrow(r),DEG_FDR05=sum(s),up=sum(s&r$log2FC>0),down=sum(s&r$log2FC<0),FDR_NA=sum(is.na(r$FDR)))}))
  DEGN<<-merge(M,D,by='analysis'); write.table(DEGN,file.path(TD,'T4_DEG_counts_all_40_analyses.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  L=do.call(rbind,lapply(AN,function(a){r=RES[[a]]; r=r[!is.na(r$FDR)&r$FDR<0.05,]; if(!nrow(r)) return(NULL); r$analysis=a; r}))
  i=match(L$gene,tg$ogs32_id); L$documented_symbol=tg$symbol[i]; L$NCBI_gene_id=tg$ncbi_gene_id[i]; L$NCBI_description=tg$ncbi_description[i]; L$project_category=tg$category[i]
  L$annotation_status=ifelse(is.na(i),'no annotation in project resources (OGS v3.2 ID only)','documented: NCBI Gene record (curated project target list)')
  L=merge(M[,c('analysis','dataset','panel','tier','method','comparison')],L,by='analysis'); L=L[order(L$analysis,L$P),]
  write.table(L,file.path(TD,'T5_DEG_FDR05_all_analyses_with_documented_annotation.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  say('  DEG rows: ',nrow(L),'; unique genes: ',length(unique(L$gene)),'; with documented annotation: ',length(unique(L$gene[!is.na(L$documented_symbol)])))})
## =============== TABLE 6/7: curated candidate panel (documented annotations only) ===============
step('T6-T7 candidate panel',{
  C=do.call(rbind,lapply(AN,function(a){r=RES[[a]]; i=match(tg$ogs32_id,r$gene)
    data.frame(analysis=a,symbol=tg$symbol,ogs32_id=tg$ogs32_id,NCBI_gene_id=tg$ncbi_gene_id,NCBI_description=tg$ncbi_description,category=tg$category,
      tested=!is.na(i),log2FC=r$log2FC[i],expression=r$expr[i],P=r$P[i],FDR=r$FDR[i])}))
  CAND<<-merge(M[,c('analysis','dataset','panel','tier','method','comparison','short')],C,by='analysis')
  write.table(CAND,file.path(TD,'T6_candidate_panel_all_analyses.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  h=CAND[!is.na(CAND$FDR)&CAND$FDR<0.05,]; write.table(h[order(h$category,h$symbol),],file.path(TD,'T7_candidate_panel_FDR05_hits.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  say('  candidate genes with FDR<0.05 in any analysis: ',nrow(h),' rows (',paste(unique(paste0(h$symbol,' [',h$short,']')),collapse='; '),')')})
## =============== TABLE 8 + FIGURE 1: per-library Vg expression ===============
step('T8 + Fig1 Vg expression per library',{
  g=read.delim(file.path(DE,'GSE168740_subsetA_lengthScaledTPM_counts.tsv'),row.names=1,check.names=FALSE)
  age=sub('.*_','',colnames(g)); trt=ifelse(substr(colnames(g),1,1)=='C','Control','Dinotefuran'); y=DGEList(as.matrix(g))
  k=filterByExpr(y,group=paste(age,trt)); y2=calcNormFactors(y[k,,keep.lib.sizes=FALSE]); eff=y2$samples$lib.size*y2$samples$norm.factors
  G=data.frame(dataset='GSE168740',library=colnames(g),route=NA,group=trt,time_or_age=age,vg_counts=as.numeric(unlist(g[VG,])),eff_lib_size=eff,
    vg_log2CPM=as.numeric(cpm(as.matrix(g[VG,,drop=FALSE]),lib.size=eff,log=TRUE,prior.count=2)),in_primary_model=TRUE,
    normalisation='TMM across all 18 libraries (display only; DE normalised within each age)')
  p=read.delim(file.path(DE,'PRJNA1017469_subsetA_lengthScaledTPM_counts.tsv'),row.names=1,check.names=FALSE)
  P=do.call(rbind,lapply(c('ingestion','contact'),function(rt){f=read.delim(file.path(DE,paste0('PRJNA1017469_sensitivity_allLibraries_',rt,'_TMM_factors.tsv')),colClasses=c(run='character'))
    e=f$lib.size*f$norm.factors; s=ss[match(f$run,ss$run),]
    data.frame(dataset='PRJNA1017469',library=f$isolate,route=rt,group=s$dose,time_or_age=s$time,vg_counts=as.numeric(unlist(p[VG,f$run])),eff_lib_size=e,
      vg_log2CPM=as.numeric(cpm(as.matrix(p[VG,f$run,drop=FALSE]),lib.size=e,log=TRUE,prior.count=2)),in_primary_model=s$decision=='RETAIN',
      normalisation='TMM factors of the all-libraries (sensitivity) model of this route, as written by the DE run')}))
  E=rbind(G,P); write.table(E,file.path(TD,'T8_Vg_per_library_expression.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  ## Fig 1a GSE168740
  G$age=factor(sub('d',' d',G$time_or_age),levels=c('1 d','5 d','10 d')); G$group=factor(G$group,levels=c('Control','Dinotefuran'))
  a=V[V$dataset=='GSE168740'&V$model=='per-age',]; a$age=factor(sub('d',' d',a$key),levels=levels(G$age))
  lab=do.call(rbind,lapply(split(a,a$age),function(z){if(!nrow(z)) return(NULL); e=z[z$method=='edgeR QL',]; d=z[z$method=='DESeq2',]
    data.frame(age=z$age[1],label=sprintf('edgeR: log2FC %.2f, FDR %s\nDESeq2: log2FC %.2f, FDR %s',e$log2FC,fmtp(e$FDR),d$log2FC,fmtp(d$FDR)))}))
  yt=max(G$vg_log2CPM)+0.9
  f1=ggplot(G,aes(group,vg_log2CPM))+stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom='crossbar',width=0.45,linewidth=0.3,colour='#52514e')+
    geom_point(aes(fill=group),shape=21,size=2.4,stroke=0.3,colour='#2b2b29',position=position_jitter(width=0.12,height=0,seed=1))+
    geom_text(data=G[G$library=='C1_1d',],aes(label='C1_1d'),nudge_x=-0.28,size=2.3,colour='#52514e')+
    geom_text(data=lab,aes(x=1.5,y=yt,label=label),inherit.aes=FALSE,size=2.25,lineheight=0.95,colour='#0b0b0b')+
    facet_wrap(~age,nrow=1)+scale_fill_manual(values=c(Control='#c9c8c2',Dinotefuran='#2a78d6'),name=NULL)+
    scale_y_continuous(expand=expansion(mult=c(0.05,0.22)))+labs(x=NULL,y='Vg (GB49544) log2 CPM',
    title='Vitellogenin (GB49544) in GSE168740 brain libraries',subtitle='Points = pooled libraries (15 brains each); bar = group mean. Labels: per-age DE results (primary edgeR QL; secondary DESeq2).',
    caption='log2 CPM from length-scaled counts with TMM normalisation across all 18 libraries (prior count 2); display only, the DE tests normalised within each age.\nPreliminary 500,000-read subsets. FDR: Benjamini-Hochberg within each contrast.')+TH+theme(legend.position='bottom')
  save_fig(f1,'Fig1a_Vg_expression_GSE168740',180,95)
  ## Fig 1b PRJNA1017469
  P$dose=factor(P$group,levels=c('CTL','SUB','D50')); P$time=factor(sub('h',' h',P$time_or_age),levels=c('1 h','4 h')); P$route=factor(P$route,levels=c('ingestion','contact'))
  P$inc=factor(ifelse(P$in_primary_model,'included in primary model','excluded from primary model (QC)'),levels=c('included in primary model','excluded from primary model (QC)'))
  b=V[V$dataset=='PRJNA1017469'&V$tier=='Primary'&grepl('_vs_CTL_',V$key),]
  b$dose=factor(sub('_.*','',b$key),levels=levels(P$dose)); b$time=factor(sub('h',' h',sub('.*_','',b$key)),levels=levels(P$time)); b$route=factor(b$route,levels=levels(P$route))
  b$label=sprintf('log2FC %.2f\nFDR %s',b$log2FC,fmtp(b$FDR)); yb=max(P$vg_log2CPM[P$in_primary_model])+1.2
  f2=ggplot(P,aes(dose,vg_log2CPM))+stat_summary(data=P[P$in_primary_model,],fun=mean,fun.min=mean,fun.max=mean,geom='crossbar',width=0.45,linewidth=0.3,colour='#52514e')+
    geom_point(aes(shape=inc,fill=dose),size=2.2,stroke=0.35,colour='#2b2b29',position=position_jitter(width=0.12,height=0,seed=2))+
    geom_text(data=b,aes(x=dose,y=yb,label=label),inherit.aes=FALSE,size=2.2,lineheight=0.95)+facet_grid(route~time)+
    scale_shape_manual(values=c(21,4),name=NULL,drop=FALSE)+scale_fill_manual(values=c(CTL='#c9c8c2',SUB='#86b6ef',D50='#1c5cab'),name='Dose group')+
    scale_y_continuous(expand=expansion(mult=c(0.05,0.18)))+labs(x=NULL,y='Vg (GB49544) log2 CPM',title='Vitellogenin (GB49544) in PRJNA1017469 head libraries',
    subtitle='Points = pooled libraries (3 heads each); bar = mean of libraries in the primary model. Labels: primary edgeR QL contrast vs CTL at the same time.',
    caption='log2 CPM from length-scaled counts with the TMM factors of the all-libraries model of each route (prior count 2). Crosses = libraries excluded from the\nprimary model by Final QC (shown for transparency; in the sensitivity model only). Preliminary 500,000-read subsets. No result reaches FDR < 0.05.')+
    TH+theme(legend.position='bottom',legend.box='vertical')+guides(fill=guide_legend(override.aes=list(shape=21)))
  save_fig(f2,'Fig1b_Vg_expression_PRJNA1017469',180,130)})
## =============== FIGURE 2: Vg log2FC across all 40 analyses ===============
step('Fig2 Vg forest (40 analyses)',{
  F=V; F$sig=ifelse(is.na(F$FDR),'FDR not available',ifelse(F$FDR<0.05,'FDR < 0.05','FDR >= 0.05'))
  F$lab=paste0(F$comparison,' [',F$method,ifelse(F$tier %in% c('Primary'),'',paste0('; ',tolower(F$tier))),']')
  F$lab=factor(F$lab,levels=rev(unique(F$lab))); F$panel=factor(F$panel,levels=c('GSE168740 (brain)','PRJNA1017469 ingestion','PRJNA1017469 contact','PRJNA1017469 route difference'))
  L=rbind(data.frame(F,measure='log2 fold change',value=F$log2FC),data.frame(F,measure='-log10 FDR',value=-log10(F$FDR)))
  L$measure=factor(L$measure,levels=c('log2 fold change','-log10 FDR')); L$fillk=ifelse(L$sig=='FDR < 0.05',L$tier,'open')
  ref=data.frame(measure=factor(c('log2 fold change','-log10 FDR'),levels=levels(L$measure)),x=c(0,-log10(0.05)),lt=c('solid','dashed'))
  na=L[L$measure=='-log10 FDR'&is.na(L$value),]
  f=ggplot(L[!is.na(L$value),],aes(value,lab))+geom_vline(data=ref,aes(xintercept=x,linetype=lt),linewidth=0.35,colour='#7d7c78',show.legend=FALSE)+
    geom_point(aes(colour=tier,fill=fillk,shape=method),size=2.1,stroke=0.6)+
    geom_text(data=na,aes(x=0.05,y=lab,label='NA (DESeq2 independent filtering)'),hjust=0,size=2.1,colour='#52514e',inherit.aes=FALSE)+
    facet_grid(panel~measure,scales='free',space='free_y',labeller=label_wrap_gen(18))+
    scale_colour_manual(values=COL_TIER,name='Analysis tier')+scale_fill_manual(values=c(COL_TIER,open='white'),guide='none')+
    scale_shape_manual(values=c('edgeR QL'=21,'DESeq2'=22),name='Method')+scale_linetype_identity()+
    labs(x=NULL,y=NULL,title='Vitellogenin (GB49544): estimated log2 fold change and FDR in all 40 analyses',
      subtitle='Filled symbol = FDR < 0.05 (only GSE168740 10 d, edgeR QL); open = not significant. Dashed line = FDR 0.05.',
      caption='Values copied from the DE run run_20261113_044101. Positive log2FC = higher in the first-named group (treated, later time, or ingestion for route differences).\nNo confidence intervals are drawn because the edgeR QL tables do not report standard errors. RNA level only; preliminary 500,000-read subsets.')+
    TH+theme(legend.position='bottom',axis.text.y=element_text(size=6.3),strip.text.y=element_text(angle=0,hjust=0))
  save_fig(f,'Fig2_Vg_log2FC_all_analyses',180,230)})
## =============== FIGURE 3: edgeR vs DESeq2 ===============
step('Fig3 edgeR vs DESeq2 (GSE168740)',{
  X=CMP; se=!is.na(X$FDR.edgeR)&X$FDR.edgeR<0.05; sd=!is.na(X$FDR.DESeq2)&X$FDR.DESeq2<0.05
  X$cat=factor(ifelse(se&sd,'FDR < 0.05 in both',ifelse(se,'edgeR only',ifelse(sd,'DESeq2 only','neither'))),levels=c('neither','edgeR only','DESeq2 only','FDR < 0.05 in both'))
  X$age=factor(c('1d'='1 d','5d'='5 d','10d'='10 d','1d_noC1_1d'='1 d without C1_1d (sensitivity)')[X$age],levels=c('1 d','5 d','10 d','1 d without C1_1d (sensitivity)'))
  X=X[order(X$cat),]; r=do.call(rbind,lapply(split(X,X$age),function(m) data.frame(age=m$age[1],lab=sprintf('r = %.2f\nn = %d genes',cor(m$log2FC.edgeR,m$log2FC.DESeq2),nrow(m)))))
  v=X[X$gene==VG,]; v$lab=sprintf('Vg  edgeR FDR %s\n      DESeq2 FDR %s',fmtp(v$FDR.edgeR),fmtp(v$FDR.DESeq2))
  f=ggplot(X,aes(log2FC.edgeR,log2FC.DESeq2))+geom_abline(slope=1,intercept=0,linewidth=0.3,colour='#7d7c78',linetype='dashed')+
    geom_point(aes(colour=cat),size=0.7,alpha=0.85)+geom_point(data=v,shape=21,size=3,stroke=0.8,colour='#0b0b0b',fill=NA)+
    geom_text(data=v,aes(label=lab),hjust=-0.12,vjust=1.1,size=2.1)+geom_text(data=r,aes(x=-Inf,y=Inf,label=lab),hjust=-0.1,vjust=1.2,size=2.3,inherit.aes=FALSE)+
    facet_wrap(~age,nrow=2)+scale_colour_manual(values=c(neither='#c9c8c2','edgeR only'='#2a78d6','DESeq2 only'='#eb6834','FDR < 0.05 in both'='#0b0b0b'),name=NULL)+
    coord_equal()+labs(x='log2 fold change, edgeR QL (primary)',y='log2 fold change, DESeq2 (secondary)',title='GSE168740: agreement between edgeR QL and DESeq2',
      subtitle='Genes tested by both methods; circle = vitellogenin (GB49544).',caption='Dinotefuran vs control within each age. Values copied from the DE run tables; r = Pearson correlation of log2FC.')+
    TH+theme(legend.position='bottom')+guides(colour=guide_legend(override.aes=list(size=2)))
  save_fig(f,'Fig3_edgeR_vs_DESeq2_GSE168740',180,190)})
## =============== FIGURE 4: DEG overview ===============
step('Fig4 DEG overview (40 analyses)',{
  D=DEGN; D$lab=paste0(D$short,' (',D$genes_tested,' tested)'); D$panel=factor(D$panel,levels=c('GSE168740 (brain)','PRJNA1017469 ingestion','PRJNA1017469 contact','PRJNA1017469 route difference'))
  D=D[order(D$panel,match(D$tier,names(COL_TIER)),D$analysis),]; D$lab=factor(D$lab,levels=rev(unique(D$lab)))
  B=rbind(data.frame(D,dir='Up',n=D$up),data.frame(D,dir='Down',n=-D$down)); B$tl=ifelse(B$n==0,'',abs(B$n))
  f=ggplot(B,aes(n,lab))+geom_vline(xintercept=0,linewidth=0.3,colour='#7d7c78')+geom_col(aes(fill=dir),width=0.7)+
    geom_text(aes(label=tl,hjust=ifelse(n>0,-0.2,1.2)),size=2.1)+geom_point(data=D,aes(x=0,y=lab,colour=tier),shape=15,size=1.6,inherit.aes=FALSE)+
    facet_grid(panel~.,scales='free_y',space='free_y',labeller=label_wrap_gen(18))+scale_fill_manual(values=COL_DIR[c('Up','Down')],name='Direction (FDR < 0.05)')+
    scale_colour_manual(values=COL_TIER,name='Analysis tier')+scale_x_continuous(labels=abs,expand=expansion(mult=0.15))+
    labs(x='Number of genes with FDR < 0.05 (down | up)',y=NULL,title='Differentially expressed genes in all 40 analyses',
      subtitle='Square at zero = analysis tier. Counts recomputed from the DE run tables; FDR is Benjamini-Hochberg within each contrast.',
      caption='Gene counts are not comparable across analyses with different numbers of tested genes, sample sizes or methods. Preliminary 500,000-read subsets.')+
    TH+theme(legend.position='bottom',legend.box='vertical',axis.text.y=element_text(size=6.3),strip.text.y=element_text(angle=0,hjust=0))
  save_fig(f,'Fig4_DEG_overview_all_analyses',180,220)})
## =============== FIGURE 5 (+S1): volcano plots ===============
volc=function(sel,name,title,w,h,ncol){
  Z=do.call(rbind,lapply(sel,function(a){r=RES[[a]]; r$analysis=a; r}))
  Z=merge(Z,M[,c('analysis','comparison','panel','tier')],by='analysis'); Z$sig=!is.na(Z$FDR)&Z$FDR<0.05
  Z$dir=factor(ifelse(Z$sig,ifelse(Z$log2FC>0,'Up','Down'),'Not significant'),levels=names(COL_DIR)); Z$cand=Z$gene %in% tg$ogs32_id
  ttl=do.call(rbind,lapply(sel,function(a){r=RES[[a]]; m=M[M$analysis==a,]; data.frame(analysis=a,facet=sprintf('%s\n%s | %s: %d of %d genes FDR < 0.05',m$panel,m$comparison,m$tier,sum(r$FDR<0.05,na.rm=TRUE),nrow(r)))}))
  Z$facet=factor(ttl$facet[match(Z$analysis,ttl$analysis)],levels=ttl$facet); Z$y=-log10(Z$P); Z=Z[order(Z$sig),]
  thr=do.call(rbind,lapply(split(Z,Z$facet,drop=TRUE),function(z) if(any(z$sig)) data.frame(facet=z$facet[1],y=-log10(max(z$P[z$sig]))) else NULL))
  v=Z[Z$gene==VG,]; v$lab=sprintf('Vg\nlog2FC %.2f, FDR %s',v$log2FC,fmtp(v$FDR)); cs=Z[Z$cand & Z$sig,]
  i=match(cs$gene,tg$ogs32_id); cs$lab=tg$symbol[i]
  f=ggplot(Z,aes(log2FC,y))+geom_vline(xintercept=0,linewidth=0.25,colour='#b9b8b3')+geom_point(aes(colour=dir),size=0.6,alpha=0.9)+
    {if(!is.null(thr)) geom_hline(data=thr,aes(yintercept=y),linetype='dashed',linewidth=0.3,colour='#52514e')}+
    geom_point(data=Z[Z$cand,],shape=21,size=1.6,stroke=0.35,colour='#52514e',fill=NA)+
    {if(nrow(cs)) geom_text(data=cs,aes(label=lab),size=2,vjust=-0.7,colour='#0b0b0b')}+
    geom_point(data=v,shape=21,size=2.6,stroke=0.8,colour='#0b0b0b',fill='#eda100')+geom_text(data=v,aes(label=lab),size=2.05,hjust=-0.12,vjust=0.5,lineheight=0.9)+
    facet_wrap(~facet,ncol=ncol,scales='free')+scale_colour_manual(values=COL_DIR,name='FDR < 0.05',drop=FALSE)+
    labs(x='log2 fold change',y='-log10 P',title=title,subtitle='Yellow = vitellogenin (GB49544); ringed = curated project target genes (documented NCBI annotation). Dashed line = P at the FDR 0.05 boundary (if any gene passes).',
      caption='Values copied from the DE run tables (edgeR QL quasi-likelihood F-tests). FDR is Benjamini-Hochberg within each contrast. Preliminary 500,000-read subsets.')+
    TH+theme(legend.position='bottom',strip.text=element_text(size=6.8))+guides(colour=guide_legend(override.aes=list(size=2)))
  save_fig(f,name,w,h)}
step('Fig5 volcano, selected primary contrasts',volc(c('GSE168740_subsetA_1d_edgeR_QLF','GSE168740_subsetA_5d_edgeR_QLF','GSE168740_subsetA_10d_edgeR_QLF',
  'PRJNA1017469_primary_ingestion_SUB_vs_CTL_1h_edgeR_QLF','PRJNA1017469_primary_ingestion_D50_vs_CTL_4h_edgeR_QLF','PRJNA1017469_primary_contact_D50_vs_CTL_4h_edgeR_QLF'),
  'Fig5_volcano_selected_primary_contrasts','Volcano plots of selected primary contrasts',180,140,3))
step('FigS1 volcano, all PRJNA1017469 primary',volc(M$analysis[M$dataset=='PRJNA1017469'&M$model=='primary'],'FigS1_volcano_PRJNA1017469_all_primary_contrasts',
  'PRJNA1017469 primary models: all 14 contrasts (contrast 7 of each route is diagnostic/QC)',180,250,3))
## =============== FIGURE 6 (+S2): curated candidate panel heatmap ===============
heat=function(an,name,title,w,h){
  H=CAND[CAND$analysis %in% an,]; keep=unique(H$ogs32_id[H$tested]); H=H[H$ogs32_id %in% keep,]
  H$row=paste0(H$symbol,' (',H$ogs32_id,')'); H$col=factor(H$short,levels=unique(M$short[match(an,M$analysis)]))
  H$tier=factor(H$tier,levels=names(COL_TIER)); H$category=factor(H$category); H$star=ifelse(!is.na(H$FDR)&H$FDR<0.05,'*','')
  H$tierpanel=factor(paste(H$dataset,H$tier,sep='\n'),levels=unique(paste(M$dataset,M$tier,sep='\n')[match(an,M$analysis)]))
  f=ggplot(H,aes(col,row))+geom_tile(aes(fill=pmax(pmin(log2FC,3),-3)),colour='white',linewidth=0.25)+geom_text(aes(label=star),size=3,vjust=0.75)+
    facet_grid(category~tierpanel,scales='free',space='free',labeller=labeller(category=label_wrap_gen(14)))+
    scale_fill_gradient2(low='#2a78d6',mid='#f0efec',high='#e34948',midpoint=0,limits=c(-3,3),na.value='#dcdbd6',name='log2FC\n(capped at +/-3)')+
    labs(x=NULL,y=NULL,title=title,subtitle='Curated project target genes with documented NCBI Gene annotation. * = FDR < 0.05; grey = not tested (removed by the expression filter).',
      caption='Categories are the project\'s curation labels (data/metadata/gene_targets_ncbi.tsv), not enrichment results. Values copied from the DE run tables.')+
    TH+theme(axis.text.x=element_text(angle=60,hjust=1,size=6),axis.text.y=element_text(size=6),strip.text.y=element_text(angle=0,hjust=0,size=6.3),strip.text.x=element_text(size=6.5),panel.grid=element_blank())
  save_fig(f,name,w,h)}
step('Fig6 candidate heatmap, primary',heat(M$analysis[M$tier=='Primary'],'Fig6_candidate_genes_primary_analyses','Curated candidate genes in the primary analyses',180,230))
step('FigS2 candidate heatmap, all 40',heat(AN[order(match(M$tier,names(COL_TIER)),M$panel,M$analysis)],'FigS2_candidate_genes_all_40_analyses','Curated candidate genes in all 40 analyses',260,250))
## ---------------- finish ----------------
writeLines(capture.output(sessionInfo()),file.path(OUT,'R_sessionInfo.txt'))
write.table(data.frame(step=names(STATUS),status=STATUS),file.path(OUT,'STEP_STATUS.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
nerr=sum(STATUS!='OK'); say(sprintf('\nSUMMARY: %d of %d steps OK. Output: %s',sum(STATUS=='OK'),length(STATUS),normalizePath(OUT)))
if(nerr) {say('STEPS WITH ERRORS: ',paste(names(STATUS)[STATUS!='OK'],collapse=', ')); quit(status=1)}
writeLines(format(Sys.time(),'%Y-%m-%dT%H:%M:%S'),file.path(OUT,'_FINISHED'))
