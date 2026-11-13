## PRELIMINARY SUBSAMPLE ANALYSIS: kallisto on the first N Read-1 reads per run (N = 500,000 by default; the actual N per run is read from
## each _COMPLETE.json and printed below). Subset A = reads 1..N; GSE168740 also subset B = reads N+1..2N. NOT full-depth quantification.
## Run manually (not from RUN_ALL_WINDOWS.bat) only after all runs of a study are COMPLETE; a study with any missing run is skipped.
## Transcript->gene: tximport (type='kallisto'), tx2gene from OGS v3.2 GFF3 (mRNA ID -> Parent gene), IDs prefixed 'gnl|Amel_4.5|' to match the
## kallisto target_ids; ignoreTxVersion=FALSE (see scripts/R/tx2gene.R). Every target_id must map, otherwise the script stops.
## Input: abundance.tsv. dropInfReps=TRUE because kallisto was run without bootstraps (-b 0; run_info n_bootstraps = 0), so no
## inferential replicates exist; test_r_environment.R R2.9-R2.11 verifies that HDF5 and TSV imports give identical counts.
##  - edgeR / limma-voom: countsFromAbundance='lengthScaledTPM' (tximport-recommended route for count-based tools);
##  - DESeq2: original estimated counts + transcript-length offsets via DESeqDataSetFromTximport.
## TPM is never used as counts. Filtering: edgeR::filterByExpr. Normalisation: TMM (edgeR/limma), median-of-ratios (DESeq2).
## Dispersion: edgeR estimateDisp + glmQLFit(robust=TRUE); DESeq2 defaults. Multiple testing: Benjamini-Hochberg; threshold FDR<0.05.
rl=file.path('envs','R-library',paste(R.version$major,sub('\\..*','',R.version$minor),sep='.')); if(dir.exists(rl)) .libPaths(c(normalizePath(rl),.libPaths()))
suppressMessages({library(edgeR);library(limma);library(DESeq2);library(tximport)}); source('scripts/R/tx2gene.R')
## Paths (defaults = the real project). Environment overrides exist ONLY for the mock test scripts/validation/test_de_prjna.R (scratch folders).
KDIR=Sys.getenv('DE_KALLISTO_DIR','results/rawreads/kallisto'); QCDEC=Sys.getenv('DE_QC_DECISIONS','data/metadata/PRJNA1017469_library_qc_decisions.tsv')
out=file.path(Sys.getenv('DE_OUT_DIR','results/rawreads/DE'),format(Sys.time(),'run_%Y%m%d_%H%M%S')); dir.create(out,recursive=TRUE,showWarnings=FALSE); sink(file.path(out,'MODEL_LOG.txt'),split=TRUE)
tx2gene=make_tx2gene()
cat('tx2gene rows:',nrow(tx2gene),' genes:',length(unique(tx2gene$gene)),'\n')
targets=read.delim('data/metadata/gene_targets_ncbi.tsv',colClasses='character'); tg=setNames(targets$symbol,targets$ogs32_id)
tg[c('GB47318','GB43727','GB43728','GB41326','GB53576')]=c('Abaecin','CYP9Q2','CYP9Q3','Acph-1','Apisimin')
ena=function(p) read.delim(file.path('data/metadata',paste0('ena_',p,'_runs.tsv')),colClasses='character')
PROJ_OF=c(GSE168740='PRJNA713798',GSE178742='PRJNA740269',PRJNA1017469='PRJNA1017469')
## approved duplicate ENA listings (not downloaded; see data/metadata/run_exclusions.tsv). Absent file = no exclusions.
EXCL=if(file.exists('data/metadata/run_exclusions.tsv')) read.delim('data/metadata/run_exclusions.tsv',colClasses='character',comment.char='#') else data.frame(study=character(),run_accession=character())
cat('approved run exclusions (duplicate ENA listings):',nrow(EXCL),'\n')
jnum=function(txt,k){m=regmatches(txt,regexpr(sprintf('"%s": *[-0-9.eE+]+',k),txt)); if(length(m)) as.numeric(sub('.*: *','',m)) else NA}
imp=function(lab,sub='A'){
  d=list.dirs(file.path(KDIR,lab,sub),recursive=FALSE);d=d[grepl('^[SED]RR[0-9]+$',basename(d))&file.exists(file.path(d,'_COMPLETE.json'))]
  f=file.path(d,'abundance.tsv');f=f[file.exists(f)];names(f)=basename(dirname(f))
  expected=setdiff(ena(PROJ_OF[[lab]])$run_accession,EXCL$run_accession[EXCL$study==lab]); miss=setdiff(expected,names(f))
  if(length(miss)){cat(sprintf('\n%s subset %s: %d of %d runs COMPLETE - study SKIPPED (missing e.g. %s). No partial designs are analysed.\n',lab,sub,length(f),length(expected),paste(head(miss,5),collapse=',')));return(NULL)}
  f=f[expected]
  qc=t(sapply(names(f),function(r){t=paste(readLines(file.path(dirname(f[[r]]),'_COMPLETE.json'),warn=FALSE),collapse=' ');c(reads_quantified=jnum(t,'n_processed'),p_pseudoaligned=jnum(t,'p_pseudoaligned'),reads_target=jnum(t,'reads_target'))}))
  cat(sprintf('\n%s subset %s: reads quantified per run (PRELIMINARY subsample): min %d, median %d, max %d; pseudoaligned %% min %.1f, median %.1f\n',lab,sub,min(qc[,1]),as.integer(median(qc[,1])),max(qc[,1]),min(qc[,2]),median(qc[,2])))
  write.table(cbind(run=rownames(qc),qc),file.path(out,paste0(lab,'_subset',sub,'_reads_per_run.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  check_tx2gene(f,tx2gene)
  list(ls=tximport(f,type='kallisto',tx2gene=tx2gene,countsFromAbundance='lengthScaledTPM',ignoreTxVersion=FALSE,dropInfReps=TRUE),
       raw=tximport(f,type='kallisto',tx2gene=tx2gene,countsFromAbundance='no',ignoreTxVersion=FALSE,dropInfReps=TRUE))}
VG_ID='GB49544'  # Vitellogenin: NCBI Gene 406088 = OGS v3.2 GB49544 (protein identity 0.994, 1770/1770 aa; see PROJECT_STATUS.md)
VG=list()
save_tt=function(tt,tag){tt$symbol=tg[rownames(tt)];write.table(tt,file.path(out,paste0(tag,'.tsv')),sep='\t',quote=FALSE,col.names=NA)
  pick=function(a) {k=intersect(a,names(tt)); if(length(k)) k[1] else NA}; cl=c(pick(c('logFC','log2FoldChange')),pick(c('PValue','P.Value','pvalue')),pick(c('FDR','adj.P.Val','padj')))
  if(VG_ID %in% rownames(tt)){r=tt[VG_ID,]; VG[[tag]]<<-data.frame(analysis=tag,gene=VG_ID,tested=TRUE,logFC=r[[cl[1]]],PValue=r[[cl[2]]],FDR=r[[cl[3]]],genes_tested=nrow(tt),
      call=ifelse(is.na(r[[cl[3]]]),'not tested (NA)',ifelse(r[[cl[3]]]<0.05,'significant (FDR<0.05)','not significant (FDR>=0.05)')))}
  else VG[[tag]]<<-data.frame(analysis=tag,gene=VG_ID,tested=FALSE,logFC=NA,PValue=NA,FDR=NA,genes_tested=nrow(tt),call='not tested (removed by expression filter)')}
qlf=function(cnt,design,coef,tag){y=DGEList(cnt);k=filterByExpr(y,design=design);y=calcNormFactors(y[k,,keep.lib.sizes=FALSE]);y=estimateDisp(y,design)
  fit=glmQLFit(y,design,robust=TRUE);tt=topTags(glmQLFTest(fit,coef=coef),n=Inf,adjust.method='BH')$table;save_tt(tt,paste0(tag,'_edgeR_QLF'))
  cat(tag,': genes tested',nrow(tt),'| FDR<0.05:',sum(tt$FDR<0.05),'| common dispersion',signif(y$common.dispersion,3),'\n');invisible(y)}
showdes=function(d,tag){cat('\nDesign matrix',tag,'\n');print(d)}
## ---------- GSE168740 (pools of 15 brains; unit = pooled library; n = 3 v 3 per age) ----------
md=read.delim('data/metadata/GSE168740_sample_metadata_verified.tsv',colClasses='character'); e=ena('PRJNA713798')
for(sub in c('A','B')){x=imp('GSE168740',sub); if(is.null(x)) next
  lab=md$title[match(e$sample_alias[match(colnames(x$ls$counts),e$run_accession)],md$gsm)]
  colnames(x$ls$counts)=colnames(x$raw$counts)=colnames(x$raw$length)=colnames(x$raw$abundance)=lab
  write.table(round(x$ls$counts),file.path(out,paste0('GSE168740_subset',sub,'_lengthScaledTPM_counts.tsv')),sep='\t',quote=FALSE,col.names=NA)
  for(age in c('1d','5d','10d')) for(drop in c(FALSE,TRUE)){ if(drop && age!='1d') next
    cs=grep(paste0('_',age,'$'),lab,value=TRUE); if(drop) cs=setdiff(cs,'C1_1d'); trt=factor(substr(cs,1,1),levels=c('C','T')); des=model.matrix(~trt)
    tag=paste0('GSE168740_subset',sub,'_',age,ifelse(drop,'_noC1_1d','')); if(sub=='A') showdes(des,tag)
    qlf(x$ls$counts[,cs],des,2,tag)
    if(sub=='A'){txi=list(counts=x$raw$counts[,cs],length=x$raw$length[,cs],abundance=x$raw$abundance[,cs],countsFromAbundance='no')
      dds=DESeqDataSetFromTximport(txi,data.frame(trt=trt,row.names=cs),~trt); dds=dds[rowSums(counts(dds)>=10)>=3,]
      r=as.data.frame(results(DESeq(dds,quiet=TRUE),name='trt_T_vs_C',alpha=0.05)); save_tt(r,paste0(tag,'_DESeq2')); cat(tag,'DESeq2 padj<0.05:',sum(r$padj<0.05,na.rm=TRUE),'\n')}}}
## ---------- GSE178742 ----------
mm=read.delim('data/metadata/GSE178742_sample_metadata_verified.tsv',colClasses='character'); e=ena('PRJNA740269'); x=imp('GSE178742','A')
if(!is.null(x)){d0=mm[match(e$sample_alias[match(colnames(x$ls$counts),e$run_accession)],mm$gsm),]; rownames(d0)=colnames(x$ls$counts)
  write.table(cbind(run=rownames(d0),d0),file.path(out,'GSE178742_run_to_sample_map.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(round(x$ls$counts),file.path(out,'GSE178742_subsetA_lengthScaledTPM_counts.tsv'),sep='\t',quote=FALSE,col.names=NA)
  for(ex in sort(unique(na.omit(d0$experiment)))){cs=rownames(d0)[which(d0$experiment==ex)]; d=d0[cs,]; trt=factor(d$treatment,levels=c('Unexposed','Exposed')); des=model.matrix(~trt); cnt=x$ls$counts[,cs]
    if(grepl('Apiary',ex)){
      ## (1) PRIMARY: colony pseudo-bulk (sum of 5 bees per colony) -> n = 3 exposed v 3 control colonies
      tc=tapply(d$treatment,d$colony,function(z) length(unique(z))); if(any(tc!=1)) stop('colony with mixed treatment labels: ',paste(names(tc)[tc!=1],collapse=','))
      pb=sapply(split(cs,d$colony),function(z) rowSums(cnt[,z,drop=FALSE])); ct=factor(sapply(colnames(pb),function(c) unique(d$treatment[d$colony==c])),levels=c('Unexposed','Exposed'))
      cat(ex,': bees per colony:',paste(names(table(d$colony)),table(d$colony),sep='=',collapse=' '),'| colonies per arm:',paste(names(table(ct)),table(ct),sep='=',collapse=' '),'\n')
      dp=model.matrix(~ct); rownames(dp)=colnames(pb); showdes(dp,paste(ex,'colony pseudo-bulk')); qlf(pb,dp,2,paste0('GSE178742_',ex,'_colonyPseudobulk'))
      ## (2) SENSITIVITY: bee-level limma-voom with colony as blocking factor (duplicateCorrelation)
      y=DGEList(cnt); y=calcNormFactors(y[filterByExpr(y,design=des),,keep.lib.sizes=FALSE]); v=voom(y,des); dc=duplicateCorrelation(v,des,block=d$colony)
      v=voom(y,des,block=d$colony,correlation=dc$consensus); dc=duplicateCorrelation(v,des,block=d$colony)
      fit=eBayes(lmFit(v,des,block=d$colony,correlation=dc$consensus)); tt=topTable(fit,coef=2,n=Inf,adjust.method='BH',confint=TRUE); save_tt(tt,paste0('GSE178742_',ex,'_voom_colonyBlock'))
      cat(ex,'voom+duplicateCorrelation: consensus intra-colony r =',round(dc$consensus,3),'| FDR<0.05:',sum(tt$adj.P.Val<0.05),'\n')
    } else {showdes(des,ex); qlf(cnt,des,2,paste0('GSE178742_',ex))}}}
## ---------- PRJNA1017469 (approved model, FINAL PROJECT APPROVAL 2026-11-12) ----------
## 47 libraries (one ENA run each; the 47 duplicate listings are excluded upstream). Isolate = <time>-<group>-<rep>; each library = pool of 3 heads.
## Two exposure routes analysed SEPARATELY (ingestion: CTL-IG / IGSUB / IG50; contact: CTL-CONT / CONTSUB / CONT50), each with its own controls.
## Model per route: cell means ~0+group, group = dose x time (6 levels); edgeR filterByExpr(group) -> TMM -> estimateDisp -> glmQLFit(robust) ->
## glmQLFTest per contrast; BH per contrast. Contrasts: each dose vs time-matched control (1 h, 4 h), dose x time interactions, control drift.
## Cross-route comparison: EXPLORATORY only (joint 12-group model, difference of dose effects between routes).
## Groups: CTL = route-specific control; SUB = sublethal dose (IGSUB/CONTSUB); D50 = the IG50/CONT50 groups (the approval calls them the lethal dose).
## Contrasts per route (mu = group means CTL_1h, SUB_1h, D50_1h, CTL_4h, SUB_4h, D50_4h; approved 2026-11-12):
##   TREATMENT (6): SUB_vs_CTL_1h = mu2-mu1; D50_vs_CTL_1h = mu3-mu1; SUB_vs_CTL_4h = mu5-mu4; D50_vs_CTL_4h = mu6-mu4;
##                  SUBxTime_interaction = (mu5-mu4)-(mu2-mu1) = SUB_vs_CTL_4h - SUB_vs_CTL_1h; D50xTime_interaction = D50_vs_CTL_4h - D50_vs_CTL_1h.
##   DIAGNOSTIC/QC (1): CTL_4h_vs_1h_drift = mu4-mu1 (change in controls between sampling times; NOT a treatment effect, never reported as one).
##   The 7 contrasts span 5 dimensions (6 groups -> 5 independent comparisons); the interactions are differences of the dose contrasts.
##   Each contrast is a separate QL F-test; Benjamini-Hochberg FDR is applied across GENES WITHIN each contrast, not across contrasts.
## Vitellogenin: GB49544 only (approved 2026-11-12). GB13999 is an unverified historical alias and is not used.
## LIMITATION: Petri-dish and colony allocation of the pooled bees is not reported; pools are treated as replicates, but dish-level independence
## is NOT established and pseudoreplication cannot be excluded. Results are labelled accordingly.
## Library QC decisions: data/metadata/PRJNA1017469_library_qc_decisions.tsv. While any decision is PENDING, the primary DE is NOT run.
x=imp('PRJNA1017469','A'); if(!is.null(x)){
  e=ena('PRJNA1017469'); a=read.delim('data/metadata/ena_PRJNA1017469_sample_attributes.tsv',colClasses='character',check.names=FALSE)
  runs=colnames(x$ls$counts); iso=a$isolate[match(runs,a$run_accession)]; libn=e$experiment_alias[match(runs,e$run_accession)]
  if(any(is.na(iso))||any(duplicated(iso))) stop('PRJNA1017469: isolate missing or duplicated among the quantified runs')
  p=strsplit(iso,'-'); tm=sapply(p,`[`,1); repl=sapply(p,function(z) z[length(z)]); grp=sapply(strsplit(iso,'-'),function(z) paste(z[2:(length(z)-1)],collapse='-'))
  route=ifelse(grepl('IG',grp),'ingestion',ifelse(grepl('CONT',grp),'contact',NA)); dose=ifelse(grepl('^CTL',grp),'CTL',ifelse(grepl('SUB$',grp),'SUB',ifelse(grepl('50$',grp),'D50',NA)))
  if(any(is.na(route))||any(is.na(dose))||!all(tm %in% c('1h','4h'))) stop('PRJNA1017469: unparsed isolate labels: ',paste(iso[is.na(route)|is.na(dose)|!tm %in% c('1h','4h')],collapse=','))
  qcj=t(sapply(runs,function(r){t=paste(readLines(file.path(KDIR,'PRJNA1017469','A',r,'_COMPLETE.json'),warn=FALSE),collapse=' ')
    c(reads_quantified=jnum(t,'n_processed'),n_pseudoaligned=jnum(t,'n_pseudoaligned'),p_pseudoaligned=jnum(t,'p_pseudoaligned'))}))
  dec=read.delim(QCDEC,colClasses='character',comment.char='#'); cat('library QC decisions file:',QCDEC,'\n')
  decision=ifelse(runs %in% dec$run_accession,dec$decision[match(runs,dec$run_accession)],'RETAIN')
  if(!all(decision %in% c('RETAIN','EXCLUDE_PRIMARY','PENDING'))) stop('PRJNA1017469: unknown decision value in library_qc_decisions.tsv')
  S=data.frame(run=runs,library=libn,isolate=iso,route=route,time=tm,dose=dose,replicate=repl,group=paste(dose,tm,sep='_'),qcj,decision=decision,
               ena_read_count=e$read_count[match(runs,e$run_accession)],stringsAsFactors=FALSE,row.names=runs)
  S$p_pseudoaligned_vs_route_median=round(S$p_pseudoaligned/ave(S$p_pseudoaligned,S$route,FUN=median),3)
  write.table(S,file.path(out,'PRJNA1017469_sample_sheet_and_QC.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(round(x$ls$counts),file.path(out,'PRJNA1017469_subsetA_lengthScaledTPM_counts.tsv'),sep='\t',quote=FALSE,col.names=NA)
  cat('\nPRJNA1017469 libraries by route x group x time (all quantified):\n'); print(table(paste(S$route,S$dose),S$time))
  cat('\nPRJNA1017469 libraries with a recorded QC decision:\n'); print(S[S$run %in% dec$run_accession,c('library','isolate','reads_quantified','n_pseudoaligned','p_pseudoaligned','p_pseudoaligned_vs_route_median','decision')],row.names=FALSE)
  LEV=c('CTL_1h','SUB_1h','D50_1h','CTL_4h','SUB_4h','D50_4h')
  CON=rbind(SUB_vs_CTL_1h=c(-1,1,0,0,0,0), D50_vs_CTL_1h=c(-1,0,1,0,0,0), SUB_vs_CTL_4h=c(0,0,0,-1,1,0), D50_vs_CTL_4h=c(0,0,0,-1,0,1),
            SUBxTime_interaction=c(1,-1,0,-1,1,0), D50xTime_interaction=c(1,0,-1,-1,0,1), CTL_4h_vs_1h_drift=c(-1,0,0,1,0,0)); colnames(CON)=LEV
  CTYPE=c(rep('treatment',6),'diagnostic/QC'); names(CTYPE)=rownames(CON)
  cat('\nContrast matrix (columns = group means, per route; last row = diagnostic/QC, not a treatment effect; BH FDR across genes within each contrast):\n'); print(cbind(as.data.frame(CON),type=CTYPE))
  ## estimability: full column rank, residual df > 0, every group present with >= 3 libraries (approved minimum), each contrast in the row space
  checkdes=function(keep,rt,tag){s=S[keep & S$route==rt,]; g=factor(s$group,levels=LEV); des=model.matrix(~0+g); colnames(des)=LEV; rownames(des)=s$isolate
    n=table(g); rk=qr(des)$rank; est=apply(CON,1,function(cc){X=rbind(des,cc); qr(X)$rank==rk})
    ok=rk==ncol(des) && nrow(des)-rk>0 && all(n>=3) && all(est) && is.null(limma::nonEstimable(des))
    data.frame(model=tag,route=rt,n_libraries=nrow(des),rank=rk,n_coef=ncol(des),residual_df=nrow(des)-rk,group_sizes=paste(names(n),n,sep='=',collapse=';'),
               min_group=min(n),contrasts_estimable=paste0(sum(est),'/',length(est)),PASS=ok)}
  fitroute=function(keep,rt,tag){s=S[keep & S$route==rt,]; g=factor(s$group,levels=LEV); des=model.matrix(~0+g); colnames(des)=LEV; rownames(des)=s$isolate
    showdes(des,paste(tag,rt)); y=DGEList(x$ls$counts[,s$run],group=g); k=filterByExpr(y,group=g); y=calcNormFactors(y[k,,keep.lib.sizes=FALSE]); y=estimateDisp(y,des)
    fit=glmQLFit(y,des,robust=TRUE); cat(sprintf('%s %s: libraries %d, genes after filterByExpr %d of %d, common dispersion %.3g, prior QL df (median) %.3g, residual df %d\n',
      tag,rt,ncol(y),nrow(y),length(k),y$common.dispersion,median(fit$df.prior),nrow(des)-ncol(des)))
    ids=colnames(y); if(!identical(ids,s$run)||anyDuplicated(ids)||!all(ids %in% S$run)) stop('PRJNA1017469 ',tag,' ',rt,': count-matrix columns do not match the sample table runs')
    write.table(data.frame(isolate=S$isolate[match(ids,S$run)],run=ids,lib.size=y$samples$lib.size,norm.factors=y$samples$norm.factors),file.path(out,paste0('PRJNA1017469_',tag,'_',rt,'_TMM_factors.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
    pdf(file.path(out,paste0('PRJNA1017469_',tag,'_',rt,'_MDS_BCV.pdf'))); plotMDS(y,labels=s$isolate,col=as.integer(g),main=paste(tag,rt,'MDS (logCPM, top 500)')); plotBCV(y,main=paste(tag,rt)); plotQLDisp(fit); dev.off()
    for(cn in rownames(CON)){tt=topTags(glmQLFTest(fit,contrast=CON[cn,]),n=Inf,adjust.method='BH')$table; save_tt(tt,paste0('PRJNA1017469_',tag,'_',rt,'_',cn,'_edgeR_QLF'))
      cat(sprintf('  %-22s [%s] FDR<0.05: %d (up %d, down %d)\n',cn,CTYPE[cn],sum(tt$FDR<0.05),sum(tt$FDR<0.05&tt$logFC>0),sum(tt$FDR<0.05&tt$logFC<0)))}}
  primary=S$decision=='RETAIN'; sens_all=rep(TRUE,nrow(S)); pend=S$decision=='PENDING'
  chk=rbind(checkdes(primary|pend,'ingestion','primary_if_PENDING_retained'),checkdes(primary|pend,'contact','primary_if_PENDING_retained'),
            checkdes(primary,'ingestion','primary_if_PENDING_excluded'),checkdes(primary,'contact','primary_if_PENDING_excluded'),
            checkdes(sens_all,'ingestion','sensitivity_all_libraries'),checkdes(sens_all,'contact','sensitivity_all_libraries'))
  cat('\nDesign checks (rank, residual df, group sizes, contrast estimability):\n'); print(chk,row.names=FALSE)
  write.table(chk,file.path(out,'PRJNA1017469_design_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  if(any(pend)){cat('\nPRJNA1017469: QC decision PENDING for',paste(S$library[pend],collapse=', '),'- primary DE NOT run.',
      'Record RETAIN or EXCLUDE_PRIMARY in data/metadata/PRJNA1017469_library_qc_decisions.tsv (with the evidence above) and rerun.\n')
  } else {
    pc=chk[chk$model=='primary_if_PENDING_excluded',]   # no PENDING left: this is the primary set (RETAIN only)
    if(!all(pc$PASS)) stop('PRJNA1017469 primary design check FAILED - see PRJNA1017469_design_checks.tsv; no results written for the primary model')
    for(rt in c('ingestion','contact')) fitroute(primary,rt,'primary')
    ## SENSITIVITY: all 47 libraries (QC exclusions reversed), only where the design check passes
    for(rt in c('ingestion','contact')) if(chk$PASS[chk$model=='sensitivity_all_libraries'&chk$route==rt]) fitroute(sens_all,rt,'sensitivity_allLibraries')
    ## EXPLORATORY cross-route: joint model with 12 route x dose x time cells; (dose effect in ingestion) - (dose effect in contact), per time
    s=S[primary,]; g=factor(paste(s$route,s$group,sep='.')); des=model.matrix(~0+g); colnames(des)=levels(g)
    y=DGEList(x$ls$counts[,s$run],group=g); y=calcNormFactors(y[filterByExpr(y,group=g),,keep.lib.sizes=FALSE]); y=estimateDisp(y,des); fit=glmQLFit(y,des,robust=TRUE)
    for(d in c('SUB','D50')) for(t in c('1h','4h')){cc=setNames(rep(0,ncol(des)),colnames(des))
      cc[paste0('ingestion.',d,'_',t)]=1; cc[paste0('ingestion.CTL_',t)]=-1; cc[paste0('contact.',d,'_',t)]=-1; cc[paste0('contact.CTL_',t)]=1
      tt=topTags(glmQLFTest(fit,contrast=cc),n=Inf,adjust.method='BH')$table; save_tt(tt,paste0('PRJNA1017469_EXPLORATORY_routeDifference_',d,'_',t,'_edgeR_QLF'))
      cat(sprintf('EXPLORATORY route difference %s %s: FDR<0.05: %d\n',d,t,sum(tt$FDR<0.05)))}}}
## ---------- Vitellogenin (GB49544) across every analysis written above ----------
if(length(VG)){v=do.call(rbind,VG); v$primary_or_secondary=ifelse(grepl('_drift_',v$analysis),'diagnostic/QC (not a treatment effect)',ifelse(grepl('EXPLORATORY',v$analysis),'exploratory',ifelse(grepl('DESeq2|voom|sensitivity|noC1_1d|subsetB',v$analysis),'secondary/sensitivity','primary')))
  write.table(v,file.path(out,'Vg_GB49544_all_analyses.tsv'),sep='\t',quote=FALSE,row.names=FALSE); cat('\nVitellogenin (GB49544) in each analysis:\n'); print(v[,c('analysis','logFC','PValue','FDR','call')],row.names=FALSE)}
cat('\nDE script finished; outputs in',out,'\n'); sink(); si=capture.output(sessionInfo()); writeLines(si,file.path(out,'R_sessionInfo.txt'))
if(Sys.getenv('DE_OUT_DIR')=='') writeLines(si,'logs/R_sessionInfo.txt')
writeLines(format(Sys.time(),'%Y-%m-%dT%H:%M:%S'),file.path(out,'_DE_FINISHED'))
