## Mock test of the PRJNA1017469 model in scripts/R/de_from_kallisto.R (FINAL PROJECT APPROVAL 2026-11-12).
## Run from the project folder:  Rscript scripts/validation/test_de_prjna.R      (or double-click RUN_DE_MOCK_TEST.bat)
## Writes ONLY to preflight_scratch/de_mock_<time>/. Reads the real metadata (ENA tables, sample attributes, exclusion list,
## QC decisions) and the OGS v3.2 GFF3 (tx2gene). Synthetic kallisto outputs with KNOWN effects are generated for the 47 kept runs:
##   Vg GB49544: +2 log2 in ingestion D50 at 4 h only (true positive); 150 genes: +/-1.5 log2 in contact SUB at 1 h; all else null.
## The real DE script is then run for scenarios A-F (A-E use a fixed PENDING fixture; F uses the LIVE decisions file) with DE_KALLISTO_DIR / DE_OUT_DIR / DE_QC_DECISIONS pointing to scratch.
rl=file.path('envs','R-library',paste(R.version$major,sub('\\..*','',R.version$minor),sep='.')); if(dir.exists(rl)) .libPaths(c(normalizePath(rl),.libPaths()))
source('scripts/R/tx2gene.R'); set.seed(20261112)
SCR=file.path('preflight_scratch',format(Sys.time(),'de_mock_%Y%m%d_%H%M%S')); dir.create(SCR,recursive=TRUE)
LOG=file.path(SCR,'test_de_prjna_log.txt'); RES=logical(0)
say=function(...) {m=paste0(...); cat(m,'\n'); cat(m,'\n',file=LOG,append=TRUE)}
chk=function(name,ok,detail='') {ok=isTRUE(ok); RES[name]<<-ok; say(if(ok) 'PASS ' else 'FAIL ',name,if(nzchar(detail)) paste0('   [',detail,']') else '')}
shaf=function(f) if(.Platform$OS.type=='windows') tryCatch(system2('certutil',c('-hashfile',shQuote(f),'SHA256'),stdout=TRUE)[2],error=function(e) 'n/a') else tryCatch(strsplit(system2('sha256sum',f,stdout=TRUE),' ')[[1]][1],error=function(e) 'n/a')
say('test_de_prjna.R  ',format(Sys.time()),'  R ',R.version.string,'  de_from_kallisto.R sha256: ',shaf('scripts/R/de_from_kallisto.R'))
## ---- synthetic kallisto outputs
t2g=make_tx2gene(); one=t2g[!duplicated(t2g$gene),]
VG='GB49544'; if(!VG %in% one$gene) stop('GB49544 not in tx2gene')
genes=unique(c(VG,sample(setdiff(one$gene,VG),2999))); tx=one$tx[match(genes,one$gene)]
e=read.delim('data/metadata/ena_PRJNA1017469_runs.tsv',colClasses='character'); a=read.delim('data/metadata/ena_PRJNA1017469_sample_attributes.tsv',colClasses='character',check.names=FALSE)
ex=read.delim('data/metadata/run_exclusions.tsv',colClasses='character',comment.char='#'); kept=setdiff(e$run_accession,ex$run_accession)
chk('M0 47 kept runs after the approved exclusions',length(kept)==47,length(kept))
iso=a$isolate[match(kept,a$run_accession)]; p=strsplit(iso,'-'); tm=sapply(p,`[`,1); grp=sapply(p,function(z) paste(z[2:(length(z)-1)],collapse='-'))
route=ifelse(grepl('IG',grp),'ingestion','contact'); dose=ifelse(grepl('^CTL',grp),'CTL',ifelse(grepl('SUB$',grp),'SUB','D50'))
mu0=exp(rnorm(length(genes),log(80),1.6)); mu0[1]=500; len=round(runif(length(genes),600,4000)); phi=0.05
de150=sample(setdiff(which(mu0>=50),1),150); lfc150=sample(c(-1.5,1.5),150,TRUE)
K=file.path(SCR,'kallisto'); dir.create(file.path(K,'PRJNA1017469','A'),recursive=TRUE)
for(i in seq_along(kept)){lf=rep(0,length(genes))
  if(route[i]=='ingestion'&&dose[i]=='D50'&&tm[i]=='4h') lf[1]=2
  if(route[i]=='contact'&&dose[i]=='SUB'&&tm[i]=='1h') lf[de150]=lfc150
  sz=exp(rnorm(1,0,0.2)); cnt=rnbinom(length(genes),mu=mu0*2^lf*sz,size=1/phi)
  d=file.path(K,'PRJNA1017469','A',kept[i]); dir.create(d)
  tpm=cnt/len; tpm=1e6*tpm/sum(tpm)
  write.table(data.frame(target_id=tx,length=len,eff_length=len-150,est_counts=cnt,tpm=tpm),file.path(d,'abundance.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  writeLines(sprintf('{"run": "%s", "n_processed": %d, "n_pseudoaligned": %d, "p_pseudoaligned": %.1f}',kept[i],round(sum(cnt)/0.6),sum(cnt),60),file.path(d,'_COMPLETE.json'))}
dec0=read.delim('data/metadata/PRJNA1017469_library_qc_decisions.tsv',colClasses='character',comment.char='#')   # LIVE decisions (used in scenario F)
## Fixed fixture for scenarios A-E, independent of later QC decisions: 10/12/33 EXCLUDE_PRIMARY, 17/28/45 PENDING (the approved state before Final QC).
decPend=dec0; decPend$decision[decPend$library %in% c('library_17','library_28','library_45')]='PENDING'; decPend$decision[decPend$library %in% c('library_10','library_12','library_33')]='EXCLUDE_PRIMARY'
wdec=function(d,f){writeLines('# mock decisions',f); suppressWarnings(write.table(d,f,sep='\t',quote=FALSE,row.names=FALSE,append=TRUE)); f}
runde=function(tag,decfile,kdir=K){o=file.path(SCR,paste0('out_',tag)); dir.create(o)
  Sys.setenv(DE_KALLISTO_DIR=kdir,DE_OUT_DIR=o,DE_QC_DECISIONS=decfile)    # inherited by the child Rscript (system2 env= is Unix-only)
  rc=system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/R/de_from_kallisto.R'),stdout=file.path(SCR,paste0('stdout_',tag,'.txt')),stderr=file.path(SCR,paste0('stderr_',tag,'.txt')))
  Sys.unsetenv(c('DE_KALLISTO_DIR','DE_OUT_DIR','DE_QC_DECISIONS'))
  r=list.dirs(o,recursive=FALSE); list(rc=rc,dir=if(length(r)) r[1] else NA,files=if(length(r)) list.files(r[1]) else character(0),
       log=if(length(r)&&file.exists(file.path(r[1],'MODEL_LOG.txt'))) readLines(file.path(r[1],'MODEL_LOG.txt'),warn=FALSE) else character(0))}
## ---- A: fixture with libraries 17, 28 and 45 PENDING -> QC + design checks only, no primary DE
A=runde('A_pending',wdec(decPend,file.path(SCR,'dec_A_fixture.tsv')))
chk('A1 runs to completion (exit 0, _DE_FINISHED written)',A$rc==0 && '_DE_FINISHED' %in% A$files,paste('rc',A$rc))
chk('A2 other two studies skipped (no outputs present), not analysed partially',any(grepl('GSE168740 subset A: 0 of 18 runs COMPLETE - study SKIPPED',A$log)) && any(grepl('GSE178742 subset A: 0 of 79 runs COMPLETE - study SKIPPED',A$log)))
chk('A3 PRJNA1017469 imported with 47 runs (duplicate listings not expected)','PRJNA1017469_sample_sheet_and_QC.tsv' %in% A$files && !any(grepl('PRJNA1017469 subset A: .* SKIPPED',A$log)))
chk('A4 PENDING -> primary DE NOT run (no primary/sensitivity/exploratory tables)',!any(grepl('_primary_|sensitivity|EXPLORATORY',A$files)) && any(grepl('primary DE NOT run',A$log)))
dc=read.delim(file.path(A$dir,'PRJNA1017469_design_checks.tsv'))
expd=data.frame(model=rep(c('primary_if_PENDING_retained','primary_if_PENDING_excluded','sensitivity_all_libraries'),each=2),route=rep(c('ingestion','contact'),3),
               n=c(22,22,20,21,24,23),df=c(16,16,14,15,18,17),minG=c(3,3,3,3,4,3))   # library_17 PENDING (2026-11-12): contact without pending = 21 libraries
got=dc[,c('model','route','n_libraries','residual_df','min_group')]; names(got)=names(expd)
chk('A5 design checks: n, residual df and minimum group size as derived by hand from the isolate table',isTRUE(all.equal(got,expd,check.attributes=FALSE)),paste(paste(got$n,got$df,got$minG,sep='/'),collapse=' '))
chk('A6 all six designs full rank (6 coefficients) with 7/7 contrasts estimable',all(dc$rank==6 & dc$contrasts_estimable=='7/7' & dc$PASS))
ss=read.delim(file.path(A$dir,'PRJNA1017469_sample_sheet_and_QC.tsv'),colClasses='character')
chk('A7 sample sheet: 47 rows; libraries 10/12/33 EXCLUDE_PRIMARY; 17/45/28 PENDING; all others RETAIN',nrow(ss)==47 &&
    setequal(ss$library[ss$decision=='EXCLUDE_PRIMARY'],c('library_10','library_12','library_33')) && setequal(ss$library[ss$decision=='PENDING'],c('library_17','library_45','library_28')) &&
    all(ss$decision[!ss$library %in% c('library_10','library_12','library_33','library_17','library_45','library_28')]=='RETAIN'))
## ---- B: PENDING resolved as RETAIN (mock only) -> full model
dB=decPend; dB$decision[dB$decision=='PENDING']='RETAIN'; B=runde('B_retain',wdec(dB,file.path(SCR,'dec_B.tsv')))
pt=grep('^PRJNA1017469_primary_(ingestion|contact)_.*_edgeR_QLF\\.tsv$',B$files,value=TRUE)
chk('B1 runs to completion; 14 primary tables (7 contrasts x 2 routes)',B$rc==0 && '_DE_FINISHED' %in% B$files && length(pt)==14,paste(length(pt),'tables'))
chk('B2 sensitivity (all libraries) 14 tables and exploratory cross-route 4 tables',sum(grepl('^PRJNA1017469_sensitivity_allLibraries_.*QLF\\.tsv$',B$files))==14 && sum(grepl('EXPLORATORY_routeDifference',B$files))==4)
rd=function(f) read.delim(file.path(B$dir,f),row.names=1)
t1=rd('PRJNA1017469_primary_ingestion_D50_vs_CTL_4h_edgeR_QLF.tsv'); t2=rd('PRJNA1017469_primary_contact_D50_vs_CTL_4h_edgeR_QLF.tsv'); t3=rd('PRJNA1017469_primary_ingestion_D50xTime_interaction_edgeR_QLF.tsv')
chk('B3 known effect recovered: Vg up in ingestion D50 vs CTL at 4 h (FDR<0.05, logFC 1.5-2.5)',t1[VG,'FDR']<0.05 && t1[VG,'logFC']>1.5 && t1[VG,'logFC']<2.5,sprintf('logFC %.2f FDR %.2g',t1[VG,'logFC'],t1[VG,'FDR']))
chk('B4 no Vg effect where none was simulated (contact D50 4 h: FDR>=0.05)',t2[VG,'FDR']>=0.05,sprintf('logFC %.2f FDR %.2g',t2[VG,'logFC'],t2[VG,'FDR']))
chk('B5 interaction contrast detects the 4 h-only Vg effect (positive)',t3[VG,'FDR']<0.05 && t3[VG,'logFC']>0,sprintf('logFC %.2f FDR %.2g',t3[VG,'logFC'],t3[VG,'FDR']))
t4=rd('PRJNA1017469_primary_contact_SUB_vs_CTL_1h_edgeR_QLF.tsv'); sig=rownames(t4)[t4$FDR<0.05]; tp=genes[de150]
fdp=if(length(sig)) mean(!sig %in% tp) else 0
chk('B6 contact SUB 1 h: >= 50% of the 150 simulated genes found, false discovery proportion < 0.15',mean(tp %in% sig)>=0.5 && fdp<0.15,sprintf('found %d/150, FDP %.3f',sum(tp %in% sig),fdp))
nul=sapply(c('ingestion_SUB_vs_CTL_1h','ingestion_SUB_vs_CTL_4h','ingestion_D50_vs_CTL_1h','contact_D50_vs_CTL_1h','contact_SUB_vs_CTL_4h','ingestion_CTL_4h_vs_1h_drift','contact_CTL_4h_vs_1h_drift'),
           function(k) sum(rd(paste0('PRJNA1017469_primary_',k,'_edgeR_QLF.tsv'))$FDR<0.05))
chk('B7 null contrasts: at most 2 FDR<0.05 genes each (none simulated)',all(nul<=2),paste(nul,collapse=','))
vg=read.delim(file.path(B$dir,'Vg_GB49544_all_analyses.tsv'))
chk('B8 Vg summary covers every primary table; labels primary/secondary/exploratory; drift labelled diagnostic/QC, never primary',all(sub('\\.tsv$','',pt) %in% vg$analysis) &&
    all(c('primary','secondary/sensitivity','exploratory') %in% vg$primary_or_secondary) &&
    all(vg$primary_or_secondary[grepl('_drift_',vg$analysis)]=='diagnostic/QC (not a treatment effect)') && !any(grepl('_drift_',vg$analysis[vg$primary_or_secondary=='primary'])))
tm1=read.delim(file.path(B$dir,'PRJNA1017469_primary_ingestion_TMM_factors.tsv'))
chk('B9 primary ingestion model excludes libraries 10 and 12 and retains 45 and 28 (mock decision)',nrow(tm1)==22 && !any(c('1h-IG50-3','1h-IGSUB-1') %in% tm1$isolate) && all(c('4h-CTL-IG-2','4h-IGSUB-1') %in% tm1$isolate))
## ---- E: primary DE blocked while ANY ONE of library_17, library_28, library_45 is still PENDING (each tested alone)
for(lib in c('library_17','library_28','library_45')){dE=dB; dE$decision[dE$library==lib]='PENDING'
  E=runde(paste0('E_pending_',lib),wdec(dE,file.path(SCR,paste0('dec_E_',lib,'.tsv'))))
  chk(paste0('E ',lib,' alone PENDING -> primary DE NOT run (no primary, sensitivity or exploratory tables; log names ',lib,'); run ends normally'),
      E$rc==0 && '_DE_FINISHED' %in% E$files && !any(grepl('_primary_|sensitivity|EXPLORATORY',E$files)) &&
      any(grepl(paste0('QC decision PENDING for ',lib,' - primary DE NOT run'),E$log)) && 'PRJNA1017469_design_checks.tsv' %in% E$files)}
## ---- C: an exclusion that leaves a group with 2 libraries -> must stop before writing primary results
dC=dB; dC=rbind(dC,data.frame(study='PRJNA1017469',run_accession=kept[iso=='1h-CONT50-2'],library='x',isolate='1h-CONT50-2',decision='EXCLUDE_PRIMARY',reason='mock',decided_by='mock',date='mock'))
C=runde('C_toofew',wdec(dC,file.path(SCR,'dec_C.tsv')))
chk('C1 group below 3 libraries -> script stops (non-zero exit), no primary tables, design checks written',C$rc!=0 && !any(grepl('_primary_.*QLF',C$files)) && 'PRJNA1017469_design_checks.tsv' %in% C$files && !'_DE_FINISHED' %in% C$files,paste('rc',C$rc))
## ---- D: one kept run missing -> study skipped (no partial design)
K2=file.path(SCR,'kallisto_missing'); dir.create(file.path(K2,'PRJNA1017469','A'),recursive=TRUE)
for(r in kept[-5]) file.copy(file.path(K,'PRJNA1017469','A',r),file.path(K2,'PRJNA1017469','A'),recursive=TRUE)
D=runde('D_missing',file.path(SCR,'dec_B.tsv'),K2)
chk('D1 a missing kept run -> PRJNA1017469 skipped, nothing analysed',D$rc==0 && any(grepl('PRJNA1017469 subset A: 46 of 47 runs COMPLETE - study SKIPPED',D$log)) && !any(grepl('^PRJNA1017469_',D$files)))
## ---- F: the LIVE decisions file (Final QC 2026-11-13: 10/12/33/45 EXCLUDE_PRIMARY, 17/28 RETAIN, nothing PENDING) -> primary model runs on exactly the eligible libraries
F=runde('F_live','data/metadata/PRJNA1017469_library_qc_decisions.tsv')
liveP=dec0$library[dec0$decision=='PENDING']; liveX=sort(dec0$isolate[dec0$decision=='EXCLUDE_PRIMARY'])
chk('F0 live decisions file: no PENDING library; primary exclusions = 1h-IG50-3, 1h-IGSUB-1, 1h-CTL-CONT-2, 4h-CTL-IG-2',length(liveP)==0 && identical(liveX,sort(c('1h-IG50-3','1h-IGSUB-1','1h-CTL-CONT-2','4h-CTL-IG-2'))),paste(liveX,collapse=','))
pf=grep('^PRJNA1017469_primary_(ingestion|contact)_.*_edgeR_QLF\\.tsv$',F$files,value=TRUE)
chk('F1 runs to completion; 14 primary tables, 14 sensitivity tables, 4 exploratory tables',F$rc==0 && '_DE_FINISHED' %in% F$files && length(pf)==14 && sum(grepl('sensitivity_allLibraries_.*QLF',F$files))==14 && sum(grepl('EXPLORATORY_routeDifference',F$files))==4,paste(length(pf),'primary'))
ti=read.delim(file.path(F$dir,'PRJNA1017469_primary_ingestion_TMM_factors.tsv')); tc=read.delim(file.path(F$dir,'PRJNA1017469_primary_contact_TMM_factors.tsv'))
chk('F2 primary ingestion model = 21 libraries: without 1h-IG50-3, 1h-IGSUB-1, 4h-CTL-IG-2; with 4h-IGSUB-1 (library_28)',nrow(ti)==21 && !any(c('1h-IG50-3','1h-IGSUB-1','4h-CTL-IG-2') %in% ti$isolate) && '4h-IGSUB-1' %in% ti$isolate,paste(nrow(ti),'libraries'))
chk('F3 primary contact model = 22 libraries: without 1h-CTL-CONT-2; with 4h-CONT50-2 (library_17)',nrow(tc)==22 && !'1h-CTL-CONT-2' %in% tc$isolate && '4h-CONT50-2' %in% tc$isolate,paste(nrow(tc),'libraries'))
dF=read.delim(file.path(F$dir,'PRJNA1017469_design_checks.tsv')); pp=dF[dF$model=='primary_if_PENDING_excluded',]
chk('F4 primary design checks: ingestion 21 libraries / 15 residual df / smallest group 3; contact 22 / 16 / 3; rank 6; 7/7 contrasts estimable',
    all(pp$PASS) && identical(as.integer(pp$n_libraries),c(21L,22L)) && identical(as.integer(pp$residual_df),c(15L,16L)) && all(pp$min_group==3) && all(pp$rank==6) && all(pp$contrasts_estimable=='7/7'),
    paste(pp$route,pp$n_libraries,pp$residual_df,pp$min_group,collapse=' '))
say(sprintf('\nSUMMARY: %d PASS, %d FAIL   (scratch: %s)',sum(RES),sum(!RES),normalizePath(SCR)))
if(any(!RES)) {say('FAILED: ',paste(names(RES)[!RES],collapse=', ')); quit(status=1)}
