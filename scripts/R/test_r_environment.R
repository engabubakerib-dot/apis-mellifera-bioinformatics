## R environment validation (NOT a differential-expression analysis of project data).
## Run from the project folder:  Rscript --vanilla scripts/R/test_r_environment.R
## Uses the project library envs/R-library/<major.minor>. Writes results/R_env/R_validation_<ts>.txt.
##  R1  core packages load (tximport, edgeR, limma, DESeq2) and versions are recorded
##  R2  tximport on REAL kallisto.exe outputs produced by the Windows preflight (scratch folder): every target_id maps,
##      counts are conserved (sum of gene counts = sum of transcript est_counts), ignoreTxVersion=TRUE is shown to be wrong
##  R3  statistical validation on SIMULATED counts (replicate datasets; rules fixed in advance, see the R3 block): null calibration,
##      FDR control and known-answer recovery for the PRIMARY methods (edgeR QL; colony pseudo-bulk edgeR); DESeq2 and bee-level
##      voom + duplicateCorrelation reported as SECONDARY; power grid for interpretation. R_VALIDATION_NSIM (default 20) sets replicates.
##  R4  de_from_kallisto.R imp() skips a study with missing runs (no partial designs) - tested on a scratch copy
SCRIPT_VERSION <- '2026-11-12d'
ts <- format(Sys.time(), '%Y%m%d_%H%M%S'); proj <- normalizePath(getwd(), winslash = '/')
self <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
rl <- file.path(proj, 'envs', 'R-library', paste(R.version$major, sub('\\..*', '', R.version$minor), sep = '.'))
if (dir.exists(rl)) .libPaths(c(rl, .libPaths()))
dir.create('results/R_env', recursive = TRUE, showWarnings = FALSE)
rep <- file.path('results/R_env', paste0('R_validation_', ts, '.txt')); RES <- logical(); RESN <- character()
say <- function(...) { m <- paste0(...); cat(m, '\n'); cat(m, '\n', file = rep, append = TRUE) }
chk <- function(name, ok, detail = '') { ok <- isTRUE(as.logical(ok)); RES <<- c(RES, ok); RESN <<- c(RESN, name); say(if (ok) 'PASS ' else 'FAIL ', name, if (nzchar(detail)) paste0('   [', detail, ']') else '') }
tryk <- function(name, expr) tryCatch(expr, error = function(e) chk(paste(name, ':: error'), FALSE, conditionMessage(e)))
say('R validation ', ts, ' | script version ', SCRIPT_VERSION, ' | ', normalizePath(self, winslash = '/', mustWork = FALSE), ' md5 ', if (!is.na(self) && file.exists(self)) unname(tools::md5sum(self)) else 'n/a')
say('R validation ', ts, ' | ', R.version.string, ' | Rscript ', file.path(R.home('bin'), 'Rscript'), ' | libPaths: ', paste(.libPaths(), collapse = ' ; '))

## R1 ---------------------------------------------------------------------------------------------------------------
for (p in c('tximport', 'rhdf5', 'edgeR', 'limma', 'DESeq2')) tryk(p, {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  chk(paste('R1 loads:', p), TRUE, paste(as.character(packageVersion(p)), 'from', dirname(find.package(p)))) })
for (p in c('GO.db', 'DOSE', 'enrichplot', 'clusterProfiler', 'org.Amellifera.eg.db')) {
  ok <- suppressWarnings(suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE)))
  say(if (ok) 'INFO ' else 'INFO ', 'enrichment package ', p, ': ', if (ok) paste('loads', as.character(packageVersion(p))) else 'NOT available (does not block core RNA-seq analysis)') }
source('scripts/R/tx2gene.R')

## R2 ---------------------------------------------------------------------------------------------------------------
tryk('R2', {
  t2g <- make_tx2gene()
  chk('R2.1 tx2gene built from GFF3', nrow(t2g) > 15000 && !any(duplicated(t2g$tx)), paste(nrow(t2g), 'transcripts,', length(unique(t2g$gene)), 'genes'))
  ab <- c(Sys.glob(file.path('preflight_scratch', '*', 'n1_ena', 'R', 'results', 'rawreads', 'kallisto', '*', '*', 'SRR*', 'abundance.tsv')),   # real ENA reads
         Sys.glob(file.path('preflight_scratch', '*', 'o9_real_kallisto_*', 'R', 'results', 'rawreads', 'kallisto', '*', '*', 'SRR*', 'abundance.tsv')))
  ab <- ab[file.exists(file.path(dirname(ab), '_COMPLETE.json'))]
  if (length(ab) < 2) { say('N/A  R2.2-2.5 need >= 2 COMPLETE kallisto.exe outputs from a passing preflight (found ', length(ab), ')') } else {
    ab <- head(ab, 2); names(ab) <- c('s1', 's2'); say('R2 real kallisto.exe outputs: ', paste(dirname(ab), collapse = ' ; '))
    chk('R2.2 every target_id of the real kallisto.exe outputs maps to tx2gene', isTRUE(check_tx2gene(ab, t2g)))
    tx <- tximport(ab, type = 'kallisto', tx2gene = t2g, countsFromAbundance = 'no', ignoreTxVersion = FALSE)
    raw <- sapply(ab, function(f) sum(read.delim(f)$est_counts))
    chk('R2.3 counts conserved: gene-level sum = transcript est_counts sum (both samples)', all(abs(colSums(tx$counts) - raw) < 1e-6 * pmax(1, raw)), paste(round(colSums(tx$counts), 2), collapse = ' / '))
    chk('R2.4 gene-level matrix has one row per gene, no NA', nrow(tx$counts) == length(unique(t2g$gene)) && !anyNA(tx$counts), paste(nrow(tx$counts), 'genes'))
    ls <- tximport(ab, type = 'kallisto', tx2gene = t2g, countsFromAbundance = 'lengthScaledTPM', ignoreTxVersion = FALSE)
    chk('R2.5 lengthScaledTPM library sizes equal the original library sizes', all(abs(colSums(ls$counts) - colSums(tx$counts)) < 1e-3 * colSums(tx$counts)))
    h5 <- sub('abundance\\.tsv$', 'abundance.h5', ab)
    chk('R2.8 abundance.h5 present next to every abundance.tsv (kallisto output complete)', all(file.exists(h5)))
    txh <- tryCatch(tximport(h5, type = 'kallisto', tx2gene = t2g, countsFromAbundance = 'no', ignoreTxVersion = FALSE), error = function(e) conditionMessage(e))
    chk('R2.9 tximport reads the real kallisto HDF5 files (abundance.h5) via rhdf5', is.list(txh), if (is.character(txh)) txh else 'ok')
    if (is.list(txh)) chk('R2.10 HDF5 and TSV imports agree (counts, lengths; max relative difference < 1e-4)',
        identical(dim(txh$counts), dim(tx$counts)) && max(abs(txh$counts - tx$counts) / pmax(1, tx$counts)) < 1e-4 && max(abs(txh$length - tx$length) / pmax(1, tx$length)) < 1e-4,
        sprintf('max rel diff counts %.2e, lengths %.2e', max(abs(txh$counts - tx$counts) / pmax(1, tx$counts)), max(abs(txh$length - tx$length) / pmax(1, tx$length))))
    txd <- tximport(ab, type = 'kallisto', tx2gene = t2g, countsFromAbundance = 'no', ignoreTxVersion = FALSE, dropInfReps = TRUE)
    chk('R2.11 dropInfReps=TRUE (used by the DE script; kallisto run with 0 bootstraps) gives identical counts to the default import', isTRUE(all.equal(txd$counts, tx$counts)) && is.null(tx$infReps))
    bad <- tryCatch(suppressWarnings({ b <- tximport(ab, type = 'kallisto', tx2gene = t2g, countsFromAbundance = 'no', ignoreTxVersion = TRUE)
                      if (nrow(b$counts) < 0.5 * nrow(tx$counts)) paste('silently collapsed to', nrow(b$counts), 'genes') else 'ok' }),
                    error = function(e) paste('error:', conditionMessage(e)))
    dds0 <- tryCatch(DESeq2::DESeqDataSetFromTximport(tx, data.frame(x = factor(c('a', 'b')), row.names = colnames(tx$counts)), ~1), error = function(e) conditionMessage(e))
    chk('R2.7 DESeqDataSetFromTximport accepts the real tximport object (length offsets present)', inherits(dds0, 'DESeqDataSet') && 'avgTxLength' %in% names(SummarizedExperiment::assays(dds0)), if (is.character(dds0)) dds0 else '')
    chk('R2.6 control: ignoreTxVersion=TRUE breaks these IDs (error or collapse) - why the script uses FALSE', bad != 'ok', substr(bad, 1, 120))
  } })

## R3 ---------------------------------------------------------------------------------------------------------------
## Statistical validation by simulation. Acceptance rules were fixed BEFORE running (2026-11-12) and are not tuned to results:
##  PRIMARY methods (gate): edgeR QL for pooled 3 v 3 designs; colony pseudo-bulk edgeR QL for bee-level apiary data.
##   null calibration : mean type-I error (p<0.05 among tested null genes) <= 0.06 AND share of null datasets with any FDR<0.05 call <= 0.15
##   FDR              : mean false discovery proportion (FDP) over replicate datasets <= 0.10 (nominal FDR 0.05)
##   known answer     : moderate dispersion (phi 0.05), |log2FC| = 2 on well-expressed genes -> power >= 0.80, all significant
##                      planted genes with the correct sign, median |estimated - true log2FC| < 0.30   (pooled design only)
##  SECONDARY methods (reported, not gating): DESeq2; bee-level limma-voom + duplicateCorrelation. A secondary method that fails
##   calibration is used for sensitivity only. Bee-level edgeR ignoring colony (pseudoreplication) is shown as a negative control.
##  POWER grid (information only): edgeR QL 3 v 3 by dispersion x effect size - to interpret non-significant results honestly.
NSIM <- as.integer(Sys.getenv('R_VALIDATION_NSIM', '20')); G <- 3000; say('R3 simulation: ', NSIM, ' replicate datasets per scenario, ', G, ' genes')
STAT <- list(); addstat <- function(...) STAT[[length(STAT) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
info <- function(...) say('INFO ', ...)
mk_mu0 <- function() exp(rnorm(G, log(80), 1.4))                      # mean counts: median 80, wide spread
mk_phi <- function(phi) phi * exp(rnorm(G, 0, 0.4))                   # gene-wise dispersion around phi
gid <- function(rn) as.integer(sub('^g', '', rn))
sim_pooled <- function(n, phi, lfc, mu0, phig) {                        # n per arm; independent libraries (pooled samples)
  arm <- rep(0:1, each = n); lib <- runif(2 * n, 0.7, 1.3)
  y <- sapply(seq_along(arm), function(j) rnbinom(G, mu = mu0 * lib[j] * 2^(lfc * arm[j]), size = 1 / phig))
  storage.mode(y) <- 'integer'; rownames(y) <- paste0('g', 1:G); colnames(y) <- paste0('s', seq_along(arm)); list(y = y, grp = factor(c('C', 'T')[arm + 1])) }
sim_colony <- function(phi, lfc, mu0, phig, sdc = 0.3, ncol = 3, nbee = 5) {   # bees nested in colonies; colony is the unit
  carm <- rep(0:1, each = ncol); col <- rep(seq_len(2 * ncol), each = nbee); ce <- matrix(rnorm(G * 2 * ncol, 0, sdc), G)
  lib <- runif(length(col), 0.7, 1.3)
  y <- sapply(seq_along(col), function(j) rnbinom(G, mu = mu0 * lib[j] * exp(ce[, col[j]]) * 2^(lfc * carm[col[j]]), size = 1 / phig))
  storage.mode(y) <- 'integer'; rownames(y) <- paste0('g', 1:G); colnames(y) <- paste0('b', seq_along(col)); list(y = y, col = col, trt = factor(c('C', 'T')[carm[col] + 1]), carm = factor(c('C', 'T')[carm + 1])) }
res_full <- function(idx, p, fdr, lfc) { o <- data.frame(p = rep(NA_real_, G), fdr = rep(1, G), lfc = rep(NA_real_, G)); o$p[idx] <- p; o$fdr[idx] <- fdr; o$lfc[idx] <- lfc; o }
f_edger <- function(y, grp) {
  d <- edgeR::DGEList(y); d <- d[edgeR::filterByExpr(d, group = grp), , keep.lib.sizes = FALSE]; d <- edgeR::calcNormFactors(d)
  des <- model.matrix(~grp); d <- edgeR::estimateDisp(d, des)
  tt <- edgeR::topTags(edgeR::glmQLFTest(edgeR::glmQLFit(d, des, robust = TRUE), coef = 2), n = Inf, sort.by = 'none')$table
  res_full(gid(rownames(tt)), tt$PValue, tt$FDR, tt$logFC) }
f_deseq <- function(y, grp) {
  dds <- DESeq2::DESeqDataSetFromMatrix(y, data.frame(grp = grp, row.names = colnames(y)), ~grp)
  r <- DESeq2::results(DESeq2::DESeq(dds, quiet = TRUE), name = 'grp_T_vs_C', alpha = 0.05)
  ok <- !is.na(r$pvalue); res_full(gid(rownames(r))[ok], r$pvalue[ok], ifelse(is.na(r$padj[ok]), 1, r$padj[ok]), r$log2FoldChange[ok]) }
f_pseudobulk <- function(s) { pb <- sapply(split(seq_along(s$col), s$col), function(z) rowSums(s$y[, z, drop = FALSE])); rownames(pb) <- rownames(s$y); f_edger(pb, s$carm) }
f_voomdc <- function(s) {
  des <- model.matrix(~s$trt); d <- edgeR::DGEList(s$y); d <- edgeR::calcNormFactors(d[edgeR::filterByExpr(d, group = s$trt), , keep.lib.sizes = FALSE])
  v <- limma::voom(d, des); dc <- limma::duplicateCorrelation(v, des, block = s$col)
  v <- limma::voom(d, des, block = s$col, correlation = dc$consensus); dc <- limma::duplicateCorrelation(v, des, block = s$col)
  tt <- limma::topTable(limma::eBayes(limma::lmFit(v, des, block = s$col, correlation = dc$consensus)), coef = 2, n = Inf, sort.by = 'none')
  o <- res_full(gid(rownames(tt)), tt$P.Value, tt$adj.P.Val, tt$logFC); attr(o, 'rho') <- dc$consensus; o }
f_beeindep <- function(s) f_edger(s$y, s$trt)
metr <- function(o, truth) { disc <- which(o$fdr < 0.05); nullp <- o$p[!truth & !is.na(o$p)]
  c(typeI = if (length(nullp)) mean(nullp < 0.05) else NA, any = as.numeric(length(disc) > 0), fdp = if (length(disc)) mean(!truth[disc]) else 0,
    power = if (any(truth)) mean(o$fdr[truth] < 0.05) else NA, ndisc = length(disc)) }
runsc <- function(label, nsim, gen, fits) {   # gen(): list(data, truth); fits: named list of functions(data)
  m <- lapply(fits, function(.) NULL); t0 <- Sys.time()
  for (r in seq_len(nsim)) { g <- gen(); for (k in names(fits)) m[[k]] <- rbind(m[[k]], metr(fits[[k]](g$data), g$truth)) }
  say(sprintf('  scenario %-48s %d datasets, %.0f s', label, nsim, as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
  lapply(m, function(x) colMeans(x, na.rm = TRUE)) }
tryk('R3', {
  set.seed(20261112)
  ## S1 null, pooled 3 v 3, phi 0.10
  s1 <- runsc('S1 null, pooled 3 v 3, phi 0.10', NSIM, function() { mu0 <- mk_mu0(); d <- sim_pooled(3, 0.1, rep(0, G), mu0, mk_phi(0.1)); list(data = d, truth = rep(FALSE, G)) },
              list(edgeR_QL = function(d) f_edger(d$y, d$grp), DESeq2 = function(d) f_deseq(d$y, d$grp)))
  ## S2 FDR, pooled 3 v 3, phi 0.10, 10% genes DE with |log2FC| uniform 0.5-2
  s2 <- runsc('S2 FDR, pooled 3 v 3, phi 0.10, 10% DE', NSIM, function() { mu0 <- mk_mu0(); tr <- seq_len(G) %in% sample(G, 300)
              lfc <- ifelse(tr, sample(c(-1, 1), G, TRUE) * runif(G, 0.5, 2), 0); list(data = sim_pooled(3, 0.1, lfc, mu0, mk_phi(0.1)), truth = tr) },
              list(edgeR_QL = function(d) f_edger(d$y, d$grp), DESeq2 = function(d) f_deseq(d$y, d$grp)))
  ## S3 known answer, pooled 3 v 3, phi 0.05, |log2FC| = 2 on genes with mean >= 50
  KA <- list(); k3 <- function(fun, nm) function(d) { o <- fun(d); KA[[nm]] <<- rbind(KA[[nm]], data.frame(est = o$lfc, true = attr(d, 'lfc'), sig = o$fdr < 0.05, de = attr(d, 'lfc') != 0)); o }
  s3 <- runsc('S3 known answer, pooled 3 v 3, phi 0.05, |LFC| 2', NSIM, function() { mu0 <- mk_mu0(); el <- which(mu0 >= 50); tr <- seq_len(G) %in% sample(el, 200)
              lfc <- ifelse(tr, sample(c(-2, 2), G, TRUE), 0); d <- sim_pooled(3, 0.05, lfc, mu0, mk_phi(0.05)); attr(d, 'lfc') <- lfc; list(data = d, truth = tr) },
              list(edgeR_QL = k3(function(d) f_edger(d$y, d$grp), 'edgeR_QL'), DESeq2 = k3(function(d) f_deseq(d$y, d$grp), 'DESeq2')))
  ## S4 colony structure: 3 v 3 colonies x 5 bees, colony effect sd 0.3 (log scale), bee phi 0.10 - null and 10% DE (|log2FC| 0.5-2)
  RHO <- c(); vdc <- function(s) { o <- f_voomdc(s); RHO <<- c(RHO, attr(o, 'rho')); o }
  cf <- list(pseudobulk_edgeR = function(s) f_pseudobulk(s), voom_dupCor = vdc, bees_independent_edgeR = function(s) f_beeindep(s))
  s4n <- runsc('S4 null, colonies 3 v 3 x 5 bees', NSIM, function() { mu0 <- mk_mu0(); list(data = sim_colony(0.1, rep(0, G), mu0, mk_phi(0.1)), truth = rep(FALSE, G)) }, cf)
  s4a <- runsc('S4 FDR, colonies 3 v 3 x 5 bees, 10% DE', NSIM, function() { mu0 <- mk_mu0(); tr <- seq_len(G) %in% sample(G, 300)
               lfc <- ifelse(tr, sample(c(-1, 1), G, TRUE) * runif(G, 0.5, 2), 0); list(data = sim_colony(0.1, lfc, mu0, mk_phi(0.1)), truth = tr) }, cf)
  for (k in names(s1)) addstat(scenario = 'S1 null pooled', method = k, typeI = s1[[k]]['typeI'], any_FDR_call = s1[[k]]['any'], FDP = NA, power = NA)
  for (k in names(s2)) addstat(scenario = 'S2 FDR pooled', method = k, typeI = s2[[k]]['typeI'], any_FDR_call = NA, FDP = s2[[k]]['fdp'], power = s2[[k]]['power'])
  for (k in names(s3)) addstat(scenario = 'S3 known answer pooled', method = k, typeI = s3[[k]]['typeI'], any_FDR_call = NA, FDP = s3[[k]]['fdp'], power = s3[[k]]['power'])
  for (k in names(s4n)) addstat(scenario = 'S4 null colony', method = k, typeI = s4n[[k]]['typeI'], any_FDR_call = s4n[[k]]['any'], FDP = NA, power = NA)
  for (k in names(s4a)) addstat(scenario = 'S4 FDR colony', method = k, typeI = s4a[[k]]['typeI'], any_FDR_call = NA, FDP = s4a[[k]]['fdp'], power = s4a[[k]]['power'])
  ## gates (primary methods)
  e1 <- s1$edgeR_QL; chk('R3.1 edgeR QL null calibration (pooled 3 v 3): mean type-I <= 0.06 and datasets with any call <= 0.15', e1['typeI'] <= 0.06 && e1['any'] <= 0.15, sprintf('type-I %.4f, any-call %.2f', e1['typeI'], e1['any']))
  e2 <- s2$edgeR_QL; chk('R3.2 edgeR QL FDR (pooled 3 v 3, 10% DE): mean FDP <= 0.10', e2['fdp'] <= 0.10, sprintf('mean FDP %.3f, power %.2f', e2['fdp'], e2['power']))
  ka <- KA$edgeR_QL; kd <- ka[ka$de & ka$sig, ]
  chk('R3.3 edgeR QL known answer: power >= 0.80, correct sign for all significant planted genes, median |error| < 0.30',
      s3$edgeR_QL['power'] >= 0.8 && all(sign(kd$est) == sign(kd$true)) && median(abs(kd$est - kd$true)) < 0.3,
      sprintf('power %.2f, wrong-sign %d, median |log2FC error| %.3f', s3$edgeR_QL['power'], sum(sign(kd$est) != sign(kd$true)), median(abs(kd$est - kd$true))))
  p4 <- s4n$pseudobulk_edgeR; chk('R3.4 colony pseudo-bulk edgeR null calibration: mean type-I <= 0.06 and any-call <= 0.15', p4['typeI'] <= 0.06 && p4['any'] <= 0.15, sprintf('type-I %.4f, any-call %.2f', p4['typeI'], p4['any']))
  p5 <- s4a$pseudobulk_edgeR; chk('R3.5 colony pseudo-bulk edgeR FDR: mean FDP <= 0.10', p5['fdp'] <= 0.10, sprintf('mean FDP %.3f, power %.2f', p5['fdp'], p5['power']))
  ## secondary methods: reported with verdict, not gating
  kd2 <- KA$DESeq2[KA$DESeq2$de & KA$DESeq2$sig, ]
  dsq_ok <- s1$DESeq2['typeI'] <= 0.06 && s1$DESeq2['any'] <= 0.15 && s2$DESeq2['fdp'] <= 0.10
  info(sprintf('R3.6 SECONDARY DESeq2 (pooled 3 v 3): null type-I %.4f, any-call %.2f; FDR mean FDP %.3f, power %.2f; known answer power %.2f, wrong-sign %d -> %s',
       s1$DESeq2['typeI'], s1$DESeq2['any'], s2$DESeq2['fdp'], s2$DESeq2['power'], s3$DESeq2['power'], sum(sign(kd2$est) != sign(kd2$true)),
       if (dsq_ok) 'calibrated: may be reported alongside edgeR' else 'NOT calibrated in this design: sensitivity only'))
  vok <- s4n$voom_dupCor['typeI'] <= 0.06 && s4n$voom_dupCor['any'] <= 0.15 && s4a$voom_dupCor['fdp'] <= 0.10
  info(sprintf('R3.7 SECONDARY bee-level voom + duplicateCorrelation: null type-I %.4f, any-call %.2f; FDR mean FDP %.3f, power %.2f; consensus intra-colony r median %.3f -> %s',
       s4n$voom_dupCor['typeI'], s4n$voom_dupCor['any'], s4a$voom_dupCor['fdp'], s4a$voom_dupCor['power'], median(RHO), if (vok) 'calibrated' else 'NOT calibrated with 6 colonies: sensitivity only'))
  info(sprintf('R3.8 NEGATIVE CONTROL bees treated as independent (pseudoreplication; never used): null type-I %.4f, any-call %.2f, FDR mean FDP %.3f',
       s4n$bees_independent_edgeR['typeI'], s4n$bees_independent_edgeR['any'], s4a$bees_independent_edgeR['fdp']))
  ## power grid (information only)
  for (phi in c(0.05, 0.10, 0.25)) for (L in c(1, 1.5, 2)) {
    pw <- runsc(sprintf('P power edgeR QL 3 v 3, phi %.2f, |LFC| %.1f', phi, L), max(3, NSIM %/% 4), function() { mu0 <- mk_mu0(); el <- which(mu0 >= 50); tr <- seq_len(G) %in% sample(el, 200)
                lfc <- ifelse(tr, sample(c(-L, L), G, TRUE), 0); list(data = sim_pooled(3, phi, lfc, mu0, mk_phi(phi)), truth = tr) }, list(edgeR_QL = function(d) f_edger(d$y, d$grp)))
    addstat(scenario = sprintf('POWER phi %.2f |LFC| %.1f (mean>=50)', phi, L), method = 'edgeR_QL', typeI = NA, any_FDR_call = NA, FDP = pw$edgeR_QL['fdp'], power = pw$edgeR_QL['power'])
    info(sprintf('power grid: phi %.2f (BCV %.2f), |log2FC| %.1f -> power %.2f, mean FDP %.3f', phi, sqrt(phi), L, pw$edgeR_QL['power'], pw$edgeR_QL['fdp'])) }
  st <- do.call(rbind, STAT); rownames(st) <- NULL
  write.table(st, file.path('results/R_env', paste0('stat_validation_', ts, '.tsv')), sep = '\t', quote = FALSE, row.names = FALSE)
  say('statistical validation table: results/R_env/stat_validation_', ts, '.tsv') })

## R4 ---------------------------------------------------------------------------------------------------------------
tryk('R4', {
  src <- readLines('scripts/R/de_from_kallisto.R', warn = FALSE)
  chk('R4.1 DE script skips any study with missing runs (no partial designs)', any(grepl('study SKIPPED', src, fixed = TRUE)) && any(grepl('setdiff(expected,names(f))', src, fixed = TRUE)))
  chk('R4.2 DE script uses ignoreTxVersion=FALSE only + check_tx2gene()', !any(grepl('ignoreTxVersion=TRUE', src, fixed = TRUE)) && any(grepl('check_tx2gene(', src, fixed = TRUE)))
  chk('R4.3 DE script stops on colonies with mixed treatment labels', any(grepl('mixed treatment labels', src, fixed = TRUE))) })

say('')
say('SUMMARY: ', sum(RES), ' PASS, ', sum(!RES), ' FAIL  -> ', if (length(RES) && all(RES)) 'R ENVIRONMENT VALIDATED (core + primary statistical methods)' else 'NOT VALIDATED')
if (any(!RES)) for (n in RESN[!RES]) say('  FAILED: ', n)
write.table(data.frame(test = RESN, result = ifelse(RES, 'PASS', 'FAIL')), file.path('results/R_env', paste0('R_validation_', ts, '.tsv')), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path('results/R_env', paste0('sessionInfo_validation_', ts, '.txt')))
quit(status = if (length(RES) && all(RES)) 0 else 1)
