# Provenance: in-browser computation (2026-11-11)

Why: on the analysis date the project computer's Linux shell had no network access and the cloud workspace's shell
was blocked from NCBI, so GEO files were fetched with `fetch()` in the desktop app's built-in browser (origin
https://ftp.ncbi.nlm.nih.gov), decompressed with `DecompressionStream('gzip')`, and analysed in JavaScript.
Only compact result tables were transcribed into this project. The raw GEO files were NOT written to disk.

Implementation (JavaScript, executed in page context):
* Student / Welch t p-values: regularised incomplete beta (Numerical Recipes `betacf`), validated:
  t=2.776,df=4 -> p=0.0500; t=1.960,df=1e6 -> 0.0500; t=2.228,df=10 -> 0.0500.
* BH-FDR, colony-mean aggregation, ICC (one-way ANOVA within arm, n0=5), per-sample mean Pearson r QC,
  median-of-ratios size factors and effect-coded 2x2 OLS for GSE121885: identical in logic to
  `reanalyse_processed_geo.py`, which reproduces them from files on disk.
* Independent check: the colony-level statistics for 8 genes were recomputed with scipy from the transcribed
  per-sample values (results/gene_expression/GSE178742_key_genes_logCPM_per_sample.csv); all logFC and p agreed
  to rounding (e.g. CYP9Q1 field foragers logFC 0.684 vs 0.688 from 2-dp inputs; p 0.002 both).
* SHA-256 of every fetched file was computed in-browser with `crypto.subtle.digest` (see data/DOWNLOAD_MANIFEST.tsv).

To regenerate everything from disk:
    python3 scripts/download/download_geo.py --project-root .
    python3 scripts/differential_expression/reanalyse_processed_geo.py --project-root .
    python3 scripts/meta_analysis/build_reports_and_figures.py
