## Transcript -> gene table for kallisto outputs quantified against data/raw/BeeBase/amel_OGSv3.2_cds.fa.gz.
## The CDS FASTA (and therefore kallisto target_id) uses headers like "gnl|Amel_4.5|GB48241-RA", whereas the OGS v3.2
## GFF3 mRNA IDs are "GB48241-RA" (Parent = gene "GB48241"). The prefix is therefore added to the GFF3 IDs.
## tximport must be called with ignoreTxVersion = FALSE: with TRUE, tximport strips everything after the first "." and
## "gnl|Amel_4.5|GB48241-RA" would collapse to "gnl|Amel_4" for every transcript.
make_tx2gene <- function(gff = 'data/raw/BeeBase/amel_OGSv3.2.gff3.gz', prefix = 'gnl|Amel_4.5|') {
  g <- read.delim(gzfile(gff), header = FALSE, comment.char = '#', quote = '', colClasses = 'character')
  m <- g[g$V3 == 'mRNA', 9]
  tx <- sub('.*ID=([^;]+).*', '\\1', m); gene <- sub('.*Parent=([^;]+).*', '\\1', m)
  if (any(duplicated(tx))) stop('duplicated mRNA IDs in GFF3')
  data.frame(tx = paste0(prefix, tx), gene = gene, stringsAsFactors = FALSE)
}
## Stops unless EVERY target_id of every abundance file is in tx2gene (no silent loss of transcripts).
check_tx2gene <- function(files, tx2gene) {
  for (f in files) {
    ids <- read.delim(f, colClasses = c('character', rep('NULL', 4)))[[1]]
    miss <- setdiff(ids, tx2gene$tx)
    if (length(miss)) stop(sprintf('%s: %d of %d target_ids not in tx2gene (e.g. %s)', f, length(miss), length(ids), paste(head(miss, 3), collapse = ', ')))
    if (any(duplicated(ids))) stop(f, ': duplicated target_ids')
  }
  invisible(TRUE)
}
