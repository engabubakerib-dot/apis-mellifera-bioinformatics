# Project R package setup - existing R only (no R installation, no admin rights, never writes system or personal libraries).
# Run from the project folder:  Rscript --vanilla scripts/R/setup_r_packages.R [--inspect-only] [--force-python]
#  1 INSPECTION   R, .libPaths(), what is installed where, packages used by the project scripts
#  2 NETWORK      exact error of every repository URL x download method (libcurl, wininet), proxy/TLS details
#  3 INSTALL      route A (R can reach the repositories): BiocManager::install(update = FALSE) into the project library
#                 route B (R cannot): scripts/R/fetch_r_packages.py (Python, which reaches the same servers) downloads the
#                 official package files, verifies MD5 against the repository index; R installs them from local files
#  4 VERIFY       requireNamespace() + library() of every package in a FRESH R session; exact errors recorded
# Outputs (results/R_env/): R_setup_report_<ts>.txt, network_diagnosis_<ts>.tsv, install_<ts>.log, R_packages_<ts>.tsv,
#                           sessionInfo_<ts>.txt; CORE_READY.txt only if all core packages load. Nothing is deleted.
SCRIPT_VERSION <- "2026-11-12f"
ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); args <- commandArgs(trailingOnly = TRUE)
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
proj <- normalizePath(getwd(), winslash = "/")
rminor <- paste(R.version$major, sub("\\..*", "", R.version$minor), sep = ".")
projlib <- file.path(proj, "envs", "R-library", rminor)
outdir <- file.path(proj, "results", "R_env"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rep <- file.path(outdir, paste0("R_setup_report_", ts, ".txt")); instlog <- file.path(outdir, paste0("install_", ts, ".log"))
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = rep, append = TRUE) }
ilog <- function(...) cat(paste0(...), "\n", file = instlog, append = TRUE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
options(timeout = max(7200, getOption("timeout")), warn = 1)
CORE <- c("BiocManager", "tximport", "rhdf5", "edgeR", "limma", "DESeq2")   # rhdf5: tximport needs it to read kallisto abundance.h5 (inferential replicates / HDF5 import)
ENRICH <- c("GO.db", "DOSE", "enrichplot", "clusterProfiler")
OPTIONAL <- "org.Amellifera.eg.db"
BIOC_FOR_R <- c("4.6" = "3.23", "4.5" = "3.22")   # official: bioconductor.org/about/release-announcements (3.23: 2026-04-29, R 4.6)
bioc <- unname(BIOC_FOR_R[rminor]); if (is.na(bioc)) stop("no Bioconductor release recorded for R ", rminor)
REPOS <- c(BioCsoft = sprintf("https://bioconductor.org/packages/%s/bioc", bioc),
           BioCann = sprintf("https://bioconductor.org/packages/%s/data/annotation", bioc),
           BioCexp = sprintf("https://bioconductor.org/packages/%s/data/experiment", bioc),
           CRAN = "https://cloud.r-project.org")

## ---------------------------------------------------------------------------------------------------- 1 INSPECTION
say("== setup_r_packages.R version ", SCRIPT_VERSION, " | ", normalizePath(self, winslash = "/", mustWork = FALSE), " | md5 ", if (!is.na(self) && file.exists(self)) unname(tools::md5sum(self)) else "n/a")
say("== 1. INSPECTION (before any change) ", format(Sys.time()))
say("R version: ", R.version.string, "  | R home: ", normalizePath(R.home(), winslash = "/"), "  | Rscript: ", normalizePath(rscript, winslash = "/", mustWork = FALSE))
say("Bioconductor release for R ", rminor, ": ", bioc, "  | repositories: ", paste(names(REPOS), REPOS, sep = "=", collapse = "  "))
lib0 <- .libPaths()
writable <- function(p) dir.exists(p) && file.access(p, 2) == 0
for (p in lib0) say(".libPaths() at start: ", p, "  writable=", writable(p), "  packages=", length(list.dirs(p, recursive = FALSE)))
say("Project library: ", projlib, "  exists before run=", dir.exists(projlib))
used <- unique(unlist(lapply(list.files(file.path(proj, "scripts", "R"), pattern = "\\.R$", full.names = TRUE), function(f) {
  x <- readLines(f, warn = FALSE)
  c(unlist(regmatches(x, gregexpr("(?<=library\\()[A-Za-z0-9.]+(?=\\))", x, perl = TRUE))), unlist(regmatches(x, gregexpr("[A-Za-z][A-Za-z0-9.]+(?=::)", x, perl = TRUE)))) })))
say("Packages referenced by project R scripts: ", paste(sort(setdiff(used, c("base", "utils", "stats", "methods", "tools", "BiocManager"))), collapse = ", "))
where <- function(pkg, libs = .libPaths()) {
  for (l in libs) { d <- file.path(l, pkg); if (file.exists(file.path(d, "DESCRIPTION")))
    return(list(lib = l, version = as.character(read.dcf(file.path(d, "DESCRIPTION"), "Version")))) }
  NULL }
for (p in c(CORE, ENRICH, OPTIONAL)) { w <- where(p, lib0); say(sprintf("  %-22s %s", p, if (is.null(w)) "NOT INSTALLED" else paste(w$version, "in", w$lib))) }
if ("--inspect-only" %in% args) { say("== INSPECT-ONLY: nothing installed. Report: ", rep); quit(status = 0) }

## ---------------------------------------------------------------------------------------------------- 2 NETWORK
say("== 2. NETWORK DIAGNOSIS (exact errors; no settings changed yet)")
lv <- tryCatch(libcurlVersion(), error = function(e) NA)
say("libcurl: ", as.character(lv), "  ssl backend: ", if (!is.null(attr(lv, "ssl_version"))) attr(lv, "ssl_version") else "?",
    "  | capabilities(libcurl)=", capabilities("libcurl"), "  | default download.file.method=", format(getOption("download.file.method")), "  url.method=", format(getOption("url.method")))
mask <- function(x) sub("//[^@/]*@", "//***@", x)
for (v in c("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY", "no_proxy", "R_LIBCURL_SSL_REVOKE_BEST_EFFORT", "CURL_CA_BUNDLE", "SSL_CERT_FILE"))
  if (nzchar(Sys.getenv(v))) say("env ", v, " = ", mask(Sys.getenv(v)))
if (requireNamespace("curl", quietly = TRUE)) {
  pr <- tryCatch(curl::ie_get_proxy_for_url("https://bioconductor.org"), error = function(e) paste("n/a:", conditionMessage(e)))
  say("Windows system proxy for https://bioconductor.org (curl::ie_get_proxy_for_url): ", if (is.null(pr)) "none (direct)" else mask(paste(pr, collapse = " ")))
}
URLS <- c(config = "https://bioconductor.org/config.yaml",
          bioc_index = sprintf("%s/bin/windows/contrib/%s/PACKAGES", REPOS[["BioCsoft"]], rminor),
          annotation_index = sprintf("%s/src/contrib/PACKAGES", REPOS[["BioCann"]]),
          cran_index = sprintf("%s/bin/windows/contrib/%s/PACKAGES", REPOS[["CRAN"]], rminor))
probe <- function(url, method) {
  tf <- tempfile(); msgs <- character(); t0 <- Sys.time()
  ok <- tryCatch(withCallingHandlers(download.file(url, tf, method = method, quiet = TRUE, mode = "wb") == 0,
                   warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") }),
                 error = function(e) { msgs <<- c(msgs, conditionMessage(e)); FALSE })
  size <- if (file.exists(tf)) file.size(tf) else 0
  data.frame(url = url, method = method, revoke_best_effort = Sys.getenv("R_LIBCURL_SSL_REVOKE_BEST_EFFORT"), ok = isTRUE(ok) && size > 0,
             bytes = size, seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
             error = gsub("[\t\r\n]+", " ", paste(unique(msgs), collapse = " | ")), stringsAsFactors = FALSE) }
methods <- c("libcurl", if (.Platform$OS.type == "windows") "wininet")
diag <- do.call(rbind, lapply(methods, function(m) do.call(rbind, lapply(URLS, probe, method = m))))
revoked <- any(grepl("revocation|CRYPT_E_|0x8009201[23]|SEC_E_|schannel", diag$error, ignore.case = TRUE))
if (revoked && !any(diag$ok[diag$method == "libcurl"])) {
  say("libcurl failures mention certificate-revocation checking -> repeating libcurl probes with R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE (certificates are still verified; only an unreachable revocation server is tolerated)")
  Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT = "TRUE"); diag <- rbind(diag, do.call(rbind, lapply(URLS, probe, method = "libcurl")))
}
write.table(diag, file.path(outdir, paste0("network_diagnosis_", ts, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in seq_len(nrow(diag))) say(sprintf("  %-4s %-8s revokeBE=%-4s %-55s %s", if (diag$ok[i]) "OK" else "FAIL", diag$method[i], diag$revoke_best_effort[i],
                                         sub("^https://", "", diag$url[i]), if (diag$ok[i]) paste(diag$bytes[i], "bytes") else diag$error[i]))
route <- "B"; chosen <- NA; chosen_rbe <- ""
combos <- unique(diag[, c("method", "revoke_best_effort")])
combos <- combos[order(combos$revoke_best_effort != "", combos$method != "libcurl"), , drop = FALSE]   # prefer default TLS checks, then libcurl
for (k in seq_len(nrow(combos))) {
  d <- diag[diag$method == combos$method[k] & diag$revoke_best_effort == combos$revoke_best_effort[k], ]
  if (nrow(d) && all(d$ok[d$url != URLS[["config"]]])) { route <- "A"; chosen <- combos$method[k]; chosen_rbe <- combos$revoke_best_effort[k]; break }
}
if (nzchar(chosen_rbe)) Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT = chosen_rbe) else Sys.unsetenv("R_LIBCURL_SSL_REVOKE_BEST_EFFORT")
## Official Bioconductor mirrors (https rows of https://bioconductor.org/BioC_mirrors.csv, 2026-11-12). Used only when
## bioconductor.org itself is unreachable (observed on this network: TCP reset during the TLS handshake for bioconductor.org
## while CRAN works). BiocManager supports mirrors via options(BioC_mirror = ...).
MIRRORS <- c("https://bioconductor.statistik.tu-dortmund.de", "https://ftp.gwdg.de/pub/misc/bioconductor", "https://bioconductor.uib.no",
             "https://mirrors.dotsrc.org/bioconductor", "https://mirror.accum.se/mirror/bioconductor.org", "https://bioconductor.unipi.it",
             "https://bioconductor.posit.co", "https://bioconductor.riken.jp", "https://free.nchc.org.tw/bioconductor",
             "https://mirror.aarnet.edu.au/pub/bioconductor", "https://mirrors.tuna.tsinghua.edu.cn/bioconductor", "https://mirrors.ustc.edu.cn/bioc",
             "https://mirrors.nju.edu.cn/bioconductor", "https://mirrors.zju.edu.cn/bioconductor", "https://mirrors.westlake.edu.cn/bioconductor")
mirror <- NA; reachable <- character()
if (route == "B") {
  cran_ok <- diag$ok[diag$url == URLS[["cran_index"]] & diag$method == "libcurl"]
  if (any(cran_ok)) {
    say("bioconductor.org unreachable but CRAN reachable -> probing official Bioconductor mirrors (libcurl)")
    for (mb in MIRRORS) {
      u <- c(bioc = sprintf("%s/packages/%s/bioc/bin/windows/contrib/%s/PACKAGES", mb, bioc, rminor), ann = sprintf("%s/packages/%s/data/annotation/src/contrib/PACKAGES", mb, bioc))
      d <- do.call(rbind, lapply(u, probe, method = "libcurl")); diag <- rbind(diag, d)
      say(sprintf("  %-4s mirror %-55s %s", if (all(d$ok)) "OK" else "FAIL", mb, if (all(d$ok)) paste(d$bytes, collapse = "+") else paste(unique(d$error[!d$ok]), collapse = " | ")))
      if (all(d$ok)) reachable <- c(reachable, mb)
      if (length(reachable) >= 2) break
    }
    write.table(diag, file.path(outdir, paste0("network_diagnosis_", ts, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    if (length(reachable)) {
      mirror <- reachable[1]; route <- "A"; chosen <- "libcurl"
      options(BioC_mirror = mirror)
      say("ROUTE A via official mirror ", mirror, " (options(BioC_mirror)); reachable mirrors: ", paste(reachable, collapse = ", "))
    } else say("no official mirror reachable either")
  }
}
if ("--force-python" %in% args) route <- "B"
if (route == "A") {
  options(download.file.method = chosen, url.method = chosen)
  cfg_ok <- is.na(mirror) && any(diag$ok[diag$url == URLS[["config"]] & diag$method == chosen & diag$revoke_best_effort == chosen_rbe])
  if (!cfg_ok) { options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = FALSE); say("config.yaml unreachable but package indexes reachable -> BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=FALSE (BiocManager uses its built-in R ", rminor, " -> Bioconductor ", bioc, " map)") }
  say("ROUTE A: R reaches the repositories with method '", chosen, "'", if (nzchar(Sys.getenv("R_LIBCURL_SSL_REVOKE_BEST_EFFORT"))) " + revocation best-effort" else "")
} else say("ROUTE B: R cannot reach the package repositories (see errors above) -> verified download with Python, local install by R")

## ---------------------------------------------------------------------------------------------------- 3 INSTALL
say("== 3. INSTALL into ", projlib, " (system and personal libraries are never written)")
dir.create(projlib, recursive = TRUE, showWarnings = FALSE)
if (!writable(projlib)) { say("FATAL: project library not writable: ", projlib); quit(status = 2) }
.libPaths(c(projlib, lib0)); LIBS <- .libPaths(); for (p in LIBS) say("search path: ", p)
cap <- function(expr) { msgs <- character()
  withCallingHandlers(tryCatch(expr, error = function(e) { msgs <<- c(msgs, paste("ERROR:", conditionMessage(e))); NULL }),
    warning = function(w) { msgs <<- c(msgs, paste("WARNING:", conditionMessage(w))); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, sub("\n$", "", conditionMessage(m))); invokeRestart("muffleMessage") })
  ilog(paste(msgs, collapse = "\n")); msgs }
test_load <- function(p) {
  f <- tempfile(fileext = ".R")
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(LIBS), collapse = "")),
               sprintf("ok <- requireNamespace('%s', quietly=TRUE); cat('REQUIRENAMESPACE', ok, '\\n')", p),
               sprintf("suppressPackageStartupMessages(library('%s', character.only=TRUE))", p),
               sprintf("cat('LOADED', as.character(packageVersion('%s')), normalizePath(dirname(find.package('%s')), winslash='/'), '\\n')", p, p)), f)
  out <- suppressWarnings(system2(rscript, c("--vanilla", shQuote(f, type = if (.Platform$OS.type == "windows") "cmd" else "sh")), stdout = TRUE, stderr = TRUE))
  list(ok = any(grepl("^LOADED ", out)) && is.null(attr(out, "status")), out = out) }
targets <- c(CORE, ENRICH)
todo <- targets[!vapply(targets, function(p) test_load(p)$ok, logical(1))]
say("Packages that do not load before installation: ", if (length(todo)) paste(todo, collapse = ", ") else "none")
kodir <- file.path(outdir, paste0("build_outputs_", ts)); dir.create(kodir, showWarnings = FALSE)   # R CMD INSTALL output of source builds
transient <- "SSL connect error|cannot open URL|unable to access index|Timeout|timed out|Could not resolve|Failed to connect"
if (length(todo) && route == "A") {
  if (is.null(where("BiocManager"))) { say("installing BiocManager from CRAN"); cap(install.packages("BiocManager", lib = projlib, repos = REPOS[["CRAN"]])) }
  library(BiocManager, lib.loc = where("BiocManager")$lib)
  say("BiocManager ", as.character(packageVersion("BiocManager")), "; BiocManager::version() = ", tryCatch(as.character(BiocManager::version()), error = function(e) conditionMessage(e)))
  say("BiocManager::repositories(): ", tryCatch(paste(names(BiocManager::repositories()), BiocManager::repositories(), sep = "=", collapse = "  "), error = function(e) conditionMessage(e)))
  depfixed <- character()
  for (p in setdiff(todo, "BiocManager")) for (attempt in 1:4) {
    say("installing ", p, " (BiocManager, update=FALSE, attempt ", attempt, ") ", format(Sys.time(), "%H:%M:%S")); ilog("######## ", p, " attempt ", attempt, " ", format(Sys.time()))
    m <- cap(BiocManager::install(p, lib = projlib, update = FALSE, ask = FALSE, keep_outputs = kodir,
                                  force = p %in% names(Filter(Negate(is.null), lapply(setNames(nm = p), where)))))
    bad <- grep("^ERROR|non-zero exit|failed|cannot open URL|not available|unable to", m, value = TRUE, ignore.case = TRUE)
    bad <- bad[!grepl("data/annotation/bin/windows|data/experiment/bin/windows|workflows/bin/windows|books/bin/windows", bad)]   # these repos have no Windows binaries: expected 404s
    if (length(bad)) say("  messages: ", paste(unique(bad), collapse = " | "))
    if (test_load(p)$ok) { say("  ", p, " loads"); break }
    ## a source build that stops with "there is no package called 'X'" names a missing indirect dependency (e.g. stringi,
    ## needed by stringr in the personal library): install X into the PROJECT library and retry; personal library untouched
    bo <- file.path(kodir, paste0(p, ".out")); miss <- character()
    if (file.exists(bo)) { bl <- paste(readLines(bo, warn = FALSE), collapse = " ")
      miss <- unique(gsub("[\u2018\u2019'\"]", "", regmatches(bl, gregexpr("there is no package called [\u2018'\"][^\u2019'\"]+[\u2019'\"]", bl))[[1]]))
      miss <- sub("^there is no package called ", "", miss) }
    miss <- setdiff(miss, c(p, depfixed))
    if (length(miss)) {
      for (q in miss) { say("  build log: missing dependency '", q, "' -> installing it into the project library"); ilog("######## dependency ", q, " for ", p)
        mq <- cap(BiocManager::install(q, lib = projlib, update = FALSE, ask = FALSE, keep_outputs = kodir)); depfixed <- c(depfixed, q)
        say("  ", q, if (test_load(q)$ok) " loads" else paste(" still not loading:", paste(tail(test_load(q)$out, 2), collapse = " "))) }
      next }
    if (!any(grepl(transient, m))) break
    say("  transient network error -> retry in ", 20 * attempt, " s"); Sys.sleep(20 * attempt)
  }
  for (f in list.files(kodir, pattern = "\\.out$", full.names = TRUE)) { say("  source-build log ", basename(f), ": ", paste(tail(grep("ERROR|error|Error", readLines(f, warn = FALSE), value = TRUE), 3), collapse = " | ")) }
  todo <- targets[!vapply(targets, function(p) test_load(p)$ok, logical(1))]
  if (length(todo)) { say("still not loading after route A: ", paste(todo, collapse = ", "), " -> falling back to route B for these"); route <- "A+B" }
}
if (length(todo) && route %in% c("B", "A+B")) {
  py <- c(Sys.getenv("PROJECT_PYTHON"), file.path(Sys.getenv("USERPROFILE"), "anaconda3", "python.exe"), Sys.which("python"))
  py <- py[nzchar(py) & file.exists(py)][1]
  if (is.na(py)) { say("FATAL: no Python found for route B"); quit(status = 4) }
  ip <- installed.packages(lib.loc = LIBS)[, c("Package", "Version")]
  inst <- file.path(outdir, paste0("installed_", ts, ".tsv")); write.table(ip, inst, sep = "\t", quote = FALSE, row.names = FALSE)
  say("Python downloader: ", py, " (requested: ", paste(setdiff(todo, "BiocManager"), collapse = ","), ")")
  o <- suppressWarnings(system2(py, c(shQuote(file.path(proj, "scripts", "R", "fetch_r_packages.py")), shQuote(proj), bioc, rminor, shQuote(inst),
                                      paste(todo, collapse = ","), "--mirrors", paste(unique(c(reachable, MIRRORS)), collapse = ",")), stdout = TRUE, stderr = TRUE))
  ilog(paste(o, collapse = "\n")); say("downloader exit status: ", if (is.null(attr(o, "status"))) 0 else attr(o, "status")); for (l in tail(o, 6)) say("  ", l)
  plan_f <- file.path(proj, "envs", "R-pkgcache", bioc, "install_plan.tsv")
  if (file.exists(plan_f)) {
    plan <- read.delim(plan_f, colClasses = "character"); say("installing ", nrow(plan), " verified files in dependency order")
    for (i in seq_len(nrow(plan))) {
      ilog("######## ", plan$package[i], " ", plan$version[i], " ", plan$file[i])
      m <- cap(install.packages(plan$file[i], repos = NULL, type = if (plan$type[i] == "win.binary") "win.binary" else "source", lib = projlib, keep_outputs = kodir))
      bad <- grep("^ERROR|non-zero exit|failed|cannot|unable", m, value = TRUE, ignore.case = TRUE)
      say(sprintf("  %-22s %-10s %s", plan$package[i], plan$version[i], if (length(bad)) paste("PROBLEM:", paste(unique(bad), collapse = " | ")) else "installed"))
    }
    un <- readLines(file.path(proj, "envs", "R-pkgcache", bioc, "unresolved.txt"), warn = FALSE)
    if (length(un)) say("unresolved (not in any repository index or MD5 mismatch): ", paste(un, collapse = "; "))
  } else say("no install plan produced - see ", instlog)
}

## ---------------------------------------------------------------------------------------------------- 4 VERIFY
say("== 4. VERIFY in a FRESH R session per package (requireNamespace + library)")
res <- do.call(rbind, lapply(c(targets, OPTIONAL), function(p) {
  r <- test_load(p); w <- where(p)
  st <- if (r$ok) "LOADS_OK" else if (is.null(w)) "NOT_INSTALLED" else "INSTALLED_BUT_FAILS_TO_LOAD"
  if (p == OPTIONAL && !r$ok) st <- "NOT_A_BIOCONDUCTOR_PACKAGE (no org.Amellifera.eg.db in Bioconductor; not substituted)"
  e <- if (r$ok) "" else gsub("[\t\r\n]+", " ", paste(tail(r$out, 6), collapse = " | "))
  say(sprintf("  %-22s %-28s %s", p, st, if (r$ok) grep("^LOADED", r$out, value = TRUE) else e))
  data.frame(package = p, group = if (p %in% CORE) "core" else if (p %in% ENRICH) "enrichment" else "optional", status = st,
             version = if (is.null(w)) "" else w$version, library = if (is.null(w)) "" else w$lib, error = e, stringsAsFactors = FALSE) }))
## full inventory of what the core packages actually load (fresh session), for reproducibility
inv <- tempfile(fileext = ".R")
writeLines(c(sprintf(".libPaths(%s)", paste(deparse(LIBS), collapse = "")), "for (p in c('tximport','edgeR','limma','DESeq2')) try(suppressPackageStartupMessages(library(p, character.only=TRUE)))",
             "s <- sessionInfo(); pk <- c(s$otherPkgs, s$loadedOnly)",
             "for (n in sort(names(pk))) { cat(n, pk[[n]]$Version, normalizePath(dirname(find.package(n)), winslash='/'), sep='\\t'); cat('\\n') }"), inv)
io <- suppressWarnings(system2(rscript, c("--vanilla", shQuote(inv, type = if (.Platform$OS.type == "windows") "cmd" else "sh")), stdout = TRUE, stderr = TRUE))
writeLines(c("package\tversion\tlibrary", grep("\t", io, value = TRUE)), file.path(outdir, paste0("R_core_loaded_inventory_", ts, ".tsv")))
write.table(res, file.path(outdir, paste0("R_packages_", ts, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, paste0("sessionInfo_", ts, ".txt")))
core_ok <- all(res$status[res$group == "core"] == "LOADS_OK"); enr_ok <- all(res$status[res$group == "enrichment"] == "LOADS_OK")
cat(R.version.string, "| route", route, "| project library", projlib, "\n", file = file.path(proj, "logs", "tool_versions.txt"), append = TRUE)
for (i in which(res$status == "LOADS_OK")) cat(res$package[i], res$version[i], res$library[i], "\n", file = file.path(proj, "logs", "tool_versions.txt"), append = TRUE)
if (core_ok) writeLines(c(ts, R.version.string, rscript, paste("route", route), paste(LIBS, collapse = " ; "),
                          paste(res$package[res$group == "core"], res$version[res$group == "core"], res$library[res$group == "core"])), file.path(outdir, "CORE_READY.txt"))
say("== SUMMARY: core ", if (core_ok) "ALL LOAD" else "NOT READY", " | enrichment ", if (enr_ok) "ALL LOAD" else paste("NOT ALL LOAD:", paste(res$package[res$group == "enrichment" & res$status != "LOADS_OK"], collapse = ",")),
    " | route ", route, " | report ", rep)
quit(status = if (core_ok) 0 else 1)
