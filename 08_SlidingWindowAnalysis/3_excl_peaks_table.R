#!/usr/bin/env Rscript
# =============================================================================
# Step 3: Build "excluding peaks" supplementary table (Supplementary Table 12)
#         from corrected sliding-window data.
#
# For each k x MAF bin, over non-overlapping 1-Mb windows, compute:
#   (a) mean log2(OR) and % positive windows, all windows
#   (b) the same, after excluding windows that overlap corrected peaks
#   (c) a genome-wide 2x2 chi-square enrichment test on the summed contingency
#       cells (gc_var, gc_novar, nongc_var, nongc_novar), both for all windows
#       and after peak exclusion, reported as Bonferroni-corrected log10(p)
#       with m = 90 (10 template lengths x 9 AF bins). A cell is called
#       "significant" only if it is both enriched (gc rate > non-gc rate) and
#       passes the Bonferroni threshold (corrected p < 0.05). This substantiates
#       the Results statement that enrichment survives peak exclusion for all
#       but one (k, AF) combination.
#
# Input:  $PROJECT_ROOT/geneconv_complete/sliding_window_corrected/
#             enrichment_corrected_combined.tsv
#             peaks_corrected.tsv
# Output: same dir / SupTable_EnrichmentExclPeaks_corrected.tsv
# =============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
SWDIR  <- file.path(project_root, "geneconv_complete", "sliding_window_corrected")
COMB   <- file.path(SWDIR, "enrichment_corrected_combined.tsv")
PEAKS  <- file.path(SWDIR, "peaks_corrected.tsv")
OUT    <- file.path(SWDIR, "SupTable_EnrichmentExclPeaks_corrected.tsv")

BONF_M       <- 90                       # 10 k x 9 AF bins
THRESH_LOG10 <- log10(0.05)              # corrected log10(p) must be below this

comb <- read.delim(COMB, stringsAsFactors = FALSE)
peaks <- read.delim(PEAKS, stringsAsFactors = FALSE)

# numeric coercion (guard against stray non-numeric tokens)
for (col in c("k", "window_start", "window_end", "gc_var", "gc_novar",
              "nongc_var", "nongc_novar", "log2_or")) {
  comb[[col]] <- suppressWarnings(as.numeric(comb[[col]]))
}

# non-overlapping 1-Mb windows only (start on the 1-Mb grid), drop incomplete rows
comb <- comb[((comb$window_start - 1) %% 1000000) == 0, ]
comb <- comb[stats::complete.cases(
  comb[, c("gc_var", "gc_novar", "nongc_var", "nongc_novar", "log2_or")]), ]

# peak-window lookup: map each peak interval to the 1-Mb grid windows it covers,
# keyed by chr:windowStart:k:maf  (identical mapping to the original awk builder)
peak_keys <- new.env(hash = TRUE, parent = emptyenv())
if (nrow(peaks) > 0) {
  for (i in seq_len(nrow(peaks))) {
    ps <- peaks$Peak_start[i]; pe <- peaks$Peak_end[i]
    s0 <- as.integer((ps - 1) %/% 1000000) * 1000000 + 1
    for (s in seq(s0, pe, by = 1000000)) {
      assign(paste(peaks$Chromosome[i], s, peaks$Template_length_k[i],
                   peaks$MAF_bin[i], sep = ":"), TRUE, envir = peak_keys)
    }
  }
}
is_peak_window <- function(chr, ws, k, maf) {
  exists(paste(chr, ws, k, maf, sep = ":"), envir = peak_keys, inherits = FALSE)
}

# 2x2 Pearson chi-square, df = 1, returns RAW log10(p)
log10p_chisq <- function(a, b, c, d) {
  n <- a + b + c + d
  if (!is.finite(n) || n == 0) return(NA_real_)
  r1 <- a + b; r2 <- c + d; c1 <- a + c; c2 <- b + d
  if (r1 == 0 || r2 == 0 || c1 == 0 || c2 == 0) return(0)  # degenerate -> p = 1
  ea <- r1 * c1 / n; eb <- r1 * c2 / n; ec <- r2 * c1 / n; ed <- r2 * c2 / n
  chi <- (a - ea)^2 / ea + (b - eb)^2 / eb + (c - ec)^2 / ec + (d - ed)^2 / ed
  pchisq(chi, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
}

ks   <- sort(unique(comb$k))
mafs <- sort(unique(comb$maf_bin))
rows <- list()
for (k in ks) for (maf in mafs) {
  w <- comb[comb$k == k & comb$maf_bin == maf, ]
  if (nrow(w) == 0) next
  pk <- vapply(seq_len(nrow(w)),
               function(j) is_peak_window(w$chr[j], w$window_start[j], k, maf),
               logical(1))
  we <- w[!pk, ]

  cells <- function(s) c(gv = sum(s$gc_var), gn = sum(s$gc_novar),
                         nv = sum(s$nongc_var), nn = sum(s$nongc_novar))
  A <- cells(w); E <- cells(we)
  enr <- function(x) (x["gv"] / (x["gv"] + x["gn"])) > (x["nv"] / (x["nv"] + x["nn"]))

  lp_all_raw  <- log10p_chisq(A["gv"], A["gn"], A["nv"], A["nn"])
  lp_excl_raw <- if (nrow(we) > 0) log10p_chisq(E["gv"], E["gn"], E["nv"], E["nn"]) else NA_real_
  lp_all_bonf  <- lp_all_raw  + log10(BONF_M)
  lp_excl_bonf <- lp_excl_raw + log10(BONF_M)

  sig_all  <- is.finite(lp_all_bonf)  && lp_all_bonf  < THRESH_LOG10 && isTRUE(unname(enr(A)))
  sig_excl <- is.finite(lp_excl_bonf) && lp_excl_bonf < THRESH_LOG10 && isTRUE(unname(enr(E)))

  rows[[length(rows) + 1]] <- data.frame(
    k = k, maf_bin = maf,
    n_windows_all = nrow(w),
    mean_log2OR_all = round(mean(w$log2_or), 4),
    pct_positive_all = round(100 * mean(w$log2_or > 0), 1),
    log10p_bonf_all = round(lp_all_bonf, 1),
    significant_all = sig_all,
    n_windows_excl = nrow(we),
    mean_log2OR_excl = if (nrow(we) > 0) round(mean(we$log2_or), 4) else NA_real_,
    pct_positive_excl = if (nrow(we) > 0) round(100 * mean(we$log2_or > 0), 1) else NA_real_,
    log10p_bonf_excl = round(lp_excl_bonf, 1),
    significant_excl = sig_excl,
    stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows)
tab <- tab[order(tab$k, tab$maf_bin), ]
write.table(tab, OUT, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Wrote %d rows to %s\n", nrow(tab), OUT))
cat(sprintf("Bonferroni m = %d; corrected log10(p) threshold = %.3f\n",
            BONF_M, THRESH_LOG10))
cat(sprintf("Significant+enriched, ALL windows:  %d / %d\n",
            sum(tab$significant_all), nrow(tab)))
cat(sprintf("Significant+enriched, EXCL peaks:   %d / %d\n",
            sum(tab$significant_excl), nrow(tab)))
cat("\nCells NOT significant+enriched after peak exclusion:\n")
print(tab[!tab$significant_excl,
          c("k", "maf_bin", "n_windows_excl", "mean_log2OR_excl",
            "log10p_bonf_excl", "significant_excl")], row.names = FALSE)
