#!/usr/bin/env Rscript
# ==============================================================================
# Supplementary Table 16: Linkage disequilibrium at gene conversion-compatible
# positions.
#
# Reads the LD csv produced by step 12 (__12_AnalyzeLd_1to1.R), formats it for
# the published TablesFinal layout (full mappability × population × window ×
# k × AF combinations with mean r², n, Δr², p-value, significance).
#
# Input:
#   <PROJECT_ROOT>/geneconv_complete/LD_analysis_resampling_1to1/ld_analysis_results_1to1.csv
#
# Output:
#   <PROJECT_ROOT>/geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_16.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx", repos = "https://cloud.r-project.org")
  library(openxlsx)
})

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))

IN_CSV   <- file.path(project_root, "geneconv_complete/LD_analysis_resampling_1to1/ld_analysis_results_1to1.csv")
OUT_XLSX <- file.path(project_root, "geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_16.xlsx")

ld <- read.csv(IN_CSV, stringsAsFactors = FALSE)

mapp_labels   <- c(allmapp = "All mappable", nosegdupmapp = "Outside segdups", segdupmapp = "Within segdups")
maf_labels    <- c("0_1e-05"="<10⁻⁵", "1e-05_0.0001"="10⁻⁵–10⁻⁴", "0.0001_0.001"="10⁻⁴–10⁻³",
                   "0.001_0.005"="0.001–0.005", "0.005_0.01"="0.005–0.01",
                   "0.01_0.05"="0.01–0.05", "0.05_0.1"="0.05–0.1",
                   "0.1_0.5"="0.1–0.5", "0.5_2"=">0.5")
window_kb     <- c(win10kb = 10, win25kb = 25, win50kb = 50, win100kb = 100)

ld$Mappability <- factor(mapp_labels[ld$mappability], levels = unname(mapp_labels))
ld$`Window (kb)` <- window_kb[ld$window]
ld$`MAF bin`   <- factor(maf_labels[ld$maf_bin], levels = unname(maf_labels))
ld$`% reduction` <- round(100 * (ld$nongc_mean - ld$gc_mean) / ld$nongc_mean, 1)
ld$`Δr² × 1000` <- round(ld$diff * 1000, 2)

out <- data.frame(
  Mappability        = ld$Mappability,
  Population         = ld$population,
  `Window (kb)`      = ld$`Window (kb)`,
  k                  = ld$k,
  `MAF bin`          = ld$`MAF bin`,
  `N GC`             = ld$gc_n,
  `N non-GC`         = ld$nongc_n,
  `Mean r² GC`       = round(ld$gc_mean, 4),
  `Mean r² non-GC`   = round(ld$nongc_mean, 4),
  `Δr² × 1000`       = ld$`Δr² × 1000`,
  `% reduction`      = ld$`% reduction`,
  `P-value`          = signif(ld$p_value, 4),
  Significance       = ld$sig,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
out <- out[order(out$Mappability, out$Population, out$`Window (kb)`, out$k, out$`MAF bin`), ]

dir.create(dirname(OUT_XLSX), showWarnings = FALSE, recursive = TRUE)
wb <- createWorkbook()
addWorksheet(wb, "LD")
writeData(wb, 1, "Supplementary Table 16: Linkage disequilibrium at gene conversion-compatible positions.",
          startRow = 1, startCol = 1)
writeData(wb, 1, out, startRow = 2, startCol = 1)
addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = 2, cols = seq_len(ncol(out)),
         gridExpand = TRUE)
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("Saved: %s (%d rows)\n", OUT_XLSX, nrow(out)))
