#!/usr/bin/env Rscript
# ==============================================================================
# Supplementary Table 15: GC-biased gene conversion — composition of W→S and
# S→W conversion-consistent variants by AF threshold.
#
# Input:
#   <PROJECT_ROOT>/geneconv_complete/distance_analysis/summary/gbgc_proportion_by_af.tsv
#
# Output:
#   <PROJECT_ROOT>/geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_15.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx", repos = "https://cloud.r-project.org")
  library(openxlsx)
})

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))

IN_TSV   <- file.path(project_root, "geneconv_complete/distance_analysis/summary/gbgc_proportion_by_af.tsv")
OUT_XLSX <- file.path(project_root, "geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_15.xlsx")

d <- read.delim(IN_TSV, stringsAsFactors = FALSE)

out <- data.frame(
  k             = d$k,
  `AF threshold`= d$af_threshold,
  `N W→S`  = d$n_WS,
  `N S→W`  = d$n_SW,
  `N total`     = d$n_total,
  `% W→S`  = round(d$pct_WS, 2),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dir.create(dirname(OUT_XLSX), showWarnings = FALSE, recursive = TRUE)
wb <- createWorkbook()
addWorksheet(wb, "Supplementary Table 15")
writeData(wb, 1, "Supplementary Table 15: GC-biased gene conversion — composition of W→S and S→W conversion-consistent variants.",
          startRow = 1, startCol = 1)
writeData(wb, 1, out, startRow = 2, startCol = 1)
addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = 2, cols = seq_len(ncol(out)),
         gridExpand = TRUE)
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("Saved: %s (%d rows)\n", OUT_XLSX, nrow(out)))
