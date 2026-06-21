#!/usr/bin/env Rscript
# ==============================================================================
# Supplementary Table 18: LD stratified by local recombination rate.
#
# Reads the ld_vs_recombrate_1to1.csv produced by step 12
# (__13_LD_vs_RecombRate_1to1.R) and formats it as the published xlsx layout.
#
# Input:
#   <PROJECT_ROOT>/geneconv_complete/LD_analysis_resampling_1to1/ld_vs_recombrate_1to1.csv
#
# Output:
#   <PROJECT_ROOT>/geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_18.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx", repos = "https://cloud.r-project.org")
  library(openxlsx)
})

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))

IN_CSV   <- file.path(project_root, "geneconv_complete/LD_analysis_resampling_1to1/ld_vs_recombrate_1to1.csv")
OUT_XLSX <- file.path(project_root, "geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_18.xlsx")

d <- read.csv(IN_CSV, stringsAsFactors = FALSE)

mapp_lab <- c(allmapp = "All mappable", nosegdupmapp = "Outside segdups", segdupmapp = "Within segdups")
win_kb   <- c(win10kb = 10, win25kb = 25, win50kb = 50, win100kb = 100)
maf_lab  <- c("0_1e-05"="<10⁻⁵", "1e-05_0.0001"="10⁻⁵–10⁻⁴", "0.0001_0.001"="10⁻⁴–10⁻³",
              "0.001_0.005"="0.001–0.005", "0.005_0.01"="0.005–0.01",
              "0.01_0.05"="0.01–0.05", "0.05_0.1"="0.05–0.1",
              "0.1_0.5"="0.1–0.5", "0.5_2"=">0.5")

out <- data.frame(
  Mappability                  = factor(mapp_lab[d$mappability], levels = unname(mapp_lab)),
  Population                   = d$population,
  `Window (kb)`                = win_kb[d$window],
  k                            = d$k,
  `MAF bin`                    = factor(maf_lab[d$maf_bin], levels = unname(maf_lab)),
  `Recomb bin (cM/Mb)`         = d$recomb_bin,
  `N GC`                       = d$gc_n,
  `Mean r² GC`                 = round(d$gc_mean, 5),
  `Mean recomb GC (cM/Mb)`     = round(d$gc_recomb, 3),
  `N non-GC`                   = d$nongc_n,
  `Mean r² non-GC`             = round(d$nongc_mean, 4),
  `Mean recomb non-GC (cM/Mb)` = round(d$nongc_recomb, 3),
  `Δr² × 1000`                 = round(d$diff * 1000, 2),
  `% reduction`                = round(100 * (d$nongc_mean - d$gc_mean) / d$nongc_mean, 1),
  `P-value`                    = signif(d$p_value, 4),
  Significance                 = d$sig,
  check.names = FALSE, stringsAsFactors = FALSE
)
out <- out[order(out$Mappability, out$Population, out$`Window (kb)`, out$k, out$`MAF bin`,
                 out$`Recomb bin (cM/Mb)`), ]

dir.create(dirname(OUT_XLSX), showWarnings = FALSE, recursive = TRUE)
wb <- createWorkbook()
addWorksheet(wb, "LD")
writeData(wb, 1, "Supplementary Table 18: LD stratified by local recombination rate.",
          startRow = 1, startCol = 1)
writeData(wb, 1, out, startRow = 2, startCol = 1)
addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = 2, cols = seq_len(ncol(out)),
         gridExpand = TRUE)
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("Saved: %s (%d rows)\n", OUT_XLSX, nrow(out)))
