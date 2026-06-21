#!/usr/bin/env Rscript
# ==============================================================================
# Supplementary Table 17: GWAS variant depletion at gene conversion-compatible
# positions.
#
# Reads the StatsSummary.R.out.txt files produced by step 06
# (__10_SummarizeContingencyStatistics.pl) on GWAS Catalog input, for each
# mappability stratum (all mappable / outside segdups / within segdups).
# Computes per-row GC and non-GC variant frequencies and writes a 3-sheet xlsx.
#
# Input:
#   <PROJECT_ROOT>/geneconv_complete/stats_{allmapp,nosegdupmapp,segdupmapp}/
#     _All_1000000_1000000_gwas~clean~All_StatsSummary.R.out.txt
#
# Output:
#   <PROJECT_ROOT>/geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_17.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx", repos = "https://cloud.r-project.org")
  library(openxlsx)
})

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))

OUT_XLSX <- file.path(project_root, "geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_17.xlsx")

af_labels <- c("0_1e-05"="<10⁻⁵", "1e-05_00001"="10⁻⁵–10⁻⁴", "00001_0001"="10⁻⁴–10⁻³",
               "0001_0005"="10⁻³–5×10⁻³", "0005_001"="5×10⁻³–10⁻²",
               "001_005"="10⁻²–5×10⁻²", "005_01"="5×10⁻²–0.1",
               "01_05"="0.1–0.5", "05_2"=">0.5")
k_values <- c(17,19,21,31,41,51,61,71,81,91)
N_TESTS  <- length(k_values) * length(af_labels)

build_panel <- function(rout_path) {
  d <- read.delim(rout_path, stringsAsFactors = FALSE)
  d <- d[d$RepLength %in% k_values & d$FrequencyInterval %in% names(af_labels), ]
  d$AF.bin       <- factor(af_labels[d$FrequencyInterval], levels = unname(af_labels))
  d$GC.frequency      <- d$NrOfVarPosConvPos   / (d$NrOfVarPosConvPos   + d$NrOfNoVarPosConvPos)
  d$Non.GC.frequency  <- d$NrOfVarPosNoConvPos / (d$NrOfVarPosNoConvPos + d$NrOfNoVarPosNoConvPos)
  d$log10p_bonf_enrichment  <- pmin(d$log_p_value_exponent           + log10(N_TESTS), 0)
  d$log10p_bonf_concordance <- pmin(d$log_varmatch_p_value_exponent  + log10(N_TESTS), 0)
  o <- data.frame(
    k = d$RepLength,
    `AF bin` = d$AF.bin,
    `N GWAS at GC pos`        = d$NrOfVarPosConvPos,
    `N GWAS not at GC pos`    = d$NrOfVarPosNoConvPos,
    `N no-var GC pos`         = d$NrOfNoVarPosConvPos,
    `N no-var non-GC pos`     = d$NrOfNoVarPosNoConvPos,
    `GC frequency`            = signif(d$GC.frequency, 6),
    `Non-GC frequency`        = signif(d$Non.GC.frequency, 6),
    `χ² p-value (enrichment)` = signif(d$chisq_test_p, 6),
    `log₁₀(p) enrichment`     = round(d$log10p_bonf_enrichment, 2),
    `Concordant variants`     = d$ConcorVar,
    `Discordant variants`     = d$DiscorVar,
    `χ² p-value (concordance)`= signif(d$chisq_test_varmatch_p, 6),
    `log₁₀(p) concordance`    = round(d$log10p_bonf_concordance, 2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  o <- o[order(d$RepLength, d$AF.bin), ]
  o
}

strata <- list(
  "All mappable"    = c("stats_allmapp",     "Supplementary Table 17a: GWAS variant depletion at gene conversion-compatible positions (all mappable regions)."),
  "Outside segdups" = c("stats_nosegdupmapp","Supplementary Table 17b: GWAS variant depletion at gene conversion-compatible positions, outside segmental duplications."),
  "Within segdups"  = c("stats_segdupmapp",  "Supplementary Table 17c: GWAS variant depletion at gene conversion-compatible positions, within segmental duplications.")
)

wb <- createWorkbook()
for (sheet in names(strata)) {
  subdir <- strata[[sheet]][1]
  title  <- strata[[sheet]][2]
  rout <- file.path(project_root, "geneconv_complete", subdir,
                    "_All_1000000_1000000_gwas~clean~All_StatsSummary.R.out.txt")
  panel <- build_panel(rout)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, title, startRow = 1, startCol = 1)
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  # Sub-panel header (matches the Enrichment/Concordance 2x2 layout of ST3/ST7):
  # cols 3-10 = enrichment 2x2 + stats, cols 11-14 = concordance 2x2 + stats.
  writeData(wb, sheet, "Enrichment (2×2 table)",  startRow = 2, startCol = 3)
  writeData(wb, sheet, "Concordance (2×2 table)", startRow = 2, startCol = 11)
  mergeCells(wb, sheet, cols = 3:10,  rows = 2)
  mergeCells(wb, sheet, cols = 11:14, rows = 2)
  writeData(wb, sheet, panel, startRow = 3, startCol = 1)
  drows <- 4:(3 + nrow(panel))
  # Colour scheme matching ST3/ST7: k/AF header blue, enrichment green,
  # concordance orange; data rows in the lighter tint of each group.
  addStyle(wb, sheet, createStyle(fgFill = "#E2EFDA", textDecoration = "bold", halign = "center"),
           rows = 2, cols = 3:10,  gridExpand = TRUE)
  addStyle(wb, sheet, createStyle(fgFill = "#FCE4D6", textDecoration = "bold", halign = "center"),
           rows = 2, cols = 11:14, gridExpand = TRUE)
  addStyle(wb, sheet, createStyle(fgFill = "#D9E1F2", textDecoration = "bold", border = "bottom"),
           rows = 3, cols = 1:2,   gridExpand = TRUE)
  addStyle(wb, sheet, createStyle(fgFill = "#E2EFDA", textDecoration = "bold", border = "bottom"),
           rows = 3, cols = 3:10,  gridExpand = TRUE)
  addStyle(wb, sheet, createStyle(fgFill = "#FCE4D6", textDecoration = "bold", border = "bottom"),
           rows = 3, cols = 11:14, gridExpand = TRUE)
  addStyle(wb, sheet, createStyle(fgFill = "#F2F7ED"), rows = drows, cols = 3:10,  gridExpand = TRUE)
  addStyle(wb, sheet, createStyle(fgFill = "#FDF0E8"), rows = drows, cols = 11:14, gridExpand = TRUE)
  setColWidths(wb, sheet, cols = 1:ncol(panel), widths = "auto")
  cat(sprintf("Sheet %s: %d rows\n", sheet, nrow(panel)))
}
dir.create(dirname(OUT_XLSX), showWarnings = FALSE, recursive = TRUE)
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("Saved: %s\n", OUT_XLSX))
