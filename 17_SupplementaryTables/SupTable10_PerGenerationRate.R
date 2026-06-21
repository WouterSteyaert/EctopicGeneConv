#!/usr/bin/env Rscript
# ==============================================================================
# Supplementary Table 10: Per-generation ectopic gene conversion rate.
#
# Reads the per_generation_gc_rate.tsv produced by step 14
# (PerGenerationConversionRate.R) and assembles the published xlsx layout:
# three stacked stratum blocks (a/b/c), each with grouped column headers
#   De novo mutation counts | Conversion rate | Fold excess | r_gc / mu_eff CI
# and the ST3/ST7-consistent colour scheme.
#
# Column layout per block (A..O):
#   A k | B Positions | C Concordant | D Discordant | E Excess
#   F r_gc | G lower | H upper | I mu_eff | J r_gc/mu_eff
#   K fold Estimate | L lower | M upper | N r_gc/mu_eff lower | O upper
# Columns N-O are the r_gc/mu_eff CI (= fold CI - 1), as referenced in the
# Methods ("Supplementary Table 10 columns N-O").
#
# Input:
#   <PROJECT_ROOT>/geneconv_complete/dnm_analysis/export/per_generation_gc_rate.tsv
#
# Output:
#   <PROJECT_ROOT>/geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_10.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx", repos = "https://cloud.r-project.org")
  library(openxlsx)
})

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))

N_TRIOS  <- 11963
IN_TSV   <- file.path(project_root, "geneconv_complete/dnm_analysis/export/per_generation_gc_rate.tsv")
OUT_XLSX <- file.path(project_root, "geneconv_complete/__PAPER_3/TablesFinal/SupplementaryTable_10.xlsx")

d <- read.delim(IN_TSV, stringsAsFactors = FALSE)

# Map the tsv stratum labels to the published block titles, in published order.
blocks <- list(
  list(stratum = "All mappable",           title = "(a) All mappable regions"),
  list(stratum = "Segmental duplications", title = "(b) Within segmental duplications"),
  list(stratum = "Non-segdup mappable",    title = "(c) Outside segmental duplications")
)
k_order <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)

# ── Colour scheme (matches ST3 / ST7 / ST17) ────────────────────────────────
C_BLUE    <- "#D9E1F2"; C_GREEN   <- "#E2EFDA"; C_MBLUE  <- "#D6ECF5"
C_ORANGE  <- "#FCE4D6"; C_PORANGE <- "#FDE8D4"
D_GREEN   <- "#F2F7ED"; D_BLUE    <- "#EBF5FB"; D_ORANGE <- "#FDF0E8"; D_PORANGE <- "#FEF4EE"
BORDER_COL <- "#B0B0B0"
AR <- function(...) createStyle(fontName = "Arial", fontSize = 9, ...)

# Group-header styles (no bottom border)
gh_green   <- AR(textDecoration = "bold", fgFill = C_GREEN,   halign = "center")
gh_mblue   <- AR(textDecoration = "bold", fgFill = C_MBLUE,   halign = "center")
gh_orange  <- AR(textDecoration = "bold", fgFill = C_ORANGE,  halign = "center")
gh_porange <- AR(textDecoration = "bold", fgFill = C_PORANGE, halign = "center")
# Column-header styles (bottom border, except N-O to match the published table)
ch_blue    <- AR(textDecoration = "bold", fgFill = C_BLUE,    halign = "center", border = "bottom", borderColour = BORDER_COL)
ch_green   <- AR(textDecoration = "bold", fgFill = C_GREEN,   halign = "center", border = "bottom", borderColour = BORDER_COL)
ch_mblue   <- AR(textDecoration = "bold", fgFill = C_MBLUE,   halign = "center", border = "bottom", borderColour = BORDER_COL)
ch_orange  <- AR(textDecoration = "bold", fgFill = C_ORANGE,  halign = "center", border = "bottom", borderColour = BORDER_COL)
ch_porange <- AR(textDecoration = "bold", fgFill = C_PORANGE, halign = "center")
# Data styles
da_k    <- AR(halign = "center")
da_pos  <- AR(numFmt = "#,##0",     halign = "center")
da_cnt  <- AR(numFmt = "#,##0",     fgFill = D_GREEN,   halign = "center")
da_exc  <- AR(numFmt = "#,##0.0",   fgFill = D_GREEN,   halign = "center")
da_sci  <- AR(numFmt = "0.00E+00",  fgFill = D_BLUE,    halign = "center")
da_rat  <- AR(numFmt = "0.00",      fgFill = D_BLUE,    halign = "center")
da_fold <- AR(numFmt = "0.00",      fgFill = D_ORANGE,  halign = "center")
da_ci   <- AR(                      fgFill = D_PORANGE, halign = "center")
# Title / block / footnote styles
st_title <- AR(textDecoration = "bold")
st_block <- AR(textDecoration = "bold", halign = "left")
st_foot  <- AR(valign = "top", wrapText = TRUE)

col_headers <- c("k", "Positions", "Concordant", "Discordant", "Excess",
                 "r_gc", "lower", "upper", "μ_eff", "r_gc / μ_eff",
                 "Estimate", "lower", "upper", "lower", "upper")

wb <- createWorkbook()
sheet <- "Supplementary Table 10"
addWorksheet(wb, sheet)

writeData(wb, sheet, "Supplementary Table 10: Per-generation ectopic gene conversion rate.",
          startRow = 1, startCol = 1)
mergeCells(wb, sheet, cols = 1:15, rows = 1)
addStyle(wb, sheet, st_title, rows = 1, cols = 1)

r <- 3
for (b in blocks) {
  sub <- d[d$stratum == b$stratum, ]
  sub <- sub[match(k_order, sub$k), ]

  # ── Block title ───────────────────────────────────────────────────────────
  writeData(wb, sheet, b$title, startRow = r, startCol = 1)
  mergeCells(wb, sheet, cols = 1:15, rows = r)
  addStyle(wb, sheet, st_block, rows = r, cols = 1)

  # ── Group-header row ──────────────────────────────────────────────────────
  gr <- r + 1
  writeData(wb, sheet, "De novo mutation counts",               startRow = gr, startCol = 3)
  writeData(wb, sheet, "Conversion rate",                       startRow = gr, startCol = 6)
  writeData(wb, sheet, "Fold excess (concordant / discordant)", startRow = gr, startCol = 11)
  writeData(wb, sheet, "r_gc / μ_eff CI",                  startRow = gr, startCol = 14)
  mergeCells(wb, sheet, cols = 3:5,   rows = gr)
  mergeCells(wb, sheet, cols = 6:10,  rows = gr)
  mergeCells(wb, sheet, cols = 11:13, rows = gr)
  mergeCells(wb, sheet, cols = 14:15, rows = gr)
  addStyle(wb, sheet, gh_green,   rows = gr, cols = 3:5,   gridExpand = TRUE)
  addStyle(wb, sheet, gh_mblue,   rows = gr, cols = 6:10,  gridExpand = TRUE)
  addStyle(wb, sheet, gh_orange,  rows = gr, cols = 11:13, gridExpand = TRUE)
  addStyle(wb, sheet, gh_porange, rows = gr, cols = 14:15, gridExpand = TRUE)

  # ── Column-header row ─────────────────────────────────────────────────────
  hr <- r + 2
  writeData(wb, sheet, as.list(col_headers), startRow = hr, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet, ch_blue,    rows = hr, cols = 1:2,   gridExpand = TRUE)
  addStyle(wb, sheet, ch_green,   rows = hr, cols = 3:5,   gridExpand = TRUE)
  addStyle(wb, sheet, ch_mblue,   rows = hr, cols = 6:10,  gridExpand = TRUE)
  addStyle(wb, sheet, ch_orange,  rows = hr, cols = 11:13, gridExpand = TRUE)
  addStyle(wb, sheet, ch_porange, rows = hr, cols = 14:15, gridExpand = TRUE)

  # ── Data ──────────────────────────────────────────────────────────────────
  # r_gc/mu_eff lower (col N) carries a footnote dagger where it was clipped to 0
  # (normal approximation negative; see footnote). Written as character so the
  # dagger can be attached; other CI values keep their 2-decimal display.
  n_lo <- as.character(round(sub$rgc_over_mu_eff_lo, 2))
  clip <- sub$rgc_over_mu_eff_lo == 0 & (sub$fold_lo - 1) < 0
  n_lo[clip] <- paste0(n_lo[clip], "†")

  blk <- data.frame(
    k        = sub$k,
    Positions= sub$ConcorTotal,
    Concor   = sub$ConcorVar,
    Discor   = sub$DiscorVar,
    Excess   = round(sub$ExcessConcordant, 1),
    r_gc     = as.numeric(sub$r_gc_direct),
    rgc_lo   = as.numeric(sub$r_gc_lo),
    rgc_hi   = as.numeric(sub$r_gc_hi),
    mu_eff   = as.numeric(sub$mu_eff_per_allele),
    ratio    = round(sub$rgc_over_mu_eff, 2),
    fold     = round(sub$fold_excess, 2),
    fold_lo  = round(sub$fold_lo, 2),
    fold_hi  = round(sub$fold_hi, 2),
    ci_lo    = n_lo,
    ci_hi    = round(sub$rgc_over_mu_eff_hi, 2),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  dr <- r + 3
  writeData(wb, sheet, blk, startRow = dr, startCol = 1, colNames = FALSE)
  drows <- dr:(dr + nrow(blk) - 1)
  addStyle(wb, sheet, da_k,    rows = drows, cols = 1,     gridExpand = TRUE)
  addStyle(wb, sheet, da_pos,  rows = drows, cols = 2,     gridExpand = TRUE)
  addStyle(wb, sheet, da_cnt,  rows = drows, cols = 3:4,   gridExpand = TRUE)
  addStyle(wb, sheet, da_exc,  rows = drows, cols = 5,     gridExpand = TRUE)
  addStyle(wb, sheet, da_sci,  rows = drows, cols = 6:9,   gridExpand = TRUE)
  addStyle(wb, sheet, da_rat,  rows = drows, cols = 10,    gridExpand = TRUE)
  addStyle(wb, sheet, da_fold, rows = drows, cols = 11:13, gridExpand = TRUE)
  addStyle(wb, sheet, da_ci,   rows = drows, cols = 14:15, gridExpand = TRUE)

  r <- r + 14   # 13 rows used (title + 2 headers + 10 data) + 1 blank spacer
}

# ── Footnotes ────────────────────────────────────────────────────────────────
foot <- c(
  "N trios = 11,963. Positions = number of gene conversion-compatible positions for each k.",
  "r_gc = excess concordant DNMs / (N × positions). μ_eff = effective per-allele point mutation rate (from discordant DNM density). 95% CI: Poisson approximation.",
  paste0("Fold excess = concordant / discordant DNMs. r_gc / μ_eff = fold excess − 1. ",
         "† Lower bound clipped to 0: normal approximation yields negative value ",
         "(k = 17, within SDs: −3.001 × 10⁻¹⁰; confirmed under exact Poisson CI).")
)
fr <- r   # r already points one row past the block's trailing spacer
for (i in seq_along(foot)) {
  writeData(wb, sheet, foot[i], startRow = fr + i - 1, startCol = 1)
  mergeCells(wb, sheet, cols = 1:15, rows = fr + i - 1)
  addStyle(wb, sheet, st_foot, rows = fr + i - 1, cols = 1)
}

setColWidths(wb, sheet, cols = 1:15, widths = "auto")

dir.create(dirname(OUT_XLSX), showWarnings = FALSE, recursive = TRUE)
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("Saved: %s (3 stratum blocks, %d k-values each)\n", OUT_XLSX, length(k_order)))
