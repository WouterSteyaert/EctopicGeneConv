#!/usr/bin/env Rscript
# ==============================================================================
# Supplementary Table: Sliding-window enrichment peaks (Supplementary Table 11)
#
# Reads the corrected peak set produced by step 08 (the canonical paper version
# with control-group correction + MIN_GC_VAR = 10) and emits a clean TSV with
# display-formatted MAF bin labels, sorted by chromosome and position.
#
# Inputs:  <PROJECT_ROOT>/geneconv_complete/sliding_window_corrected/peaks_corrected.tsv
# Outputs: <PROJECT_ROOT>/geneconv_complete/Tables/SupTable_SlidingWindowPeaks.tsv
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))

SW_DIR  <- file.path(project_root, "geneconv_complete/sliding_window_corrected")
OUT_DIR <- file.path(project_root, "geneconv_complete/Tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

peaks <- read.delim(file.path(SW_DIR, "peaks_corrected.tsv"),
                    header = TRUE, stringsAsFactors = FALSE)

# ---- MAF bin display labels ----
maf_labels <- c(
  "0_1e-05"     = "0 - 1e-5",
  "1e-05_00001" = "1e-5 - 1e-4",
  "00001_0001"  = "1e-4 - 1e-3",
  "0001_0005"   = "1e-3 - 5e-3",
  "0005_001"    = "5e-3 - 0.01",
  "001_005"     = "0.01 - 0.05",
  "005_01"      = "0.05 - 0.1",
  "01_05"       = "0.1 - 0.5",
  "05_2"        = "0.5 - 2.0"
)

# ---- Build clean table ----
sup <- data.frame(
  Chromosome          = peaks$Chromosome,
  Peak_start          = peaks$Peak_start,
  Peak_end            = peaks$Peak_end,
  Width_Mb            = peaks$Width_Mb,
  Template_length_k   = peaks$Template_length_k,
  MAF_bin             = maf_labels[peaks$MAF_bin],
  Max_log2_OR         = round(peaks$Max_log2_OR, 3),
  Mean_log2_OR        = round(peaks$Mean_log2_OR, 3),
  Chromosomal_region  = peaks$Chromosomal_region,
  Centromere_dist_Mb  = peaks$Centromere_dist_Mb,
  Telomere_dist_Mb    = peaks$Telomere_dist_Mb,
  Segdup_overlap      = peaks$Segdup_overlap,
  Segdup_fraction     = round(peaks$Segdup_fraction, 3),
  stringsAsFactors    = FALSE
)

# ---- Sort by chromosome (numeric order) then position ----
chr_order <- c(as.character(1:22), "X", "Y")
sup$chr_rank <- match(sup$Chromosome, chr_order)
sup <- sup[order(sup$chr_rank, sup$Peak_start, sup$Template_length_k), ]
sup$chr_rank <- NULL
rownames(sup) <- NULL

# ---- Write ----
outfile <- file.path(OUT_DIR, "SupTable_SlidingWindowPeaks.tsv")
write.table(sup, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Written: %s\n", outfile))
cat(sprintf("Rows: %d\n", nrow(sup)))
cat(sprintf("Columns: %s\n", paste(names(sup), collapse = ", ")))
