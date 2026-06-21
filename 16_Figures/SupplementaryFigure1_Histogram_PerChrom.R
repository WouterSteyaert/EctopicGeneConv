#!/usr/bin/env Rscript
# ==============================================================================
# Supplementary Figure 1: per-chromosome log2(OR) histograms
#
# For each chromosome, histograms of log2(OR) across non-overlapping 1-Mb
# windows, faceted by template length k (rows) and allele frequency bin
# (columns). Four representative k values are shown (k = 21, 41, 61, 91).
#
# Reference lines:
#   grey dashed at log2(OR) = 0  (null)
#   red solid    at the per-panel median
#
# Input:
#   <PROJECT_ROOT>/geneconv_complete/stats_allmapp/{k}/
#     {chr}_1000000_1000000_gnomad~genome_{maf}.1.{chromsize}.txt
#
# Output:
#   <PROJECT_ROOT>/geneconv_complete/figures/SupplementaryFigure1_Histogram_PerChrom.pdf
#     (24 pages, one per chromosome)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- Paths ----
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
STATS_DIR <- file.path(project_root, "geneconv_complete/stats_allmapp")
FIG_DIR   <- file.path(project_root, "geneconv_complete/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Parameters ----
k_values <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)

maf_bins <- c("0_1e-05", "1e-05_00001", "00001_0001", "0001_0005",
              "0005_001", "001_005", "005_01", "01_05", "05_2")

maf_labels <- c(
  "0_1e-05"     = "AF < 1e-5",
  "1e-05_00001" = "1e-5 - 1e-4",
  "00001_0001"  = "1e-4 - 1e-3",
  "0001_0005"   = "1e-3 - 5e-3",
  "0005_001"    = "5e-3 - 0.01",
  "001_005"     = "0.01 - 0.05",
  "005_01"      = "0.05 - 0.1",
  "01_05"       = "0.1 - 0.5",
  "05_2"        = "> 0.5"
)

chrom_sizes <- c(
  "1"=248956422, "2"=242193529, "3"=198295559, "4"=190214555,
  "5"=181538259, "6"=170805979, "7"=159345973, "8"=145138636,
  "9"=138394717, "10"=133797422, "11"=135086622, "12"=133275309,
  "13"=114364328, "14"=107043718, "15"=101991189, "16"=90338345,
  "17"=83257441, "18"=80373285, "19"=58617616, "20"=64444167,
  "21"=46709983, "22"=50818468, "X"=156040895, "Y"=57227415
)
chroms <- names(chrom_sizes)

# ---- Read all data ----
cat("Reading data...\n")

all_data <- list()
for (k in k_values) {
  k_dir <- file.path(STATS_DIR, k)
  for (maf in maf_bins) {
    rows <- list()
    for (chr in chroms) {
      fname <- sprintf("%s_1000000_1000000_gnomad~genome_%s.1.%d.txt",
                       chr, maf, chrom_sizes[chr])
      fpath <- file.path(k_dir, fname)
      if (!file.exists(fpath)) next
      d <- read.delim(fpath, header = TRUE, stringsAsFactors = FALSE)
      d$gc_var      <- d$NrOfVarPosConvPos
      d$nongc_var   <- d$NrOfVarPosNoConvPos
      d$gc_novar    <- d$NrOfNoVarPosConvPos
      d$nongc_novar <- d$NrOfNoVarPosNoConvPos
      d$or <- (d$gc_var * d$nongc_novar) / (d$nongc_var * d$gc_novar)
      d$log2_or <- log2(d$or)
      d <- d[is.finite(d$log2_or), ]
      if (nrow(d) > 0) {
        rows[[length(rows) + 1]] <- data.frame(
          chr     = chr,
          window  = d$WindowStart,
          log2_or = d$log2_or,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(rows) > 0) {
      df <- do.call(rbind, rows)
      df$k   <- k
      df$maf <- maf
      all_data[[length(all_data) + 1]] <- df
    }
  }
  cat(sprintf("  k = %d done\n", k))
}

dat <- do.call(rbind, all_data)
cat(sprintf("Total windows: %d\n", nrow(dat)))

# ---- Factor labels ----
dat$k_label <- factor(paste0("k = ", dat$k),
                      levels = paste0("k = ", k_values))
dat$maf_label <- factor(maf_labels[dat$maf],
                        levels = maf_labels)

# ---- Helper: histogram plot function ----
# drop = FALSE so chromosomes with no variants in a given panel (e.g. chrY
# at AF < 1e-5) still produce an empty facet, preserving a uniform grid.
make_histogram <- function(data, base_size = 9) {
  medians <- aggregate(log2_or ~ k_label + maf_label, data = data, FUN = median)

  ggplot(data, aes(x = log2_or)) +
    geom_histogram(bins = 50, fill = "steelblue", colour = NA, alpha = 0.7) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = "dashed") +
    geom_vline(data = medians, aes(xintercept = log2_or),
               colour = "red", linewidth = 0.5) +
    facet_grid(k_label ~ maf_label, scales = "free_y", drop = FALSE) +
    labs(x = expression(log[2]~"odds ratio"),
         y = "Number of 1-Mb windows") +
    theme_bw(base_size = base_size) +
    theme(
      strip.text       = element_text(size = base_size - 1, face = "bold"),
      axis.text        = element_text(size = base_size - 2),
      axis.title       = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(5, 5, 5, 5)
    )
}

# ===========================================================================
# Per-chromosome supplement: 1 page per chromosome, 4 k x 9 MAF
# Representative template lengths (k = 21, 41, 61, 91)
# ===========================================================================
cat("\n=== Per-chromosome supplement (up to 24 pages) ===\n")

chr_k <- c(21, 41, 61, 91)
chrom_order <- c(as.character(1:22), "X", "Y")

chr_file <- file.path(FIG_DIR, "SupplementaryFigure1_Histogram_PerChrom.pdf")
pdf(chr_file, width = 18, height = 10, family = "sans")

for (ch in chrom_order) {
  dat_chr <- dat[dat$chr == ch & dat$k %in% chr_k, ]
  if (nrow(dat_chr) == 0) {
    cat(sprintf("  chr%s: no data, skipping\n", ch))
    next
  }
  dat_chr$k_label   <- factor(dat_chr$k_label,   levels = paste0("k = ", chr_k))
  dat_chr$maf_label <- factor(dat_chr$maf_label, levels = maf_labels)

  p <- make_histogram(dat_chr, base_size = 10) +
    ggtitle(paste0("Chromosome ", ch)) +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
  print(p)
  cat(sprintf("  chr%s: %d windows\n", ch, nrow(dat_chr)))
}

dev.off()
cat(sprintf("  Saved: %s\n", chr_file))

cat("\nDone.\n")
