#!/usr/bin/env Rscript
# ==============================================================================
# Supplementary Figure 2: per-chromosome sliding-window enrichment
#
# 24 pages (one per chromosome), each a 10 x 9 facet grid:
#   rows    = 10 template lengths k = 17, 19, 21, 31, 41, 51, 61, 71, 81, 91
#   columns = 9 allele frequency bins
# Each panel: blue line of log2(OR) vs genomic position (1-Mb sliding window,
# 50-kb step), with a horizontal black reference at log2(OR) = 0 and a red
# dashed vertical line at the GRCh38 centromere start.
#
# Input:
#   <PROJECT_ROOT>/geneconv_complete/sliding_window_analysis/
#     enrichment_allmapp_k{K}_maf{MAF}.tsv
#
# Output:
#   <PROJECT_ROOT>/geneconv_complete/figures/SupplementaryFigure2_KxMAF.pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- Paths ----
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
SW_DIR  <- file.path(project_root, "geneconv_complete/sliding_window_analysis")
OUT_DIR <- file.path(project_root, "geneconv_complete/figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- GRCh38 centromere boundaries ----
centromere_start <- c(
  "1"=121700000, "2"=91800000, "3"=87800000, "4"=48200000,
  "5"=46100000, "6"=58500000, "7"=58100000, "8"=43100000,
  "9"=42200000, "10"=38000000, "11"=51000000, "12"=33200000,
  "13"=16500000, "14"=16100000, "15"=17500000, "16"=35300000,
  "17"=22700000, "18"=15400000, "19"=24200000, "20"=26200000,
  "21"=10900000, "22"=13700000, "X"=58100000, "Y"=10300000
)
centromere_end <- c(
  "1"=125100000, "2"=96000000, "3"=93900000, "4"=51800000,
  "5"=51400000, "6"=62600000, "7"=62100000, "8"=47200000,
  "9"=45500000, "10"=41600000, "11"=55800000, "12"=37800000,
  "13"=18900000, "14"=18200000, "15"=20500000, "16"=38400000,
  "17"=27400000, "18"=21500000, "19"=28100000, "20"=29900000,
  "21"=13000000, "22"=17400000, "X"=63800000, "Y"=10600000
)

chrom_sizes <- c(
  "1"=248956422, "2"=242193529, "3"=198295559, "4"=190214555,
  "5"=181538259, "6"=170805979, "7"=159345973, "8"=145138636,
  "9"=138394717, "10"=133797422, "11"=135086622, "12"=133275309,
  "13"=114364328, "14"=107043718, "15"=101991189, "16"=90338345,
  "17"=83257441, "18"=80373285, "19"=58617616, "20"=64444167,
  "21"=46709983, "22"=50818468, "X"=156040895, "Y"=57227415
)

chrom_order <- c(as.character(1:22), "X", "Y")

# ---- Selected k values and MAF bins ----
k_select   <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)
maf_select <- c("0_1e-05", "1e-05_00001", "00001_0001", "0001_0005",
                "0005_001", "001_005", "005_01", "01_05", "05_2")

maf_labels <- c(
  "0_1e-05"     = "<1e-5",
  "1e-05_00001" = "1e-5 - 1e-4",
  "00001_0001"  = "1e-4 - 1e-3",
  "0001_0005"   = "1e-3 - 5e-3",
  "0005_001"    = "5e-3 - 0.01",
  "001_005"     = "0.01 - 0.05",
  "005_01"      = "0.05 - 0.1",
  "01_05"       = "0.1 - 0.5",
  "05_2"        = ">0.5"
)

# ---- Read selected files ----
cat("=== Reading enrichment files ===\n")

all_data <- list()
for (k_val in k_select) {
  for (maf_bin in maf_select) {
    infile <- file.path(SW_DIR,
                        sprintf("enrichment_allmapp_k%d_maf%s.tsv", k_val, maf_bin))
    if (!file.exists(infile)) {
      cat("  MISSING:", basename(infile), "\n")
      next
    }
    cat(sprintf("  k=%d MAF=%s\n", k_val, maf_bin))
    d <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE)
    d$chr <- as.character(d$chr)
    d <- d[, c("chr", "window_mid", "log2_or", "region_type")]
    d$k   <- k_val
    d$maf <- maf_bin
    all_data <- c(all_data, list(d))
  }
}

dat <- do.call(rbind, all_data)
rm(all_data)
gc(verbose = FALSE)

# Facet labels
dat$k_label   <- factor(paste0("k = ", dat$k),
                         levels = paste0("k = ", k_select))
dat$maf_label <- factor(maf_labels[dat$maf],
                         levels = maf_labels[maf_select])

cat(sprintf("Total rows: %d\n\n", nrow(dat)))

# ---- Generate one PDF with 24 pages ----
pdf_file <- file.path(OUT_DIR, "SupplementaryFigure2_KxMAF.pdf")
cat(sprintf("Writing: %s\n", pdf_file))
pdf(pdf_file, width = 16, height = 11)

for (chr_name in chrom_order) {

  cat(sprintf("  chr%s\n", chr_name))

  d_chr <- dat[dat$chr == chr_name, ]
  if (nrow(d_chr) == 0) next

  chr_size <- chrom_sizes[chr_name]
  centro_s <- centromere_start[chr_name]
  centro_e <- centromere_end[chr_name]

  # Centromere annotation needs facet-compatible data frame
  facet_grid <- expand.grid(
    k_label   = levels(dat$k_label),
    maf_label = levels(dat$maf_label),
    stringsAsFactors = FALSE
  )
  facet_grid$k_label   <- factor(facet_grid$k_label, levels = levels(dat$k_label))
  facet_grid$maf_label <- factor(facet_grid$maf_label, levels = levels(dat$maf_label))

  # Sort within panel so geom_line connects in genomic order
  d_chr <- d_chr[order(d_chr$k, d_chr$maf, d_chr$window_mid), ]

  p <- ggplot(d_chr, aes(x = window_mid / 1e6, y = log2_or)) +
    # Reference
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    # Centromere line
    geom_vline(xintercept = centro_s / 1e6,
               colour = "red", alpha = 0.5, linetype = "dashed", linewidth = 0.25) +
    # Line plot (blue, thin)
    geom_line(colour = "steelblue", linewidth = 0.25, na.rm = TRUE) +
    # Facet grid: rows = k, columns = MAF
    facet_grid(k_label ~ maf_label, scales = "free_y") +
    scale_x_continuous(limits = c(0, chr_size / 1e6), expand = c(0.01, 0)) +
    labs(
      title = sprintf("Chromosome %s", chr_name),
      x = "Position (Mb)",
      y = expression(log[2]~"odds ratio")
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
      strip.text       = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "grey80"),
      axis.text        = element_text(size = 6.5, color = "black"),
      axis.title       = element_text(size = 9),
      panel.spacing    = unit(0.4, "lines"),
      plot.margin      = margin(8, 8, 5, 5)
    )

  print(p)
}

dev.off()
cat(sprintf("\nDone! 24 pages: %s\n", pdf_file))
