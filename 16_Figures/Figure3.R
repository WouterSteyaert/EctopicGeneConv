#!/usr/bin/env Rscript
# ==============================================================================
# Figure 3: Sliding window enrichment
# ==============================================================================
# Panels:
#   (a) 3x3 histogram grid: log2(OR) distribution
#       k = {21, 51, 91} x MAF = {<1e-5, 0.01-0.05, 0.1-0.5}
#   (b) Ideogram: unique peaks across all chromosomes, colored by segdup
#
# Data sources:
#   (a) stats_allmapp/{k}/ — non-overlapping 1Mb windows
#   (b) sliding_window_analysis/peaks_summary.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

# ---- Paths ----
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
STATS_DIR <- file.path(project_root, "geneconv_complete/stats_allmapp")
SW_DIR    <- file.path(project_root, "geneconv_complete/sliding_window_corrected")
FIG_DIR   <- file.path(project_root, "geneconv_complete/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- GRCh38 chromosome info ----
chrom_sizes <- c(
  "1"=248956422, "2"=242193529, "3"=198295559, "4"=190214555,
  "5"=181538259, "6"=170805979, "7"=159345973, "8"=145138636,
  "9"=138394717, "10"=133797422, "11"=135086622, "12"=133275309,
  "13"=114364328, "14"=107043718, "15"=101991189, "16"=90338345,
  "17"=83257441, "18"=80373285, "19"=58617616, "20"=64444167,
  "21"=46709983, "22"=50818468, "X"=156040895, "Y"=57227415
)

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

chrom_order <- c(as.character(1:22), "X", "Y")

# ===========================================================================
# Panel (a): 3x3 histogram grid
# ===========================================================================
cat("=== Panel (a): 3x3 histogram ===\n")

main_k   <- c(21, 51, 91)
main_maf <- c("0_1e-05", "001_005", "01_05")
maf_labels <- c(
  "0_1e-05" = "AF < 1e-5",
  "001_005"  = "AF 0.01 - 0.05",
  "01_05"    = "AF 0.1 - 0.5"
)

autosomes <- names(chrom_sizes)[1:22]

# Read corrected enrichment data (non-overlapping 1Mb windows only)
corrected_all <- read.delim(file.path(SW_DIR, "enrichment_corrected_combined.tsv"),
                            header = TRUE, stringsAsFactors = FALSE)
# Non-overlapping windows only
corrected_all <- corrected_all[(corrected_all$window_start - 1) %% 1000000 == 0, ]
# Autosomes only
corrected_all <- corrected_all[corrected_all$chr %in% autosomes, ]

hist_data <- list()
for (k in main_k) {
  for (maf in main_maf) {
    df <- corrected_all[corrected_all$k == k & corrected_all$maf_bin == maf &
                        is.finite(corrected_all$log2_or), ]
    if (nrow(df) > 0) {
      hist_data[[length(hist_data) + 1]] <- data.frame(
        log2_or = df$log2_or, k = k, maf = maf, stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("  k = %d done\n", k))
}

hdat <- do.call(rbind, hist_data)
hdat$k_label   <- factor(paste0("k = ", hdat$k), levels = paste0("k = ", main_k))
hdat$maf_label <- factor(maf_labels[hdat$maf], levels = maf_labels)

cat(sprintf("  Windows: %d\n", nrow(hdat)))

medians <- aggregate(log2_or ~ k_label + maf_label, data = hdat, FUN = median)

p_a <- ggplot(hdat, aes(x = log2_or)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = NA, alpha = 0.7) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = "dashed") +
  geom_vline(data = medians, aes(xintercept = log2_or),
             colour = "red", linewidth = 0.5) +
  facet_grid(k_label ~ maf_label, scales = "free_y") +
  labs(x = expression(log[2]~"odds ratio"),
       y = "Number of 1-Mb windows") +
  ggtitle("a") +
  theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0),
    strip.text       = element_text(size = 9, face = "bold"),
    axis.text        = element_text(size = 7),
    axis.title       = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(5, 10, 5, 5)
  )

# ===========================================================================
# Panel (b): Ideogram of sliding window peaks
# ===========================================================================
cat("\n=== Panel (b): Ideogram ===\n")

peaks <- read.delim(file.path(SW_DIR, "peaks_corrected.tsv"),
                    header = TRUE, stringsAsFactors = FALSE)
# Normalize column names (corrected peaks use capitalized names)
names(peaks) <- tolower(names(peaks))
names(peaks) <- gsub("\\.", "_", names(peaks))
cat(sprintf("Total peak calls: %d\n", nrow(peaks)))

peaks$chr <- as.character(peaks$chromosome)

# Collapse the per-(k, AF) peak calls into DISTINCT GENOMIC REGIONS by an
# interval-merge across strata: overlapping or contiguous (adjacent on the
# 1-Mb grid) peak intervals on the same chromosome merge into a single region.
# A naive group-by-exact-coordinate over-counts loci that recur with slightly
# different boundaries across strata (it would report 262 rather than the
# 193 distinct regions).
classify_region <- function(chr, mid) {
  cd <- if (mid >= centromere_start[chr] && mid <= centromere_end[chr]) 0
        else min(abs(mid - centromere_start[chr]), abs(mid - centromere_end[chr]))
  td <- min(mid, chrom_sizes[chr] - mid)
  if (cd < 5e6) "pericentromeric" else if (td < 5e6) "subtelomeric" else "interstitial"
}
peaks <- peaks[order(match(peaks$chr, chrom_order), peaks$peak_start), ]
loc_list <- list()
for (cc in unique(peaks$chr)) {
  pc <- peaks[peaks$chr == cc, ]
  flush <- function(s, e, idx) {
    mid    <- (s + e) / 2
    sd_set <- pc$segdup_overlap[idx]
    seg    <- if (any(sd_set == "segdup")) "segdup" else
              if (any(sd_set == "partial")) "partial" else "none"
    data.frame(chr = cc, peak_start = s, peak_end = e, mid = mid,
               n_strata = length(idx), max_log2_or = max(pc$max_log2_or[idx]),
               segdup = seg, region = classify_region(cc, mid),
               stringsAsFactors = FALSE)
  }
  cur_s <- pc$peak_start[1]; cur_e <- pc$peak_end[1]; mem <- 1
  if (nrow(pc) > 1) for (i in 2:nrow(pc)) {
    if (pc$peak_start[i] <= cur_e + 1) {           # overlapping or contiguous
      cur_e <- max(cur_e, pc$peak_end[i]); mem <- c(mem, i)
    } else {
      loc_list[[length(loc_list) + 1]] <- flush(cur_s, cur_e, mem)
      cur_s <- pc$peak_start[i]; cur_e <- pc$peak_end[i]; mem <- i
    }
  }
  loc_list[[length(loc_list) + 1]] <- flush(cur_s, cur_e, mem)
}
loc_agg <- do.call(rbind, loc_list)
rownames(loc_agg) <- NULL
cat(sprintf("Total peak calls: %d ; distinct genomic regions (interval-merged): %d\n",
            nrow(peaks), nrow(loc_agg)))
cat("  Region class:\n"); print(table(loc_agg$region))
cat(sprintf("  SD-overlapping regions: %d (%.0f%%); entirely outside SDs: %d (%.0f%%)\n",
            sum(loc_agg$segdup != "none"), 100 * mean(loc_agg$segdup != "none"),
            sum(loc_agg$segdup == "none"), 100 * mean(loc_agg$segdup == "none")))

chr_df <- data.frame(
  chr      = chrom_order,
  size     = chrom_sizes[chrom_order],
  centro_s = centromere_start[chrom_order],
  centro_e = centromere_end[chrom_order],
  y        = rev(seq_along(chrom_order)),
  stringsAsFactors = FALSE
)

loc_agg$y <- chr_df$y[match(loc_agg$chr, chr_df$chr)]
loc_agg <- loc_agg[!is.na(loc_agg$y), ]

loc_agg$segdup <- factor(loc_agg$segdup,
                          levels = c("segdup", "partial", "none"))

bar_h <- 0.3
tick_max <- 0.55
loc_agg$tick_h <- tick_max * pmin(loc_agg$n_strata, 13) / 13

loc_segdup   <- loc_agg[loc_agg$segdup != "none", ]
loc_nosegdup <- loc_agg[loc_agg$segdup == "none", ]

p_b <- ggplot() +
  geom_rect(data = chr_df,
            aes(xmin = 0, xmax = centro_s / 1e6,
                ymin = y - bar_h, ymax = y + bar_h),
            fill = "grey88", colour = "grey60", linewidth = 0.25) +
  geom_rect(data = chr_df,
            aes(xmin = centro_e / 1e6, xmax = size / 1e6,
                ymin = y - bar_h, ymax = y + bar_h),
            fill = "grey88", colour = "grey60", linewidth = 0.25) +
  geom_rect(data = chr_df,
            aes(xmin = centro_s / 1e6, xmax = centro_e / 1e6,
                ymin = y - bar_h * 0.4, ymax = y + bar_h * 0.4),
            fill = "grey55", colour = "grey60", linewidth = 0.25) +
  geom_segment(data = loc_segdup,
               aes(x    = mid / 1e6, xend = mid / 1e6,
                   y    = y + bar_h,
                   yend = y + bar_h + tick_h,
                   colour = segdup),
               linewidth = 0.9, lineend = "round") +
  geom_segment(data = loc_nosegdup,
               aes(x    = mid / 1e6, xend = mid / 1e6,
                   y    = y + bar_h,
                   yend = y + bar_h + tick_h,
                   colour = segdup),
               linewidth = 1.2, lineend = "round") +
  geom_segment(data = loc_nosegdup,
               aes(x    = mid / 1e6, xend = mid / 1e6,
                   y    = y - bar_h,
                   yend = y - bar_h - tick_h),
               colour = "#E41A1C", linewidth = 1.2, lineend = "round",
               show.legend = FALSE) +
  scale_colour_manual(
    values = c("segdup"  = "#93C47D",
               "partial" = "#F6B26B",
               "none"    = "#E41A1C"),
    labels = c("Segdup (>50%)   ",
               "Partial (<50%)   ",
               "No segdup"),
    name = NULL,
    drop = FALSE,
    guide = guide_legend(override.aes = list(linewidth = 3))
  ) +
  scale_y_continuous(
    breaks = chr_df$y,
    labels = chr_df$chr,
    expand = c(0.03, 0)
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = "Position (Mb)", y = NULL) +
  ggtitle("b") +
  theme_classic(base_size = 9) +
  theme(
    plot.title         = element_text(face = "bold", size = 11, hjust = 0),
    axis.text.y        = element_text(size = 7, face = "bold"),
    axis.text.x        = element_text(size = 7),
    axis.ticks.y       = element_blank(),
    legend.position    = "bottom",
    legend.key.size    = unit(0.5, "cm"),
    legend.text        = element_text(size = 7),
    legend.spacing.x   = unit(1.0, "cm"),
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(5, 10, 5, 5)
  )

# ===========================================================================
# Combined Figure 3
# ===========================================================================
cat("\n=== Combining panels ===\n")

fig3 <- plot_grid(p_a, p_b,
                  ncol = 1,
                  rel_heights = c(1, 2))

ggsave(file.path(FIG_DIR, "Figure3.pdf"), fig3,
       width = 10, height = 14, device = cairo_pdf)

# Also save panels separately
ggsave(file.path(FIG_DIR, "Figure3a_Histogram.pdf"), p_a,
       width = 9, height = 6, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure3b_Ideogram.pdf"), p_b,
       width = 10, height = 7.5, device = cairo_pdf)

cat(sprintf("\nFigure 3 saved to: %s\n", FIG_DIR))
cat("  Figure3.pdf             (combined)\n")
cat("  Figure3a_Histogram.pdf  (panel a)\n")
cat("  Figure3b_Ideogram.pdf   (panel b)\n")

# ===========================================================================
# Source data
# ===========================================================================
TSV_DIR <- file.path(FIG_DIR, "source_data")
dir.create(TSV_DIR, showWarnings = FALSE, recursive = TRUE)

# Panel (a): per-window log2(OR) values
df_a <- hdat[, c("k", "maf", "log2_or")]
write.table(df_a, file.path(TSV_DIR, "Figure3a_sliding_window_log2OR.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Panel (b): peak locations
df_b <- loc_agg[, c("chr", "peak_start", "peak_end", "n_strata",
                     "max_log2_or", "segdup", "region")]
write.table(df_b, file.path(TSV_DIR, "Figure3b_peak_locations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nSource data written to:\n")
cat("  source_data/Figure3a_sliding_window_log2OR.tsv\n")
cat("  source_data/Figure3b_peak_locations.tsv\n")
