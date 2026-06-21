#!/usr/bin/env Rscript
# =============================================================================
# Step 7b: gBGC (GC-biased gene conversion) plots
# =============================================================================
#
# Reads summary tables from 7_gbgc_summary.sh and produces publication-quality
# plots showing the gBGC signal in allele frequency, not counts.
#
# Key insight: WS/SW count ratio is always <1 due to mutation bias (~2x more
# S→W mutations). The real gBGC signal is that WS variants reach 2-3x higher
# allele frequency than SW variants.
#
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
basedir <- file.path(project_root, "geneconv_complete/distance_analysis/summary")
plotdir <- file.path(basedir, "plots")
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

# Custom theme (matching 6_plots.R)
theme_gc <- theme_bw(base_size = 12) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "right",
        strip.background = element_rect(fill = "grey90")
    )

dist_levels <- c("<1kb", "1-10kb", "10-100kb", "100kb-1Mb",
                 "1-10Mb", "10-100Mb", ">100Mb", "inter")

k_values <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)

# =============================================================================
# Load data
# =============================================================================
cat("Loading summary tables...\n")

df_af <- read.delim(file.path(basedir, "gbgc_af_by_k.tsv"))
df_af_dist <- read.delim(file.path(basedir, "gbgc_af_by_k_distance.tsv"))
df_prop <- read.delim(file.path(basedir, "gbgc_proportion_by_af.tsv"))
df_spec <- read.delim(file.path(basedir, "gbgc_af_spectrum.tsv"))

# Also load the existing transitions table for plot 15
df_trans <- read.delim(file.path(basedir, "transitions_by_k.tsv"))

# =============================================================================
# 11: Core gBGC figure — Mean AF of WS vs SW per k
# =============================================================================
cat("11: Mean AF by k...\n")

df_af$k <- as.numeric(df_af$k)

# Left panel: mean AF for WS and SW
p11a <- ggplot(df_af, aes(x = k, y = mean_AF, colour = direction)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("WS" = "#2166ac", "SW" = "#b2182b")) +
    scale_x_continuous(breaks = k_values) +
    labs(title = "Mean allele frequency",
         x = "k-mer length", y = "Mean AF",
         colour = "Direction") +
    theme_gc

# Right panel: AF ratio WS/SW
df_af_wide <- df_af %>%
    select(k, direction, mean_AF) %>%
    pivot_wider(names_from = direction, values_from = mean_AF) %>%
    mutate(AF_ratio = WS / SW)

p11b <- ggplot(df_af_wide, aes(x = k, y = AF_ratio)) +
    geom_line(linewidth = 0.9, colour = "#1a1a1a") +
    geom_point(size = 3, colour = "#1a1a1a") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = k_values) +
    labs(title = "WS / SW allele frequency ratio",
         subtitle = "Ratio > 1 = gBGC favoring GC alleles",
         x = "k-mer length", y = "AF ratio (WS / SW)") +
    theme_gc

title11 <- ggdraw() +
    draw_label("GC-biased gene conversion: allele frequency signal",
               fontface = "bold", size = 14, x = 0.02, hjust = 0)
p11 <- plot_grid(title11,
                 plot_grid(p11a, p11b, ncol = 2, labels = c("(a)", "(b)")),
                 ncol = 1, rel_heights = c(0.08, 1))

ggsave(file.path(plotdir, "11_gbgc_mean_af_by_k.pdf"), p11, width = 14, height = 5.5)

# =============================================================================
# 12: Rising %WS with AF threshold
# =============================================================================
cat("12: %WS by AF threshold...\n")

df_prop$k <- factor(df_prop$k)
df_prop$af_threshold <- as.numeric(df_prop$af_threshold)

# Baseline: %WS at lowest AF threshold (=0, all variants) per k
baseline <- df_prop %>%
    filter(af_threshold == 0) %>%
    summarise(baseline_pct = mean(pct_WS)) %>%
    pull(baseline_pct)

# Only plot thresholds > 0 on log scale, plus 0 as a separate point
df_prop_plot <- df_prop %>%
    filter(n_total >= 10)  # require minimum N

p12 <- ggplot(df_prop_plot %>% filter(af_threshold > 0),
              aes(x = af_threshold, y = pct_WS, colour = k, group = k)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    geom_hline(yintercept = baseline, linetype = "dashed", colour = "grey40",
               linewidth = 0.5) +
    geom_hline(yintercept = 50, linetype = "dotted", colour = "grey60",
               linewidth = 0.5) +
    annotate("text", x = 0.3, y = baseline + 1.5, label = "Mutation baseline",
             colour = "grey40", size = 3, hjust = 1) +
    annotate("text", x = 0.3, y = 51.5, label = "No bias (50%)",
             colour = "grey60", size = 3, hjust = 1) +
    scale_x_log10(labels = scales::label_scientific()) +
    scale_colour_viridis_d(option = "turbo") +
    labs(title = "Fraction of W→S variants increases with allele frequency",
         subtitle = "gBGC pushes WS alleles to higher frequency → enrichment at high AF",
         x = "AF threshold (variants with AF ≥ threshold)",
         y = "% W→S among W→S + S→W",
         colour = "k-mer\nlength") +
    theme_gc

ggsave(file.path(plotdir, "12_gbgc_ws_proportion_by_af.pdf"), p12, width = 10, height = 6)

# =============================================================================
# 13: AF distribution comparison WS vs SW
# =============================================================================
cat("13: AF spectrum...\n")

af_bin_levels <- c("<1e-4", "1e-4_1e-3", "1e-3_0.01",
                   "0.01_0.05", "0.05_0.1", "0.1_0.5", ">=0.5")
af_bin_labels <- c("<0.01%", "0.01-0.1%", "0.1-1%",
                   "1-5%", "5-10%", "10-50%", "≥50%")

df_spec$af_bin <- factor(df_spec$af_bin, levels = af_bin_levels, labels = af_bin_labels)
df_spec$k <- as.numeric(df_spec$k)

# Subset of k values for clarity
k_subset <- c(17, 31, 51, 71, 91)
df_spec_sub <- df_spec %>% filter(k %in% k_subset)
df_spec_sub$k <- factor(df_spec_sub$k)

p13 <- ggplot(df_spec_sub, aes(x = af_bin, y = pct_of_direction,
                                fill = direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ k, ncol = 5, labeller = labeller(k = function(x) paste0("k = ", x))) +
    scale_fill_manual(values = c("WS" = "#2166ac", "SW" = "#b2182b")) +
    labs(title = "Allele frequency spectrum: W→S vs S→W variants",
         subtitle = "WS distribution shifted toward higher AF = gBGC signature",
         x = "Allele frequency bin", y = "% of variants in direction",
         fill = "Direction") +
    theme_gc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(plotdir, "13_gbgc_af_spectrum.pdf"), p13, width = 16, height = 5)

# =============================================================================
# 14: gBGC signal vs paralog distance
# =============================================================================
cat("14: gBGC by distance...\n")

df_af_dist$dist_bin <- factor(df_af_dist$dist_bin, levels = dist_levels)
df_af_dist$k <- as.numeric(df_af_dist$k)

# Compute AF ratio per k × distance bin
df_ratio_dist <- df_af_dist %>%
    select(k, dist_bin, direction, mean_AF) %>%
    pivot_wider(names_from = direction, values_from = mean_AF) %>%
    filter(!is.na(WS) & !is.na(SW) & SW > 0) %>%
    mutate(AF_ratio = WS / SW)

df_ratio_dist$k <- factor(df_ratio_dist$k)

p14 <- ggplot(df_ratio_dist %>% filter(dist_bin != "inter"),
              aes(x = dist_bin, y = AF_ratio, colour = k, group = k)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_colour_viridis_d(option = "turbo") +
    labs(title = "gBGC signal (WS/SW AF ratio) vs paralog distance",
         subtitle = "Ratio > 1 indicates GC-biased gene conversion",
         x = "Paralog distance bin", y = "Mean AF ratio (WS / SW)",
         colour = "k-mer\nlength") +
    theme_gc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plotdir, "14_gbgc_by_distance.pdf"), p14, width = 10, height = 6)

# =============================================================================
# 15: Mutation-bias corrected WS/SW count ratio
# =============================================================================
cat("15: Corrected WS/SW ratio...\n")

# Get WS and SW counts from transitions table
df_wsw <- df_trans %>%
    filter(direction %in% c("WS", "SW")) %>%
    group_by(k, direction) %>%
    summarise(n_total = sum(n), .groups = "drop") %>%
    pivot_wider(names_from = direction, values_from = n_total) %>%
    mutate(
        raw_ratio = WS / SW,
        pct_WS = 100 * WS / (WS + SW)
    )
df_wsw$k <- as.numeric(as.character(df_wsw$k))

# Baseline: use the proportion table at af_threshold = 0 and lowest k (k=17)
# This represents the mutation-driven baseline before gBGC selection
baseline_ratio <- df_wsw %>%
    filter(k == min(k)) %>%
    pull(raw_ratio)

df_wsw <- df_wsw %>%
    mutate(corrected_ratio = raw_ratio / baseline_ratio)

# Two-panel plot
p15a <- ggplot(df_wsw, aes(x = k, y = raw_ratio)) +
    geom_line(linewidth = 0.9, colour = "#b2182b") +
    geom_point(size = 3, colour = "#b2182b") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = k_values) +
    labs(title = "Raw WS/SW count ratio",
         subtitle = "Always < 1 due to ~2x higher S→W mutation rate",
         x = "k-mer length", y = "WS / SW count ratio") +
    theme_gc

p15b <- ggplot(df_wsw, aes(x = k, y = corrected_ratio)) +
    geom_line(linewidth = 0.9, colour = "#2166ac") +
    geom_point(size = 3, colour = "#2166ac") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = k_values) +
    labs(title = "Corrected WS/SW count ratio",
         subtitle = sprintf("Normalized by k=17 baseline (%.3f)", baseline_ratio),
         x = "k-mer length", y = "Corrected ratio (raw / baseline)") +
    theme_gc

title15 <- ggdraw() +
    draw_label("Mutation-bias correction of WS/SW count ratio",
               fontface = "bold", size = 14, x = 0.02, hjust = 0)
p15 <- plot_grid(title15,
                 plot_grid(p15a, p15b, ncol = 2, labels = c("A", "B")),
                 ncol = 1, rel_heights = c(0.08, 1))

ggsave(file.path(plotdir, "15_gbgc_corrected_ratio.pdf"), p15, width = 14, height = 5.5)

# =============================================================================
# Done
# =============================================================================
cat("\nAll gBGC plots written to:", plotdir, "\n")
cat("Files:\n")
system(paste("ls -lh", file.path(plotdir, "1[1-5]_*.pdf")))
