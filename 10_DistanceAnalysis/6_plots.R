#!/usr/bin/env Rscript
# =============================================================================
# Step 6: Distance analysis plots
# =============================================================================
#
# Reads summary tables from step 5 and raw gnomad files for distributions.
# All plot filenames are numbered for easy ordering.
#
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
basedir <- file.path(project_root, "geneconv_complete/distance_analysis/summary")
rawdir  <- file.path(project_root, "geneconv_complete/distance_analysis")
plotdir <- file.path(basedir, "plots")
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

# Custom theme
theme_gc <- theme_bw(base_size = 12) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "right",
        strip.background = element_rect(fill = "grey90")
    )

# Distance bin order
dist_levels <- c("<1kb", "1-10kb", "10-100kb", "100kb-1Mb",
                 "1-10Mb", "10-100Mb", ">100Mb", "inter")

k_values <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)

# Load summary tables
df1 <- read.delim(file.path(basedir, "distance_af_by_k.tsv"))
df1$dist_bin <- factor(df1$dist_bin, levels = dist_levels)
df1$k <- factor(df1$k)

df1_intra <- df1 %>% filter(dist_bin != "inter")
df1_inter <- df1 %>% filter(dist_bin == "inter")

# =============================================================================
# 01: AF vs distance (intra)
# =============================================================================
cat("01: AF vs distance...\n")

p01 <- ggplot(df1_intra, aes(x = dist_bin, y = mean_AF, colour = k, group = k)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_colour_viridis_d(option = "turbo") +
    labs(title = "Mean allele frequency vs donor-acceptor distance",
         x = "Distance bin", y = "Mean AF (max of pair)",
         colour = "k-mer\nlength") +
    theme_gc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plotdir, "01_distance_vs_af.pdf"), p01, width = 10, height = 6)

# =============================================================================
# 02: % pairs with variant vs distance (intra)
# =============================================================================
cat("02: % pairs with variant...\n")

p02 <- ggplot(df1_intra, aes(x = dist_bin, y = 100 * n_any_var / n_total,
                              colour = k, group = k)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_colour_viridis_d(option = "turbo") +
    labs(title = "% pairs with GC-consistent gnomAD variant vs distance",
         x = "Distance bin", y = "% pairs with variant",
         colour = "k-mer\nlength") +
    theme_gc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plotdir, "02_distance_vs_pct_variant.pdf"), p02, width = 10, height = 6)

# =============================================================================
# 03: GC-bias (% W→S) vs distance (intra)
# =============================================================================
cat("03: GC-bias vs distance...\n")

p03 <- ggplot(df1_intra, aes(x = dist_bin, y = pct_WS, colour = k, group = k)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey50") +
    scale_colour_viridis_d(option = "turbo") +
    labs(title = "GC-bias (% W→S) vs distance",
         subtitle = "Dashed line = 50% (no bias)",
         x = "Distance bin", y = "% W→S transitions",
         colour = "k-mer\nlength") +
    theme_gc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plotdir, "03_distance_vs_gc_bias.pdf"), p03, width = 10, height = 6)

# =============================================================================
# 04: Intra vs inter comparison
# =============================================================================
cat("04: Intra vs inter...\n")

# Combine intra (all bins aggregated) vs inter
df1_intra_agg <- df1_intra %>%
    group_by(k) %>%
    summarise(
        type = "intra",
        mean_AF = weighted.mean(mean_AF, pmax(n_any_var, 1)),
        pct_WS = weighted.mean(pct_WS, pmax(n_any_var + n_both_var, 1), na.rm = TRUE),
        n_total = sum(n_total),
        n_any_var = sum(n_any_var),
        n_both_var = sum(n_both_var),
        .groups = "drop"
    ) %>%
    mutate(pct_var = 100 * n_any_var / n_total)

df1_inter_agg <- df1_inter %>%
    mutate(
        type = "inter",
        pct_var = 100 * n_any_var / n_total
    ) %>%
    select(k, type, n_total, n_any_var, n_both_var, mean_AF, pct_WS, pct_var)

df_compare <- bind_rows(df1_intra_agg, df1_inter_agg)
df_compare$k_num <- as.numeric(as.character(df_compare$k))

# 04a: AF intra vs inter
p04a <- ggplot(df_compare, aes(x = k_num, y = mean_AF, colour = type)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    labs(title = "Mean AF: intra- vs inter-chromosomal pairs",
         x = "k-mer length", y = "Mean AF (max of pair)",
         colour = "Type") +
    scale_x_continuous(breaks = k_values) +
    scale_colour_manual(values = c("intra" = "#d95f02", "inter" = "#7570b3")) +
    theme_gc

ggsave(file.path(plotdir, "04a_intra_vs_inter_af.pdf"), p04a, width = 8, height = 5)

# 04b: % with variant intra vs inter
p04b <- ggplot(df_compare, aes(x = k_num, y = pct_var, colour = type)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    labs(title = "% pairs with variant: intra- vs inter-chromosomal",
         x = "k-mer length", y = "% pairs with GC-consistent variant",
         colour = "Type") +
    scale_x_continuous(breaks = k_values) +
    scale_colour_manual(values = c("intra" = "#d95f02", "inter" = "#7570b3")) +
    theme_gc

ggsave(file.path(plotdir, "04b_intra_vs_inter_pct_variant.pdf"), p04b, width = 8, height = 5)

# 04c: GC-bias intra vs inter
p04c <- ggplot(df_compare, aes(x = k_num, y = pct_WS, colour = type)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey50") +
    labs(title = "GC-bias (% W→S): intra- vs inter-chromosomal",
         x = "k-mer length", y = "% W→S transitions",
         colour = "Type") +
    scale_x_continuous(breaks = k_values) +
    scale_colour_manual(values = c("intra" = "#d95f02", "inter" = "#7570b3")) +
    theme_gc

ggsave(file.path(plotdir, "04c_intra_vs_inter_gc_bias.pdf"), p04c, width = 8, height = 5)

# =============================================================================
# 05: Pair counts (total, with variant, without)
# =============================================================================
cat("05: Pair counts...\n")

df_counts <- df1 %>%
    group_by(k) %>%
    summarise(
        total = sum(n_total),
        any_variant = sum(n_any_var),
        both_variant = sum(n_both_var),
        .groups = "drop"
    ) %>%
    mutate(
        no_variant = total - any_variant,
        one_variant = any_variant - both_variant
    )
df_counts$k_num <- as.numeric(as.character(df_counts$k))

df_counts_long <- df_counts %>%
    select(k, k_num, no_variant, one_variant, both_variant) %>%
    pivot_longer(cols = c(no_variant, one_variant, both_variant),
                 names_to = "category", values_to = "count") %>%
    mutate(category = factor(category,
        levels = c("both_variant", "one_variant", "no_variant"),
        labels = c("Both positions", "One position", "Neither")))

p05a <- ggplot(df_counts_long, aes(x = k, y = count / 1e6, fill = category)) +
    geom_col(position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = c("Both positions" = "#2166ac",
                                  "One position" = "#92c5de",
                                  "Neither" = "#f4a582")) +
    scale_y_log10(labels = scales::comma) +
    labs(title = "Pair counts by gnomAD variant status",
         x = "k-mer length", y = "Number of pairs (millions, log scale)",
         fill = "GC-consistent\nvariant found") +
    theme_gc

ggsave(file.path(plotdir, "05a_pair_counts_absolute.pdf"), p05a, width = 10, height = 6)

# Proportional
p05b <- ggplot(df_counts_long, aes(x = k, y = count, fill = category)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("Both positions" = "#2166ac",
                                  "One position" = "#92c5de",
                                  "Neither" = "#f4a582")) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Pair counts by gnomAD variant status (proportional)",
         x = "k-mer length", y = "Proportion",
         fill = "GC-consistent\nvariant found") +
    theme_gc

ggsave(file.path(plotdir, "05b_pair_counts_proportional.pdf"), p05b, width = 10, height = 6)

# =============================================================================
# 06: AF distribution per k (from raw gnomad files)
# =============================================================================
cat("06: AF distribution...\n")

# Read a sample of AF values from each k (subsample for large files)
af_data <- data.frame()
for (kk in k_values) {
    infile <- file.path(rawdir, kk, "distances.diff.gnomad.txt.gz")
    if (!file.exists(infile)) next
    cat("  Reading k=", kk, "...\n")

    # Use awk to extract non-zero AFs (max of pair), pipe through head for large files
    cmd <- sprintf(
        "zcat %s | awk -F'\\t' '{a1=$9+0; a2=$12+0; if(a1>0||a2>0){m=(a1>a2)?a1:a2; print m}}' | shuf -n 100000 2>/dev/null || zcat %s | awk -F'\\t' '{a1=$9+0; a2=$12+0; if(a1>0||a2>0){m=(a1>a2)?a1:a2; print m}}' | head -100000",
        infile, infile
    )
    afs <- tryCatch(
        as.numeric(readLines(pipe(cmd))),
        error = function(e) numeric(0)
    )
    if (length(afs) > 0) {
        af_data <- bind_rows(af_data, data.frame(k = kk, AF = afs))
    }
}

if (nrow(af_data) > 0) {
    af_data$k <- factor(af_data$k)

    # 06a: Density plot
    p06a <- ggplot(af_data, aes(x = AF, colour = k)) +
        geom_density(linewidth = 0.6) +
        scale_colour_viridis_d(option = "turbo") +
        scale_x_log10(labels = scales::comma) +
        labs(title = "AF distribution of GC-consistent variants (log scale)",
             subtitle = "Max AF per pair, up to 100K sampled per k",
             x = "Allele frequency", y = "Density",
             colour = "k-mer\nlength") +
        theme_gc

    ggsave(file.path(plotdir, "06a_af_distribution_density.pdf"), p06a, width = 10, height = 6)

    # 06b: AF binned into categories
    af_data <- af_data %>%
        mutate(af_bin = cut(AF,
            breaks = c(0, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0),
            labels = c("<0.1%", "0.1-1%", "1-5%", "5-10%", "10-50%", ">50%"),
            include.lowest = TRUE))

    p06b <- ggplot(af_data, aes(x = k, fill = af_bin)) +
        geom_bar(position = "fill") +
        scale_fill_brewer(palette = "YlOrRd", direction = 1) +
        scale_y_continuous(labels = scales::percent_format()) +
        labs(title = "AF distribution by category per k",
             x = "k-mer length", y = "Proportion",
             fill = "AF bin") +
        theme_gc

    ggsave(file.path(plotdir, "06b_af_distribution_binned.pdf"), p06b, width = 10, height = 6)
}

# =============================================================================
# 07: Transition spectrum by k
# =============================================================================
cat("07: Transition spectrum...\n")

df2 <- read.delim(file.path(basedir, "transitions_by_k.tsv"))
df2$k <- factor(df2$k)

ws_trans <- c("A>G", "A>C", "T>G", "T>C")
sw_trans <- c("G>A", "C>A", "G>T", "C>T")
ss_trans <- c("C>G", "G>C")
ww_trans <- c("A>T", "T>A")
trans_order <- c(ws_trans, sw_trans, ss_trans, ww_trans)

df2$transition <- factor(df2$transition, levels = trans_order)
df2 <- df2 %>%
    group_by(k) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

p07 <- ggplot(df2, aes(x = transition, y = pct, fill = direction)) +
    geom_col() +
    facet_wrap(~ k, ncol = 5) +
    scale_fill_manual(values = c("WS" = "#2166ac", "SW" = "#b2182b",
                                 "SS" = "#762a83", "WW" = "#1b7837")) +
    labs(title = "Transition spectrum per k-mer length",
         x = "Nucleotide transition", y = "% of all transitions",
         fill = "Direction") +
    theme_gc +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))

ggsave(file.path(plotdir, "07_transition_spectrum.pdf"), p07, width = 14, height = 8)

# =============================================================================
# 08: Complementary transition pairs
# =============================================================================
cat("08: Complementary pairs...\n")

df3 <- read.delim(file.path(basedir, "complementary_pairs.tsv"))
df3$k <- as.numeric(df3$k)

p08 <- ggplot(df3, aes(x = k, y = ratio_fwd_rev, colour = pair)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    labs(title = "Complementary transition ratios across k",
         subtitle = "Forward/Reverse count ratio (dashed = 1.0 = symmetric)",
         x = "k-mer length", y = "Forward / Reverse ratio",
         colour = "Transition\npair") +
    scale_x_continuous(breaks = unique(df3$k)) +
    theme_gc

ggsave(file.path(plotdir, "08_complementary_pairs.pdf"), p08, width = 10, height = 6)

# =============================================================================
# 09: GC-bias (W→S / S→W) by k
# =============================================================================
cat("09: GC-bias by k...\n")

df4 <- df2 %>%
    group_by(k, direction) %>%
    summarise(n_total = sum(n), .groups = "drop") %>%
    pivot_wider(names_from = direction, values_from = n_total, values_fill = 0) %>%
    mutate(
        WS_SW_ratio = WS / SW,
        pct_WS = 100 * WS / (WS + SW + SS + WW)
    )
df4$k_num <- as.numeric(as.character(df4$k))

p09 <- ggplot(df4, aes(x = k_num, y = WS_SW_ratio)) +
    geom_line(linewidth = 0.8, colour = "#2166ac") +
    geom_point(size = 3, colour = "#2166ac") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    labs(title = "GC-bias: W→S / S→W ratio by k-mer length",
         subtitle = "Ratio > 1 indicates GC-biased gene conversion",
         x = "k-mer length", y = "W→S / S→W ratio") +
    scale_x_continuous(breaks = unique(df4$k_num)) +
    theme_gc

ggsave(file.path(plotdir, "09_gc_bias_by_k.pdf"), p09, width = 8, height = 5)

# =============================================================================
# 10: Symmetry (|AF1 - AF2|) by distance and k
# =============================================================================
cat("10: Symmetry...\n")

df5 <- read.delim(file.path(basedir, "symmetry_by_k.tsv"))
df5$dist_bin <- factor(df5$dist_bin, levels = dist_levels)
df5$k <- factor(df5$k)

df5_intra <- df5 %>% filter(dist_bin != "inter")

p10 <- ggplot(df5_intra, aes(x = dist_bin, y = mean_abs_diff,
                              colour = k, group = k)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_colour_viridis_d(option = "turbo") +
    labs(title = "AF symmetry vs distance",
         subtitle = "Lower |AF1-AF2| = more symmetric",
         x = "Distance bin", y = "Mean |AF1 - AF2|",
         colour = "k-mer\nlength") +
    theme_gc +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plotdir, "10_symmetry_by_distance.pdf"), p10, width = 10, height = 6)

cat("\nAll plots written to:", plotdir, "\n")
cat("Files:\n")
system(paste("ls -lh", plotdir))
