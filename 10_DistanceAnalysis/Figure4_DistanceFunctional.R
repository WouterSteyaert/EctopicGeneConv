#!/usr/bin/env Rscript
# ==============================================================================
# Figure 4: Distance analysis
# ==============================================================================
# Panels:
#   (a) Mean AF vs distance between GC pairs (per k)
#   (b) Intra- vs inter-chromosomal mean AF comparison
#   (c) Predicted vs observed AF scatter (leave-one-k-out CV)
#
# Supplementary: Heatmap (k x distance -> mean AF)
#
# Reads: distance_analysis/summary/distance_af_by_k.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  if (!requireNamespace("cowplot", quietly = TRUE))
    install.packages("cowplot", repos = "https://cloud.r-project.org")
  library(cowplot)
})

# ---- Paths ----
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
DATA_DIR <- file.path(project_root, "geneconv_complete")
FIG_DIR  <- file.path(project_root, "geneconv_complete/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Theme ----
theme_fig <- theme_classic(base_size = 9) +
  theme(
    text         = element_text(family = "sans"),
    axis.text    = element_text(color = "black"),
    axis.title   = element_text(size = 9),
    plot.title   = element_text(face = "bold", size = 11, hjust = 0),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    plot.margin  = margin(5, 10, 5, 5)
  )

# ---- k-value color palette ----
k_values <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)
k_colors <- colorRampPalette(c("#9ECAE1", "#2171B5", "#08306B"))(10)
names(k_colors) <- as.character(k_values)

# ===========================================================================
# Read distance data
# ===========================================================================
cat("Reading distance data...\n")
dist_df <- read.delim(file.path(DATA_DIR, "distance_analysis/summary/distance_af_by_k.tsv"),
                      header = TRUE, stringsAsFactors = FALSE)

# Distance bin ordering
dist_order <- c("<1kb", "1-10kb", "10-100kb", "100kb-1Mb",
                "1-10Mb", "10-100Mb", ">100Mb", "inter")
dist_labels <- c("<1kb", "1\u201310kb", "10\u2013100kb", "100kb\u20131Mb",
                 "1\u201310Mb", "10\u2013100Mb", ">100Mb", "Inter-chr")

dist_df$dist_bin <- factor(dist_df$dist_bin, levels = dist_order, labels = dist_labels)
dist_df$k <- factor(dist_df$k, levels = k_values)

# ===========================================================================
# Panel (a): Mean AF vs distance
# ===========================================================================
# Intra-chromosomal only for main distance plot
intra <- dist_df[dist_df$dist_bin != "Inter-chr", ]

p_a <- ggplot(intra, aes(x = dist_bin, y = mean_AF, colour = k, group = k)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = k_colors, name = "k") +
  scale_y_continuous(limits = c(0, 0.13), expand = c(0, 0)) +
  labs(x = "Distance between GC pairs",
       y = "Mean allele frequency") +
  ggtitle("a") +
  theme_fig +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = c(0.85, 0.75),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))

# ===========================================================================
# Panel (b): Intra- vs inter-chromosomal mean AF
# ===========================================================================
# Compare mean AF for intra-chromosomal (all bins combined) vs inter-chromosomal
intra_all <- dist_df[dist_df$dist_bin != "Inter-chr", ]
inter_all <- dist_df[dist_df$dist_bin == "Inter-chr", ]

# Weighted mean AF across intra bins (weighted by n_any_var: pairs with at
# least one variant. n_total counts all pairs including those with zero
# variants, which would dilute the mean inappropriately.)
intra_agg <- aggregate(cbind(n_any_var, af_sum = mean_AF * n_any_var) ~ k,
                       data = intra_all, FUN = sum)
intra_agg$mean_AF <- intra_agg$af_sum / intra_agg$n_any_var

inter_agg <- inter_all[, c("k", "mean_AF")]

chr_comp <- rbind(
  data.frame(k = intra_agg$k, mean_AF = intra_agg$mean_AF,
             type = "Intra-chromosomal"),
  data.frame(k = inter_agg$k, mean_AF = inter_agg$mean_AF,
             type = "Inter-chromosomal")
)
chr_comp$type <- factor(chr_comp$type,
  levels = c("Intra-chromosomal", "Inter-chromosomal"))

p_b <- ggplot(chr_comp, aes(x = factor(k, levels = k_values), y = mean_AF,
                             fill = type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Intra-chromosomal" = "#E41A1C",
                                "Inter-chromosomal" = "#377EB8"),
                    name = NULL) +
  labs(x = "k-mer length", y = "Mean allele frequency") +
  ggtitle("b") +
  theme_fig +
  theme(legend.position = c(0.65, 0.90),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size = unit(0.35, "cm"),
        axis.text.x = element_text(size = 7))

# ===========================================================================
# Panel (c): Predicted vs Observed AF scatter (leave-one-k-out CV)
# ===========================================================================
cat("Panel c: Predicted vs Observed AF (leave-one-k-out CV)...\n")

# Map distance bins to approximate log-midpoints
dist_midpoints <- setNames(
  c(500, 5000, 50000, 500000, 5e6, 5e7, 2e8),
  c("<1kb", "1\u201310kb", "10\u2013100kb", "100kb\u20131Mb",
    "1\u201310Mb", "10\u2013100Mb", ">100Mb")
)
intra_model <- intra %>%
  mutate(
    k_num = as.numeric(as.character(k)),
    log_dist = log10(dist_midpoints[as.character(dist_bin)])
  )

# Leave-one-k-out cross-validation: for each k, fit on the other 9, predict that k
k_levels <- unique(intra_model$k_num)
intra_model$predicted_cv <- NA_real_

for (held_out in k_levels) {
  train <- intra_model[intra_model$k_num != held_out, ]
  test  <- intra_model[intra_model$k_num == held_out, ]
  fit_cv <- lm(mean_AF ~ k_num + log_dist, data = train)
  intra_model$predicted_cv[intra_model$k_num == held_out] <- predict(fit_cv, newdata = test)
}

# Cross-validated R²
ss_res <- sum((intra_model$mean_AF - intra_model$predicted_cv)^2)
ss_tot <- sum((intra_model$mean_AF - mean(intra_model$mean_AF))^2)
r2_cv <- 1 - ss_res / ss_tot

# Also report in-sample R² for comparison
fit_full <- lm(mean_AF ~ k_num + log_dist, data = intra_model)
r2_insample <- summary(fit_full)$r.squared

cat(sprintf("  In-sample R² = %.3f\n", r2_insample))
cat(sprintf("  Leave-one-k-out CV R² = %.3f\n", r2_cv))

p_c <- ggplot(intra_model, aes(x = predicted_cv, y = mean_AF, colour = k)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2) +
  scale_colour_manual(values = k_colors, name = "k") +
  labs(x = "Predicted mean AF\n(leave-one-k-out CV)",
       y = "Observed mean AF") +
  ggtitle(sprintf("c (CV R\u00B2 = %.2f)", r2_cv)) +
  coord_equal(xlim = c(0, 0.13), ylim = c(0, 0.13)) +
  theme_fig +
  theme(legend.position = c(0.85, 0.35),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))

# ===========================================================================
# Combine
# ===========================================================================
fig4 <- plot_grid(p_a, p_b, p_c,
                  ncol = 1, nrow = 3,
                  align = "v",
                  rel_heights = c(1, 1, 1.5))

ggsave(file.path(FIG_DIR, "Figure4_DistanceFunctional.pdf"), fig4,
       width = 6.5, height = 10.5, device = cairo_pdf)

ggsave(file.path(FIG_DIR, "Figure4a_distance_AF.pdf"), p_a,
       width = 4, height = 3.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure4b_intra_inter_AF.pdf"), p_b,
       width = 4, height = 3.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure4c_scatter_predicted_vs_observed.pdf"), p_c,
       width = 4, height = 3.5, device = cairo_pdf)

# ===========================================================================
# Supplementary: Heatmap of mean AF by k and distance
# ===========================================================================
cat("Supplementary: Heatmap...\n")
heat_df <- dist_df %>%
  mutate(k_num = as.numeric(as.character(k)))

p_heat <- ggplot(heat_df, aes(x = dist_bin, y = factor(k_num), fill = mean_AF)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", mean_AF * 100)),
            size = 2, colour = "black") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = "Mean AF",
                       labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Distance between GC pairs",
       y = "k-mer length") +
  theme_fig +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.3, "cm"))

ggsave(file.path(FIG_DIR, "FigureS_heatmap_k_distance_AF.pdf"), p_heat,
       width = 5, height = 3.5, device = cairo_pdf)

cat("Figure 4 saved to:", FIG_DIR, "\n")
