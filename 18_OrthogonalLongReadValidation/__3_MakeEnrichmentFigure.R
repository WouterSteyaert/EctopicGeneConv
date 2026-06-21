#!/usr/bin/env Rscript
# Variant-density enrichment comparison: gnomAD vs GA4K HiFi calls per
# (mappability stratum, k, AF bin), using identical positions, denominators
# and formula as the gnomAD contingency analysis. Plus the LRS/gnomAD ratio
# panel.
#
# Input  : <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/v2_enrichment_compare.tsv
# Output : <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/figures/
#            Enrichment_gnomAD_vs_LRS.{pdf,png}
#            EnrichmentRatio_LRS_over_gnomAD.{pdf,png}
#
# Environment: PROJECT_ROOT must point to the project root.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (PROJECT_ROOT == "" || !dir.exists(PROJECT_ROOT))
  stop("PROJECT_ROOT env var must point to a directory")

BASE   <- file.path(PROJECT_ROOT, "geneconv_complete", "artefact_validation_ga4k")
RES    <- file.path(BASE, "results")
FIGDIR <- file.path(BASE, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

d <- read.delim(file.path(RES, "v2_enrichment_compare.tsv"), stringsAsFactors = FALSE)

# Restrict to AF bins both technologies can see (>= 0.005)
d <- d |> filter(af_bin %in% c("0.005_0.01","0.01_0.05","0.05_0.1","0.1_0.5","0.5_2"))

strat_lvls <- c("allmapp", "nosegdupmapp", "segdupmapp")
strat_lbl  <- c("All mappable", "Outside SDs", "Within SDs")
d$stratum  <- factor(d$stratum, levels = strat_lvls, labels = strat_lbl)
d$k        <- factor(d$k, levels = c(17,19,21,31,41,51,61,71,81,91))
af_lvls    <- c("0.005_0.01","0.01_0.05","0.05_0.1","0.1_0.5","0.5_2")
af_lbl     <- c("0.005-0.01","0.01-0.05","0.05-0.1","0.1-0.5","0.5-1")
d$af_bin   <- factor(d$af_bin, levels = af_lvls, labels = af_lbl)

d_long <- d |>
  select(stratum, k, af_bin, gnomAD = enr_gnomad, `GA4K (HiFi)` = enr_lrs) |>
  pivot_longer(c(gnomAD, `GA4K (HiFi)`), names_to = "source", values_to = "enrichment")
d_long$source <- factor(d_long$source, levels = c("gnomAD", "GA4K (HiFi)"))

p1 <- ggplot(d_long, aes(x = k, y = enrichment, color = source, group = source)) +
  geom_hline(yintercept = 1, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  facet_grid(af_bin ~ stratum, scales = "free_y") +
  scale_color_manual(values = c("gnomAD" = "#377EB8", "GA4K (HiFi)" = "#E41A1C"),
                     name = NULL) +
  scale_y_log10() +
  labs(x = "Template length k",
       y = "Variant-density enrichment (GC / non-GC), log scale") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIGDIR, "Enrichment_gnomAD_vs_LRS.pdf"), p1,
       width = 12, height = 11, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Enrichment_gnomAD_vs_LRS.png"), p1,
       width = 12, height = 11, dpi = 150)

p2 <- ggplot(d, aes(x = k, y = ratio_lrs_over_gnomad, color = af_bin, group = af_bin)) +
  geom_hline(yintercept = 1, color = "black", linewidth = 0.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  facet_wrap(~ stratum, ncol = 3) +
  scale_color_brewer(palette = "Spectral", direction = -1, name = "AF bin") +
  labs(x = "Template length k",
       y = "LRS enrichment / gnomAD enrichment") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIGDIR, "EnrichmentRatio_LRS_over_gnomAD.pdf"), p2,
       width = 11, height = 4.5, device = cairo_pdf)
ggsave(file.path(FIGDIR, "EnrichmentRatio_LRS_over_gnomAD.png"), p2,
       width = 11, height = 4.5, dpi = 150)

cat("Enrichment figures written\n")
