#!/usr/bin/env Rscript
# Headline figure: GC vs non-GC HiFi confirmation rates per (mappability
# stratum, k), pooled over the reliable AF zone (AF >= 0.01).
#
# Input  : <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/v2_conf_rates_reliableAF.tsv
# Output : <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/figures/Confirmation_GC_vs_nonGC.{pdf,png}
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

d <- read.delim(file.path(RES, "v2_conf_rates_reliableAF.tsv"), stringsAsFactors = FALSE)

strat_lvls <- c("allmapp", "nosegdupmapp", "segdupmapp")
strat_lbl  <- c("All mappable", "Outside SDs", "Within SDs")
d$stratum  <- factor(d$stratum, levels = strat_lvls, labels = strat_lbl)
d$k        <- factor(d$k, levels = c(17,19,21,31,41,51,61,71,81,91))

d_long <- d |>
  select(stratum, k,
         GC_rate = conf_GC, GC_lo, GC_hi,
         nonGC_rate = conf_nonGC, nonGC_lo, nonGC_hi) |>
  pivot_longer(cols = -c(stratum, k),
               names_to = c("class", ".value"),
               names_pattern = "(GC|nonGC)_?(.+)")
d_long$class <- factor(d_long$class, levels = c("nonGC", "GC"),
                       labels = c("non-GC", "GC"))

yr <- range(c(d$GC_lo, d$nonGC_lo, d$GC_hi, d$nonGC_hi))

p <- ggplot(d_long, aes(x = k, y = rate, color = class, group = class)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = class), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  facet_wrap(~ stratum, ncol = 3) +
  scale_color_manual(values = c("non-GC" = "#377EB8", "GC" = "#E41A1C"), name = NULL) +
  scale_fill_manual(values = c("non-GC" = "#377EB8", "GC" = "#E41A1C"), name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     limits = c(max(0.90, yr[1] - 0.005), min(1.0, yr[2] + 0.005))) +
  labs(x = "Template length k",
       y = "HiFi confirmation rate (Wilson 95% CI)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIGDIR, "Confirmation_GC_vs_nonGC.pdf"), p,
       width = 11, height = 4.8, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Confirmation_GC_vs_nonGC.png"), p,
       width = 11, height = 4.8, dpi = 150)
cat("Confirmation figure written\n")
