#!/usr/bin/env Rscript
# AF-adjusted Delta-mismatch forest plot per (mappability stratum, k),
# and worst-case enrichment-inflation bound at the upper 95% CI.
#
# Input  : <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/results/v2_mh_pooled.tsv
# Output : <PROJECT_ROOT>/geneconv_complete/artefact_validation_ga4k/figures/
#            Forest_DeltaMismatch_MH.{pdf,png}
#            MaxEnrichmentInflation.{pdf,png}
#
# Environment: PROJECT_ROOT must point to the project root.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (PROJECT_ROOT == "" || !dir.exists(PROJECT_ROOT))
  stop("PROJECT_ROOT env var must point to a directory")

BASE   <- file.path(PROJECT_ROOT, "geneconv_complete", "artefact_validation_ga4k")
RES    <- file.path(BASE, "results")
FIGDIR <- file.path(BASE, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

d <- read.delim(file.path(RES, "v2_mh_pooled.tsv"), stringsAsFactors = FALSE)

strat_lvls <- c("allmapp", "nosegdupmapp", "segdupmapp")
strat_lbl  <- c("All mappable", "Outside SDs", "Within SDs")
d$stratum  <- factor(d$stratum, levels = strat_lvls, labels = strat_lbl)
d$k        <- factor(d$k, levels = c(17,19,21,31,41,51,61,71,81,91))

p <- ggplot(d, aes(x = k, y = rd_pooled)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.3, linewidth = 0.5, color = "#444444") +
  geom_point(size = 2.5, color = "#E41A1C") +
  facet_wrap(~ stratum, ncol = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(x = "Template length k",
       y = "AF-adjusted Δ mismatch rate (GC − non-GC), 95% CI") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIGDIR, "Forest_DeltaMismatch_MH.pdf"), p,
       width = 11, height = 4.5, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Forest_DeltaMismatch_MH.png"), p,
       width = 11, height = 4.5, dpi = 150)

# Worst-case enrichment inflation (cap visualisation at 1.05x; max actual ~1.01x)
d$max_infl_capped <- pmin(d$max_inflation_at_ci_hi, 1.05)

p2 <- ggplot(d, aes(x = k, y = max_infl_capped)) +
  geom_hline(yintercept = 1, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = c(1.05, 1.10), linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5, color = "#377EB8") +
  geom_segment(aes(x = k, xend = k, y = 1, yend = max_infl_capped),
               color = "#377EB8", linewidth = 0.4) +
  facet_wrap(~ stratum, ncol = 3) +
  scale_y_continuous(breaks = c(1, 1.005, 1.01, 1.02, 1.05),
                     labels = c("1.00x", "1.005x", "1.01x", "1.02x", "1.05x")) +
  labs(x = "Template length k",
       y = "Maximum enrichment-inflation factor (upper 95% CI)") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIGDIR, "MaxEnrichmentInflation.pdf"), p2,
       width = 11, height = 4.5, device = cairo_pdf)
ggsave(file.path(FIGDIR, "MaxEnrichmentInflation.png"), p2,
       width = 11, height = 4.5, dpi = 150)

cat("Forest + max inflation figures written\n")
