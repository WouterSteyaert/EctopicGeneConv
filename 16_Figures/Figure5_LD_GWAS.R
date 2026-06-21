#!/usr/bin/env Rscript
# ==============================================================================
# Figure 5 — NEW 2-panel version (clean LD + GWAS depletion)
# Panel a: paired bar plot mean r^2 GC vs Non-GC (k = 21, 41, 71) — clean data
# Panel b: GWAS-depletion heatmap (k x MAF), all mappable, block 1 PercDiff_Pos
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(cowplot)
})

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
BASE   <- file.path(project_root, "geneconv_complete")
OUTDIR <- file.path(BASE, "LD_analysis_resampling_1to1/figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

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

# =============================================================================
# Panel a: clean LD comparison (mean r^2 GC vs Non-GC)
# =============================================================================
ld <- read.csv(file.path(BASE, "LD_analysis_resampling_1to1",
                          "ld_analysis_results_1to1.csv"),
               header = TRUE, stringsAsFactors = FALSE)

key_k <- c(21, 41, 71)
ld_a <- ld %>%
  filter(mappability == "allmapp", population == "EUR",
         window == "win10kb", k %in% key_k)

ld_a$maf_bin <- factor(ld_a$maf_bin,
  levels = c("0_1e-05", "1e-05_0.0001", "0.0001_0.001", "0.001_0.005",
             "0.005_0.01", "0.01_0.05", "0.05_0.1", "0.1_0.5", "0.5_2"),
  labels = c("<1e-5", "1e-5–1e-4", "1e-4–1e-3",
             "1e-3–5e-3", "5e-3–0.01",
             "0.01–0.05", "0.05–0.1", "0.1–0.5", ">0.5"))
ld_a$k <- factor(ld_a$k, levels = key_k)

ld_a_long <- bind_rows(
  ld_a %>% mutate(type = "GC",     mean_r2 = gc_mean),
  ld_a %>% mutate(type = "Non-GC", mean_r2 = nongc_mean)
)
ld_a_long$type <- factor(ld_a_long$type, levels = c("Non-GC", "GC"))

p_a <- ggplot(ld_a_long, aes(x = maf_bin, y = mean_r2, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_wrap(~ paste0("k = ", k), ncol = 1) +
  scale_fill_manual(values = c("Non-GC" = "#377EB8", "GC" = "#E41A1C"),
                     name = NULL) +
  labs(x = "Allele frequency",
       y = expression(paste("Mean ", r^2))) +
  ggtitle("a") +
  theme_fig +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = c(0.12, 0.97),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size = unit(0.35, "cm"),
        strip.text = element_text(face = "bold", size = 9))

# =============================================================================
# Panel b: GWAS-depletion heatmap (block 1 PercDiff_Pos, all mappable)
# =============================================================================
read_block1 <- function(fp) {
  lines <- readLines(fp)
  blank <- which(trimws(lines) == "")[1]
  if (is.na(blank)) blank <- length(lines) + 1
  txt <- lines[1:(blank - 1)]
  header_raw <- strsplit(txt[1], "\t", fixed = TRUE)[[1]]
  header_raw <- header_raw[nchar(header_raw) > 0]
  mat <- do.call(rbind, lapply(txt[-1], function(l) {
    as.numeric(strsplit(l, "\t", fixed = TRUE)[[1]])
  }))
  df <- as.data.frame(mat)
  colnames(df) <- c("k", header_raw)
  df
}

gwas_fp <- file.path(BASE, "stats_allmapp",
                      "_All_1000000_1000000_gwas~clean~All_StatsSummary.pos.txt")
gwas <- read_block1(gwas_fp)

# Read R.out.txt for log10(p) values (Bonferroni significance flag)
gwas_rout_fp <- file.path(BASE, "stats_allmapp",
                          "_All_1000000_1000000_gwas~clean~All_StatsSummary.R.out.txt")
gwas_rout <- read.delim(gwas_rout_fp, stringsAsFactors = FALSE)
N_TESTS_GWAS <- 90  # 10 k x 9 AF bins
SIG_THRESHOLD_LOG10P <- log10(0.05 / N_TESTS_GWAS)  # = -3.255; Bonferroni p < 0.05

maf_codes_gwas <- c("0_1e-05", "1e-05_00001", "00001_0001",
                     "0001_0005", "0005_001", "001_005",
                     "005_01", "01_05", "05_2")
maf_labels_pretty <- c("<1e-5", "1e-5–1e-4", "1e-4–1e-3",
                       "1e-3–5e-3", "5e-3–0.01",
                       "0.01–0.05", "0.05–0.1",
                       "0.1–0.5", ">0.5")

gwas_long <- gwas %>%
  tidyr::pivot_longer(-k, names_to = "maf_bin", values_to = "depletion")
gwas_long$maf_bin_code <- gwas_long$maf_bin
gwas_long$maf_bin <- factor(gwas_long$maf_bin,
                             levels = maf_codes_gwas,
                             labels = maf_labels_pretty)
# Attach Bonferroni-corrected significance flag
gwas_long$log10p <- gwas_rout$log_p_value_exponent[
  match(paste(gwas_long$k, gwas_long$maf_bin_code),
        paste(gwas_rout$RepLength, gwas_rout$FrequencyInterval))]
gwas_long$log10p_bonf <- gwas_long$log10p + log10(N_TESTS_GWAS)
gwas_long$sig <- gwas_long$log10p_bonf < log10(0.05)
gwas_long$label <- ifelse(gwas_long$sig,
                          sprintf("%.0f", gwas_long$depletion),
                          sprintf("%.0f ns", gwas_long$depletion))
gwas_long$k <- factor(gwas_long$k, levels = sort(unique(gwas_long$k)))

p_b <- ggplot(gwas_long, aes(x = maf_bin, y = k, fill = depletion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 2.4,
            color = "black") +
  scale_fill_gradient2(low = "#E41A1C", mid = "white", high = "#4DAF4A",
                        midpoint = 0,
                        limits = c(-60, 60), oob = scales::squish,
                        name = "% deviation") +
  labs(x = "Allele frequency", y = "k") +
  ggtitle("b") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(face = "bold", size = 11, hjust = 0))

# =============================================================================
# Compose 2-panel figure
# =============================================================================
fig5 <- plot_grid(p_a, p_b, ncol = 2, rel_widths = c(1, 1.2),
                   align = "h", axis = "tb")

ggsave(file.path(OUTDIR, "Figure5_new.pdf"), fig5,
       width = 10, height = 7, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Figure5_new.png"), fig5,
       width = 10, height = 7, dpi = 150)

ggsave(file.path(OUTDIR, "Figure5a_clean.pdf"), p_a,
       width = 4.8, height = 6.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Figure5b_GWAS_depletion.pdf"), p_b,
       width = 5.5, height = 5, device = cairo_pdf)

# Source data export for paper deposit (Figure 5b)
write.table(
  data.frame(k = as.integer(as.character(gwas_long$k)),
             maf_bin = gwas_long$maf_bin_code,
             pct_depletion = round(gwas_long$depletion, 2),
             log10p_bonferroni = round(gwas_long$log10p_bonf, 3),
             significant_bonferroni = gwas_long$sig),
  file.path(OUTDIR, "Figure5b_GWAS_depletion.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

cat("Wrote:\n")
cat("  ", file.path(OUTDIR, "Figure5_new.pdf"), "\n")
cat("  ", file.path(OUTDIR, "Figure5_new.png"), "\n")
