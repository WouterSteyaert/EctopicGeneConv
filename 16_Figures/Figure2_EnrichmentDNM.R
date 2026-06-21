#!/usr/bin/env Rscript
# ==============================================================================
# Figure 2: Gene conversion enrichment in gnomAD variants and de novo mutations
# ==============================================================================
# Panels:
#   (a) Heatmap: Positional enrichment (gnomAD)
#   (b) Heatmap: Concordance excess (gnomAD)
#   (c) Barplot: BLAST-based tract length validation (excess % per AF bin)
#   (d) Line: Enrichment across k for gnomAD vs DNM vs Random
#   (e) Heatmap: DNM concordance rate (ConcorVar / (ConcorVar + DiscorVar))
#   (f) Line: Parent-of-origin enrichment (Father vs Mother, AF < 1e-5)
#
# Reads: stats_allmapp/*_StatsSummary.pos.txt
#        dnm_analysis/export/stats_allmapp/*_StatsSummary.pos.txt
#        tractlengthsSumWgs/_TractLengthAnalysis.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  if (!requireNamespace("cowplot", quietly = TRUE))
    install.packages("cowplot", repos = "https://cloud.r-project.org")
  library(cowplot)
})

# ---- Helper: simple melt (base R, no reshape2 needed) ----
simple_melt <- function(df, id.var = "k") {
  val_cols <- setdiff(names(df), id.var)
  do.call(rbind, lapply(val_cols, function(col) {
    data.frame(k = df[[id.var]], maf = col, val = df[[col]],
               stringsAsFactors = FALSE)
  }))
}

# ---- Paths ----
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
DATA_DIR <- file.path(project_root, "geneconv_complete")
FIG_DIR  <- file.path(project_root, "geneconv_complete/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Constants ----
k_values <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)

maf_codes <- c("0_1e-05", "1e-05_00001", "00001_0001", "0001_0005",
               "0005_001", "001_005", "005_01", "01_05", "05_2")

maf_labels <- c("<1e-5", "1e-5\u20131e-4", "1e-4\u20131e-3",
                "1e-3\u20135e-3", "5e-3\u20130.01",
                "0.01\u20130.05", "0.05\u20130.1", "0.1\u20130.5", ">0.5")

# ---- Theme ----
theme_fig <- theme_classic(base_size = 9) +
  theme(
    text             = element_text(family = "sans"),
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(size = 9),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    plot.margin      = margin(5, 10, 5, 5)
  )

# ---- Helper: read StatsSummary.pos.txt ----
# Format: 3 tab-separated blocks separated by blank lines
#   Block 1: PercDiff_Pos  (relative positional enrichment %)
#   Block 2: ExcessPercPos (absolute excess as fraction of all variants %)
#   Block 3: ExcessPercVar (concordance excess %)
read_pos_txt <- function(filepath) {
  lines <- readLines(filepath)
  header <- strsplit(lines[1], "\t")[[1]]
  header <- header[header != ""]

  data_lines <- lines[-1]
  # Trim trailing empties
  while (length(data_lines) > 0 && trimws(data_lines[length(data_lines)]) == "")
    data_lines <- data_lines[-length(data_lines)]

  blank_idx <- which(trimws(data_lines) == "")

  parse_block <- function(ll) {
    mat <- do.call(rbind, lapply(ll, function(l) as.numeric(strsplit(l, "\t")[[1]])))
    df <- as.data.frame(mat)
    colnames(df) <- c("k", header)
    df
  }

  blocks <- list()
  if (length(blank_idx) == 0) {
    blocks[[1]] <- parse_block(data_lines)
  } else if (length(blank_idx) == 1) {
    blocks[[1]] <- parse_block(data_lines[1:(blank_idx[1] - 1)])
    rem <- data_lines[(blank_idx[1] + 1):length(data_lines)]
    rem <- rem[trimws(rem) != ""]
    if (length(rem) > 0) blocks[[2]] <- parse_block(rem)
  } else {
    blocks[[1]] <- parse_block(data_lines[1:(blank_idx[1] - 1)])
    blocks[[2]] <- parse_block(data_lines[(blank_idx[1] + 1):(blank_idx[2] - 1)])
    rem <- data_lines[(blank_idx[2] + 1):length(data_lines)]
    rem <- rem[trimws(rem) != ""]
    if (length(rem) > 0) blocks[[3]] <- parse_block(rem)
  }
  blocks
}

# ---- Helper: format cell values ----
fmt <- function(x) {
  ifelse(is.na(x), "",
    ifelse(abs(x) >= 100, sprintf("%.0f", x),
      ifelse(abs(x) >= 10, sprintf("%.0f", x),
        ifelse(abs(x) >= 1, sprintf("%.1f", x),
          sprintf("%.2f", x)))))
}

# ---- Read data ----
gnomad   <- read_pos_txt(file.path(DATA_DIR, "stats_allmapp",
              "_All_1000000_1000000_gnomad~genome_StatsSummary.pos.txt"))
dnm_all  <- read_pos_txt(file.path(DATA_DIR, "dnm_analysis/export/stats_allmapp",
              "_All_1000000_1000000_ALL_StatsSummary.pos.txt"))
random3 <- read_pos_txt(file.path(DATA_DIR, "stats_allmapp",
             "_All_1000000_1000000_random3_StatsSummary.pos.txt"))

# DNM concordance rates from R.out.txt
dnm_rout <- read.delim(file.path(DATA_DIR, "dnm_analysis/export/stats_allmapp",
              "_All_1000000_1000000_ALL_StatsSummary.R.out.txt"),
              header = TRUE, stringsAsFactors = FALSE)

# Tract length analysis (BLAST-based validation)
tract <- read.delim(file.path(DATA_DIR, "tractlengthsSumWgs",
           "_TractLengthAnalysis.txt"),
           header = TRUE, stringsAsFactors = FALSE)

# Parent-of-origin DNM data — use export2 (chrX/Y excluded to prevent
# hemizygosity bias in parent-of-origin analysis)
dnm_father <- read_pos_txt(file.path(DATA_DIR, "dnm_analysis/export2/stats_allmapp",
               "_All_1000000_1000000_FATHER_ORIGIN_StatsSummary.pos.txt"))
dnm_mother <- read_pos_txt(file.path(DATA_DIR, "dnm_analysis/export2/stats_allmapp",
               "_All_1000000_1000000_MOTHER_ORIGIN_StatsSummary.pos.txt"))

cat("All data loaded.\n")

# ===========================================================================
# Panel (a): gnomAD positional enrichment heatmap
# ===========================================================================
enrich <- gnomad[[1]]  # PercDiff_Pos

df_a <- simple_melt(enrich)
df_a$maf <- factor(df_a$maf, levels = maf_codes, labels = maf_labels)
df_a$k   <- factor(df_a$k,   levels = rev(k_values))
df_a$lab <- fmt(df_a$val)
# Text white on dark tiles (sqrt-scaled: val > 100 is dark blue)
df_a$tcol <- ifelse(df_a$val > 100, "white", "grey20")

p_a <- ggplot(df_a, aes(maf, k, fill = val)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = lab, colour = tcol), size = 2.1, show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1",
                "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"),
    trans  = "sqrt",
    name   = "Enrichment (%)",
    breaks = c(0, 5, 25, 100, 250, 450),
    limits = c(0, NA)
  ) +
  labs(x = "Allele frequency", y = "k-mer length") +
  ggtitle("a") +
  theme_fig +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        panel.grid  = element_blank(),
        legend.key.height = unit(0.9, "cm"),
        legend.key.width  = unit(0.35, "cm"))

# ===========================================================================
# Panel (b): gnomAD concordance excess heatmap
# ===========================================================================
concor <- gnomad[[3]]  # ExcessPercVar

df_b <- simple_melt(concor)
df_b$maf <- factor(df_b$maf, levels = maf_codes, labels = maf_labels)
df_b$k   <- factor(df_b$k,   levels = rev(k_values))
df_b$lab <- fmt(df_b$val)
df_b$tcol <- ifelse(df_b$val > 3.5, "white", "grey20")
df_b$ns_mark <- ifelse(as.character(df_b$k) == "17" & df_b$maf == "<1e-5", "\u00d7", "")

p_b <- ggplot(df_b, aes(maf, k, fill = val)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = lab, colour = tcol), size = 2.1, show.legend = FALSE) +
  geom_text(aes(label = ns_mark), size = 5, colour = "grey40", show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B",
                "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"),
    name   = "Concordance\nexcess (%)",
    breaks = c(0, 1, 2, 4, 6),
    limits = c(0, NA)
  ) +
  labs(x = "Allele frequency", y = "k-mer length") +
  ggtitle("b") +
  theme_fig +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        panel.grid  = element_blank(),
        legend.key.height = unit(0.9, "cm"),
        legend.key.width  = unit(0.35, "cm"))

# ===========================================================================
# Panel (c): BLAST-based tract length validation (excess % per AF bin)
# ===========================================================================
tract_maf_labels <- c("<1e-5", "1e-5\u20131e-4", "1e-4\u20131e-3",
                       "1e-3\u20135e-3", "5e-3\u20130.01",
                       "0.01\u20130.05", "0.05\u20130.1", "0.1\u20130.5", ">0.5")

df_c <- data.frame(
  maf     = factor(tract_maf_labels, levels = tract_maf_labels),
  excess  = tract$RelDifference
)

p_c <- ggplot(df_c, aes(maf, excess)) +
  geom_col(width = 0.7, fill = "#4292C6", colour = "grey30", linewidth = 0.2) +
  labs(x = "Allele frequency", y = "Variant excess (%)") +
  ggtitle("c") +
  theme_fig +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5))

# ===========================================================================
# Panel (d): gnomAD vs DNM vs Random — mean enrichment across k
# ===========================================================================
avg_cols <- maf_codes  # all frequency bins

make_avg <- function(block, label) {
  data.frame(
    k      = block$k,
    excess = rowMeans(block[, avg_cols, drop = FALSE], na.rm = TRUE),
    source = label
  )
}

gnomad_avg <- make_avg(gnomad[[1]], "gnomAD")
dnm_avg    <- make_avg(dnm_all[[1]], "De novo mutations")

# Random control (random3 set — single "all" column)
rand_avg  <- data.frame(k      = random3[[1]]$k,
                         excess = random3[[1]][, 2],
                         source = "Random control")

df_d <- rbind(gnomad_avg, dnm_avg, rand_avg)
df_d$source <- factor(df_d$source,
  levels = c("De novo mutations", "gnomAD", "Random control"))
df_d$k <- factor(df_d$k, levels = k_values)

p_d <- ggplot(df_d, aes(k, excess, colour = source, shape = source, group = source)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("#E41A1C", "#377EB8", "#999999")) +
  scale_shape_manual(values  = c(16, 17, 15)) +
  labs(x = "k-mer length", y = "Mean positional enrichment (%)",
       colour = NULL, shape = NULL) +
  ggtitle("d") +
  theme_fig +
  theme(legend.position    = c(0.35, 0.88),
        legend.background  = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size    = unit(0.4, "cm"),
        legend.spacing.y   = unit(0.1, "cm"))

# ===========================================================================
# Panel (e): DNM concordance rate heatmap
# ===========================================================================
# Concordance rate = ConcorVar / (ConcorVar + DiscorVar) × 100
dnm_rout$concor_rate <- dnm_rout$ConcorVar /
  (dnm_rout$ConcorVar + dnm_rout$DiscorVar) * 100

# Map FrequencyInterval to MAF codes
dnm_rout$maf <- dnm_rout$FrequencyInterval
dnm_rout$k   <- dnm_rout$RepLength

# Build matrix
df_e <- dnm_rout[, c("k", "maf", "concor_rate")]
df_e$maf <- factor(df_e$maf, levels = maf_codes, labels = maf_labels)
df_e$k   <- factor(df_e$k,   levels = rev(k_values))
df_e$lab <- fmt(df_e$concor_rate)
df_e$tcol <- ifelse(df_e$concor_rate > 90, "white", "grey20")

p_e <- ggplot(df_e, aes(maf, k, fill = concor_rate)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = lab, colour = tcol), size = 2.1, show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B",
                "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),
    name   = "Concordance\nrate (%)",
    breaks = c(50, 60, 70, 80, 90, 100),
    limits = c(45, 100)
  ) +
  labs(x = "Allele frequency", y = "k-mer length") +
  ggtitle("e") +
  theme_fig +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        panel.grid  = element_blank(),
        legend.key.height = unit(0.9, "cm"),
        legend.key.width  = unit(0.35, "cm"))

# ===========================================================================
# Panel (f): Parent-of-origin enrichment (Father vs Mother, AF < 1e-5)
# ===========================================================================
df_f <- rbind(
  data.frame(k = dnm_father[[1]]$k,
             enrichment = dnm_father[[1]][, "0_1e-05"],
             origin = "Paternal"),
  data.frame(k = dnm_mother[[1]]$k,
             enrichment = dnm_mother[[1]][, "0_1e-05"],
             origin = "Maternal")
)
df_f$origin <- factor(df_f$origin, levels = c("Maternal", "Paternal"))
df_f$k <- factor(df_f$k, levels = k_values)

p_f <- ggplot(df_f, aes(k, enrichment, colour = origin, shape = origin, group = origin)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("Maternal" = "#D62728", "Paternal" = "#1F77B4")) +
  scale_shape_manual(values  = c("Maternal" = 16, "Paternal" = 17)) +
  labs(x = "k-mer length", y = "Positional enrichment (%)",
       colour = NULL, shape = NULL) +
  ggtitle("f") +
  theme_fig +
  theme(legend.position    = c(0.35, 0.88),
        legend.background  = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size    = unit(0.4, "cm"),
        legend.spacing.y   = unit(0.1, "cm"))

# ===========================================================================
# Combine: 3 rows × 2 columns
# ===========================================================================
row1 <- plot_grid(p_a, p_b, ncol = 2, align = "hv", labels = NULL)
row2 <- plot_grid(p_c, p_d, ncol = 2, align = "hv",
                  rel_widths = c(1, 1.15), labels = NULL)
row3 <- plot_grid(p_e, p_f, ncol = 2, align = "hv",
                  rel_widths = c(1, 1.15), labels = NULL)

fig2 <- plot_grid(row1, row2, row3,
                  ncol = 1,
                  rel_heights = c(1.15, 0.85, 1.15))

ggsave(file.path(FIG_DIR, "Figure2_EnrichmentDNM.pdf"), fig2,
       width = 7.5, height = 10, device = cairo_pdf)

# Individual panels
ggsave(file.path(FIG_DIR, "Figure2a_gnomAD_enrichment.pdf"), p_a,
       width = 4.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure2b_concordance.pdf"), p_b,
       width = 4.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure2c_tractlength.pdf"), p_c,
       width = 4, height = 3.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure2d_gnomAD_vs_DNM.pdf"), p_d,
       width = 4, height = 3.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure2e_DNM_concordance.pdf"), p_e,
       width = 4.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure2f_parent_of_origin.pdf"), p_f,
       width = 4, height = 3.5, device = cairo_pdf)

# Source data TSVs
TSV_DIR <- file.path(FIG_DIR, "source_data")
dir.create(TSV_DIR, showWarnings = FALSE, recursive = TRUE)

write.table(df_a[, c("k", "maf", "val")], file.path(TSV_DIR, "Figure2a_gnomAD_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_b[, c("k", "maf", "val")], file.path(TSV_DIR, "Figure2b_concordance_excess.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_c, file.path(TSV_DIR, "Figure2c_tractlength.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_d, file.path(TSV_DIR, "Figure2d_gnomAD_vs_DNM.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_e[, c("k", "maf", "concor_rate")], file.path(TSV_DIR, "Figure2e_DNM_concordance.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_f, file.path(TSV_DIR, "Figure2f_parent_of_origin.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Figure 2 saved to:", FIG_DIR, "\n")
cat("Source data TSVs saved to:", TSV_DIR, "\n")
