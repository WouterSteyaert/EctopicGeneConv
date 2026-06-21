#!/usr/bin/env Rscript
# ==============================================================================
# Supplementary Table: Recombination hotspot enrichment (MAF x k)
#
# For each k x MAF combination:
#   - Pooled log2(OR) inside and outside hotspots
#   - Delta = log2(OR_inside) - log2(OR_outside)
#   - P-value (interaction z-test on log OR ratio)
#   - BH-corrected p-value
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
STATS_DIR <- file.path(project_root, "geneconv_complete/stats_allmapp")
OUT_DIR   <- file.path(project_root, "geneconv_complete/Tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

chrom_sizes <- c(
  "1"=248956422, "2"=242193529, "3"=198295559, "4"=190214555,
  "5"=181538259, "6"=170805979, "7"=159345973, "8"=145138636,
  "9"=138394717, "10"=133797422, "11"=135086622, "12"=133275309,
  "13"=114364328, "14"=107043718, "15"=101991189, "16"=90338345,
  "17"=83257441, "18"=80373285, "19"=58617616, "20"=64444167,
  "21"=46709983, "22"=50818468
)
chroms <- names(chrom_sizes)

maf_bins <- c("0_1e-05", "1e-05_00001", "00001_0001", "0001_0005",
              "0005_001", "001_005", "005_01", "01_05", "05_2")
maf_labels <- c(
  "0_1e-05"     = "0 - 1e-5",
  "1e-05_00001" = "1e-5 - 1e-4",
  "00001_0001"  = "1e-4 - 1e-3",
  "0001_0005"   = "1e-3 - 5e-3",
  "0005_001"    = "5e-3 - 0.01",
  "001_005"     = "0.01 - 0.05",
  "005_01"      = "0.05 - 0.1",
  "01_05"       = "0.1 - 0.5",
  "05_2"        = "0.5 - 2.0"
)

k_values <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)

# ---- Compute all results ----
cat("Computing hotspot enrichment per k x MAF...\n")

results <- list()
for (k in k_values) {
  for (mi in seq_along(maf_bins)) {
    maf <- maf_bins[mi]
    a_in <- 0; b_in <- 0; c_in <- 0; d_in <- 0
    a_out <- 0; b_out <- 0; c_out <- 0; d_out <- 0
    conc_in <- 0; disc_in <- 0; conc_out <- 0; disc_out <- 0

    for (chr in chroms) {
      fg <- sprintf("%s/%d/%s_1000000_1000000_gnomad~genome_%s.1.%d.txt",
                    STATS_DIR, k, chr, maf, chrom_sizes[chr])
      fh <- sprintf("%s/%d/%s_1000000_1000000_gnomad~genome_%s.1.%d.RecombHotspots.txt",
                    STATS_DIR, k, chr, maf, chrom_sizes[chr])
      if (!file.exists(fg) || !file.exists(fh)) next
      g <- read.delim(fg, header = TRUE, stringsAsFactors = FALSE)
      h <- read.delim(fh, header = TRUE, stringsAsFactors = FALSE)

      a_in <- a_in + sum(h$NrOfVarPosConvPos)
      b_in <- b_in + sum(h$NrOfVarPosNoConvPos)
      c_in <- c_in + sum(h$NrOfNoVarPosConvPos)
      d_in <- d_in + sum(h$NrOfNoVarPosNoConvPos)

      a_out <- a_out + sum(g$NrOfVarPosConvPos) - sum(h$NrOfVarPosConvPos)
      b_out <- b_out + sum(g$NrOfVarPosNoConvPos) - sum(h$NrOfVarPosNoConvPos)
      c_out <- c_out + sum(g$NrOfNoVarPosConvPos) - sum(h$NrOfNoVarPosConvPos)
      d_out <- d_out + sum(g$NrOfNoVarPosNoConvPos) - sum(h$NrOfNoVarPosNoConvPos)

      conc_in <- conc_in + sum(h$ConcorVar)
      disc_in <- disc_in + sum(h$DiscorVar)
      conc_out <- conc_out + sum(g$ConcorVar) - sum(h$ConcorVar)
      disc_out <- disc_out + sum(g$DiscorVar) - sum(h$DiscorVar)
    }

    or_in  <- (a_in * d_in) / (b_in * c_in)
    or_out <- (a_out * d_out) / (b_out * c_out)
    log2_or_in  <- log2(or_in)
    log2_or_out <- log2(or_out)
    delta <- log2_or_in - log2_or_out

    # Interaction z-test
    if (all(c(a_in, b_in, c_in, d_in, a_out, b_out, c_out, d_out) > 0)) {
      log_ratio <- log(or_in / or_out)
      se <- sqrt(1/a_in + 1/b_in + 1/c_in + 1/d_in +
                 1/a_out + 1/b_out + 1/c_out + 1/d_out)
      z <- log_ratio / se
      p <- 2 * pnorm(-abs(z))
    } else {
      z <- NA; p <- NA
    }

    # Concordance
    n_conc_in  <- conc_in + disc_in
    n_conc_out <- conc_out + disc_out
    concordance_inside  <- if (n_conc_in > 0)  conc_in / n_conc_in   else NA
    concordance_outside <- if (n_conc_out > 0) conc_out / n_conc_out else NA

    results[[length(results) + 1]] <- data.frame(
      Template_length_k    = k,
      MAF_bin              = maf_labels[maf],
      log2_OR_inside       = round(log2_or_in, 4),
      log2_OR_outside      = round(log2_or_out, 4),
      Delta_log2_OR        = round(delta, 4),
      Direction            = ifelse(delta > 0, "+", "-"),
      Z_score              = round(z, 2),
      P_value              = p,
      Concordance_inside   = round(concordance_inside, 4),
      Concordance_outside  = round(concordance_outside, 4),
      N_variants_inside    = n_conc_in,
      N_variants_outside   = n_conc_out,
      stringsAsFactors     = FALSE
    )
  }
  cat(sprintf("  k = %d done\n", k))
}

tab <- do.call(rbind, results)
rownames(tab) <- NULL

# ---- Bonferroni correction across all k x MAF tests ----
tab$P_adjusted <- p.adjust(tab$P_value, method = "bonferroni")

# ---- Round p-values for display ----
tab$P_value    <- signif(tab$P_value, 3)
tab$P_adjusted <- signif(tab$P_adjusted, 3)

# ---- Sort ----
tab <- tab[order(tab$Template_length_k, match(tab$MAF_bin, maf_labels)), ]

# ---- Write ----
outfile <- file.path(OUT_DIR, "SupTable_RecombHotspots.tsv")
write.table(tab, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nWritten: %s\n", outfile))
cat(sprintf("Rows: %d (%d k x 9 MAF)\n", nrow(tab), length(k_values)))
cat(sprintf("Columns: %s\n", paste(names(tab), collapse = ", ")))

# ---- Summary ----
cat(sprintf("\n=== Summary (Bonferroni, %d tests) ===\n", nrow(tab)))
for (k in k_values) {
  sub <- tab[tab$Template_length_k == k, ]
  n_sig_pos <- sum(sub$P_adjusted < 0.05 & sub$Delta_log2_OR > 0, na.rm = TRUE)
  n_sig_neg <- sum(sub$P_adjusted < 0.05 & sub$Delta_log2_OR < 0, na.rm = TRUE)
  n_ns      <- sum(is.na(sub$P_adjusted) | sub$P_adjusted >= 0.05)
  mean_delta <- mean(sub$Delta_log2_OR, na.rm = TRUE)
  mean_conc  <- mean(sub$Concordance_inside, na.rm = TRUE)
  cat(sprintf("k=%d:  %d/9 sig+ (Bonf<0.05),  %d/9 sig-,  %d/9 NS  | mean delta=%.3f  mean concordance=%.1f%%\n",
              k, n_sig_pos, n_sig_neg, n_ns, mean_delta, 100 * mean_conc))
}
