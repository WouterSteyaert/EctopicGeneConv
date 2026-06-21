#!/usr/bin/env Rscript
# =============================================================================
# Per-generation gene conversion rate (r_gc) and confidence intervals.
#
# Estimates the per-site per-generation rate of ectopic gene conversion-driven
# allele change from the DNM contingency tables (step 06 + step 07), for each
# combination of template length k and three mappability strata
# (all mappable, segmental duplications, non-SD mappable).
#
# Reads:
#   <denovo_work_dir>/export/<stratum>/_All_<W>_<S>_ALL_StatsSummary.R.out.txt
#
# Configuration is read from config.GRCh38.ini via the --ConfigFile= argument
# (see 00_Configuration/README.md).  Set DNM_EXPORT_DIR via ProjectConfig:
#   paths.denovo_work_dir       (top-level DNM work directory)
#   enrichment.window_size      (Mb worker output window)
#   enrichment.step_size        (Mb worker step)
# Strata names follow the step-06 convention: stats_<allmapp|segdupmapp|nosegdupmapp>.
#
# Outputs:
#   <DNM_EXPORT_DIR>/per_generation_gc_rate.tsv     full table with CIs
#   <DNM_EXPORT_DIR>/per_generation_gc_rate.pdf     two-panel figure
#
# Usage:
#   Rscript PerGenerationConversionRate.R \
#       --config=../00_Configuration/config.GRCh38.ini \
#       --project-root=$PROJECT_ROOT
# =============================================================================

# ── Argument parsing (lightweight; no external packages) ────────────────────
args         <- commandArgs(trailingOnly = TRUE)
config_file  <- NA
project_root <- Sys.getenv("PROJECT_ROOT")
stats_root   <- NA   # optional override: where to find stats_<context>/ dirs
for (a in args) {
    if (grepl("^--config=", a))       config_file  <- sub("^--config=",       "", a)
    if (grepl("^--project-root=", a)) project_root <- sub("^--project-root=", "", a)
    if (grepl("^--stats-root=", a))   stats_root   <- sub("^--stats-root=",   "", a)
}
stopifnot(!is.na(config_file), nzchar(project_root))

# ── Minimal config parser (Config::General sections we need) ────────────────
read_config <- function(path, project_root) {
    lines <- readLines(path, warn = FALSE)
    section <- ""
    cfg <- list()
    for (ln in lines) {
        ln <- sub("#.*$", "", ln)
        ln <- trimws(ln)
        if (!nchar(ln)) next
        if (grepl("^<(/?)([A-Za-z_]+)>$", ln)) {
            m <- regmatches(ln, regexec("^<(/?)([A-Za-z_]+)>$", ln))[[1]]
            section <- if (m[2] == "/") "" else m[3]
            if (!nchar(m[2])) cfg[[section]] <- list()
            next
        }
        if (grepl("=", ln) && nchar(section)) {
            kv <- strsplit(ln, "=", fixed = FALSE)[[1]]
            key <- trimws(kv[1])
            val <- trimws(paste(kv[-1], collapse = "="))
            val <- gsub("__PROJECT_ROOT__", project_root, val, fixed = TRUE)
            cfg[[section]][[key]] <- val
        }
    }
    cfg
}
cfg <- read_config(config_file, project_root)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Parameters ──────────────────────────────────────────────────────────────
N_TRIOS    <- 11963
MU_DNM     <- 1.2e-8           # Jónsson et al. 2017 (per bp per generation)
FREQ_BIN   <- "0_1e-05"        # truly de novo (gnomAD AF < 1e-5)
W          <- as.integer(cfg$enrichment$window_size %||% 1000000)
S          <- as.integer(cfg$enrichment$step_size   %||% 1000000)

# Where to find the DNM-specific stats_<context>/ summary tables:
#   1. --stats-root= override (highest priority)
#   2. <denovo_work_dir>/export   (the convention used in the working pipeline)
#   3. <stats_dir>                (matches the step-06 default output location)
BASE_DIR <- if (!is.na(stats_root)) {
    stats_root
} else if (dir.exists(file.path(cfg$paths$denovo_work_dir, "export"))) {
    file.path(cfg$paths$denovo_work_dir, "export")
} else {
    cfg$paths$stats_dir
}
OUTPUT_DIR <- BASE_DIR

STRATA <- data.frame(
  dir   = c("stats_allmapp", "stats_segdupmapp", "stats_nosegdupmapp"),
  label = c("All mappable", "Segmental duplications", "Non-segdup mappable"),
  col   = c("steelblue", "firebrick", "grey40"),
  stringsAsFactors = FALSE
)

# ── Read and compute for each stratum ───────────────────────────────────────
all_results <- list()

for (i in seq_len(nrow(STRATA))) {
  infile <- file.path(BASE_DIR, STRATA$dir[i],
                      sprintf("_All_%d_%d_ALL_StatsSummary.R.out.txt", W, S))
  if (!file.exists(infile)) {
    warning(sprintf("Missing summary file for stratum '%s': %s", STRATA$label[i], infile))
    next
  }
  dt <- read.delim(infile, check.names = TRUE)
  dt <- dt[dt$FrequencyInterval == FREQ_BIN, ]
  dt <- dt[order(as.integer(dt$RepLength)), ]
  dt$k <- as.integer(dt$RepLength)

  dt$ConcorVar   <- as.numeric(dt$ConcorVar)
  dt$ConcorNoVar <- as.numeric(dt$ConcorNoVar)
  dt$DiscorVar   <- as.numeric(dt$DiscorVar)
  dt$DiscorNoVar <- as.numeric(dt$DiscorNoVar)
  dt$ConcorTotal <- dt$ConcorVar + dt$ConcorNoVar
  dt$DiscorTotal <- dt$DiscorVar + dt$DiscorNoVar
  dt$ConcorFreq  <- dt$ConcorVar / dt$ConcorTotal
  dt$DiscorFreq  <- dt$DiscorVar / dt$DiscorTotal

  dt$ExcessConcordant <- dt$ConcorVar - dt$DiscorFreq * dt$ConcorTotal
  dt$fold_excess      <- dt$ConcorFreq / dt$DiscorFreq
  dt$r_gc_direct      <- dt$ExcessConcordant / (N_TRIOS * dt$ConcorTotal)
  dt$r_gc_lit         <- (MU_DNM / 3) * ((dt$ConcorFreq - dt$DiscorFreq) / dt$DiscorFreq)
  dt$mu_effective     <- dt$DiscorFreq / N_TRIOS

  # ── 95% CIs ─────────────────────────────────────────────────────────────
  # Fold excess via delta method on the log scale:
  #   log(fold) ± 1.96 * sqrt(1/ConcorVar + 1/DiscorVar)
  se_log_fold   <- sqrt(1 / dt$ConcorVar + 1 / dt$DiscorVar)
  dt$fold_lo    <- exp(log(dt$fold_excess) - 1.96 * se_log_fold)
  dt$fold_hi    <- exp(log(dt$fold_excess) + 1.96 * se_log_fold)

  # r_gc / mu_eff = fold - 1 (linear shift of fold CI)
  dt$rgc_over_mu_eff    <- dt$fold_excess - 1
  dt$rgc_over_mu_eff_lo <- dt$fold_lo     - 1
  dt$rgc_over_mu_eff_hi <- dt$fold_hi     - 1
  dt$rgc_over_mu_eff_lo[dt$rgc_over_mu_eff_lo < 0] <- 0   # clip negative lower bounds

  # r_gc itself: Poisson variance on Excess
  #   Var(Excess) = ConcorVar + (ConcorTotal/DiscorTotal)^2 * DiscorVar
  var_excess     <- dt$ConcorVar + (dt$ConcorTotal / dt$DiscorTotal)^2 * dt$DiscorVar
  se_excess      <- sqrt(var_excess)
  dt$r_gc_lo_raw <- (dt$ExcessConcordant - 1.96 * se_excess) / (N_TRIOS * dt$ConcorTotal)
  dt$r_gc_hi     <- (dt$ExcessConcordant + 1.96 * se_excess) / (N_TRIOS * dt$ConcorTotal)
  dt$r_gc_lo     <- pmax(dt$r_gc_lo_raw, 0)   # clip negative lower bounds (see Methods)

  # Diagnostic columns (referenced in Methods rationale for clipping at k=17 SDs):
  #   - exact Poisson 2.5% lower bound on ConcorVar (qchisq method)
  #   - the background expected count = DiscorFreq * ConcorTotal
  #   - sqrt(ConcorVar) (count-noise scale)
  dt$ConcorVar_PoisLo  <- qchisq(0.025, 2 * dt$ConcorVar) / 2
  dt$BackgroundExpect  <- dt$DiscorFreq * dt$ConcorTotal
  dt$SqrtConcorVar     <- sqrt(dt$ConcorVar)

  dt$stratum <- STRATA$label[i]
  all_results[[i]] <- dt
}

combined <- do.call(rbind, all_results)

# ── Output table ────────────────────────────────────────────────────────────
out <- data.frame(
  stratum             = combined$stratum,
  k                   = combined$k,
  ConcorVar           = combined$ConcorVar,
  ConcorTotal         = combined$ConcorTotal,
  DiscorVar           = combined$DiscorVar,
  DiscorTotal         = combined$DiscorTotal,
  ExcessConcordant    = round(combined$ExcessConcordant, 1),
  ExcessPerTrio       = round(combined$ExcessConcordant / N_TRIOS, 3),
  fold_excess         = round(combined$fold_excess, 2),
  fold_lo             = round(combined$fold_lo, 2),
  fold_hi             = round(combined$fold_hi, 2),
  r_gc_direct         = sprintf("%.3e", combined$r_gc_direct),
  r_gc_lo             = sprintf("%.3e", combined$r_gc_lo),
  r_gc_lo_raw         = sprintf("%.3e", combined$r_gc_lo_raw),
  r_gc_hi             = sprintf("%.3e", combined$r_gc_hi),
  r_gc_lit            = sprintf("%.3e", combined$r_gc_lit),
  ConcorVar_PoisLo    = round(combined$ConcorVar_PoisLo, 1),
  BackgroundExpect    = round(combined$BackgroundExpect, 1),
  SqrtConcorVar       = round(combined$SqrtConcorVar, 1),
  rgc_over_mu_eff     = round(combined$rgc_over_mu_eff,    2),
  rgc_over_mu_eff_lo  = round(combined$rgc_over_mu_eff_lo, 2),
  rgc_over_mu_eff_hi  = round(combined$rgc_over_mu_eff_hi, 2),
  mu_eff_per_allele   = sprintf("%.3e", combined$mu_effective)
)

cat("\n=== Per-generation GC-driven allele change rate (gnomAD AF < 1e-5) ===\n")
cat(sprintf("N trios = %d | mu_dnm (lit.) = %.1e\n\n", N_TRIOS, MU_DNM))
for (s in STRATA$label) {
  cat(sprintf("── %s ──\n", s))
  sub <- out[out$stratum == s, -1]
  print(sub, row.names = FALSE)
  cat("\n")
}

out_path <- file.path(OUTPUT_DIR, "per_generation_gc_rate.tsv")
write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Table written to: %s\n", out_path))

# ── Figure ──────────────────────────────────────────────────────────────────
fig_path <- file.path(OUTPUT_DIR, "per_generation_gc_rate.pdf")
pdf(fig_path, width = 8, height = 8.5)
par(mfrow = c(2, 1), mar = c(4.5, 5.5, 3, 1), mgp = c(3.5, 0.8, 0))

k_vals <- sort(unique(combined$k))

# Panel A: fold excess + 95% CI bars
ylim_fold <- range(c(combined$fold_lo, combined$fold_hi, 1), na.rm = TRUE)
plot(NA, xlim = range(k_vals), ylim = ylim_fold,
     xlab = "Minimum repeat length k (bp)",
     ylab = "Fold over background\n(concordant / discordant rate)",
     main = "A", xaxt = "n", las = 1)
axis(1, at = k_vals)
abline(h = 1, lty = 2, col = "grey70")
for (i in seq_len(nrow(STRATA))) {
  sub <- combined[combined$stratum == STRATA$label[i], ]
  arrows(sub$k, sub$fold_lo, sub$k, sub$fold_hi,
         length = 0.04, angle = 90, code = 3, col = STRATA$col[i])
  lines(sub$k, sub$fold_excess, type = "b", pch = 19, cex = 1.3, lwd = 2,
        col = STRATA$col[i])
}
legend("topleft", legend = STRATA$label, col = STRATA$col,
       pch = 19, lwd = 2, bty = "n", cex = 0.85)
mtext(sprintf("N = %s trios | gnomAD AF < 1e-5", format(N_TRIOS, big.mark = ",")),
      side = 3, line = 0.2, cex = 0.8, col = "grey40")

# Panel B: r_gc per site per generation (log scale)
pos_mask  <- combined$r_gc_direct > 0
ylim_rate <- range(combined$r_gc_direct[pos_mask], na.rm = TRUE)
plot(NA, xlim = range(k_vals), ylim = ylim_rate,
     xlab = "Minimum repeat length k (bp)",
     ylab = expression("r"["gc"] ~ "(per site per generation)"),
     main = "B", xaxt = "n", las = 1, log = "y")
axis(1, at = k_vals)
for (i in seq_len(nrow(STRATA))) {
  sub <- combined[combined$stratum == STRATA$label[i], ]
  pos <- sub$r_gc_direct > 0
  lines(sub$k[pos], sub$r_gc_direct[pos], type = "b", pch = 19, cex = 1.3,
        lwd = 2, col = STRATA$col[i])
}
abline(h = MU_DNM / 3, lty = 3, col = "grey50")
text(max(k_vals), MU_DNM / 3, expression(mu / 3), pos = 1, cex = 0.8, col = "grey50")
legend("bottomright", legend = STRATA$label, col = STRATA$col,
       pch = 19, lwd = 2, bty = "n", cex = 0.85)

dev.off()
cat(sprintf("Figure written to: %s\n", fig_path))
cat("\nDone.\n")
