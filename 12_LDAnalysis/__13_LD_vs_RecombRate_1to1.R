#!/usr/bin/env Rscript
#===============================================================================
# LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 13
# LD vs recombination rate (1:1 paralogs only)
#
# Input:
#   - Per-variant LD files  (LD_analysis_resampling_1to1/{mapp}/{pop}/{win}/)
#   - deCODE recombAvg.bedGraph
#
# Output:
#   - ld_vs_recombrate_1to1.csv
#   - ld_recomb_heatmap_1to1_*.png
#   - ld_recomb_trend_1to1_*.png
#   - ld_recomb_confounding_1to1_*.png
#   - ld_recomb_matrix_1to1_*.tsv
#===============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

#===============================================================================
# Configuration
#===============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
BASE      <- file.path(project_root, "geneconv_complete")
LD_BASE   <- file.path(BASE, "LD_analysis_resampling_1to1")
# Override via RECOMB_BEDGRAPH env var if your bedGraph lives elsewhere.
RECOMB_BG <- Sys.getenv("RECOMB_BEDGRAPH",
                       file.path(BASE, "bed/deCODE/recombAvg.bedGraph"))

MAPP_CATEGORIES <- c("allmapp", "nosegdupmapp", "segdupmapp")
POPULATIONS     <- c("ALL", "EUR")
WINDOW_SIZES    <- c("win10kb", "win25kb", "win50kb", "win100kb")

K_VALUES <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)
MAF_BINS <- c("0.001_0.005", "0.005_0.01", "0.01_0.05", "0.05_0.1", "0.1_0.5", "0.5_2")

RECOMB_BREAKS <- c(0, 0.5, 1, 2, 5, 10, Inf)
RECOMB_LABELS <- c("0-0.5", "0.5-1", "1-2", "2-5", "5-10", ">=10")

MIN_N <- 10

#===============================================================================
# 1. Load recombination rate map
#===============================================================================

cat("Loading deCODE recombination map...\n")

bg <- read.table(RECOMB_BG, header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("chr", "start", "end", "rate"))
bg$chr <- sub("^chr", "", bg$chr)
bg_by_chr <- split(bg, bg$chr)

cat(sprintf("  %d intervals, %d chromosomes\n", nrow(bg), length(bg_by_chr)))

#===============================================================================
# 2. Helper functions
#===============================================================================

annotate_recomb <- function(df) {
    df$recomb_rate <- NA_real_
    for (ch in unique(df$chr)) {
        ch_str <- as.character(ch)
        if (!ch_str %in% names(bg_by_chr)) next
        rc <- bg_by_chr[[ch_str]]
        idx <- which(df$chr == ch)
        ii  <- findInterval(df$pos[idx], rc$start)
        ok  <- ii > 0 & ii <= nrow(rc) & df$pos[idx] <= rc$end[ii]
        df$recomb_rate[idx[ok]] <- rc$rate[ii[ok]]
    }
    df
}

load_pervar <- function(path) {
    if (!file.exists(path) || file.info(path)$size < 50) return(NULL)
    tryCatch({
        d <- read.table(path, header = TRUE, stringsAsFactors = FALSE,
                        col.names = c("chr", "pos", "mean_r2", "n_pairs"))
        if (nrow(d) == 0) return(NULL)
        annotate_recomb(d)
    }, error = function(e) NULL)
}

compare_ld <- function(gc, nongc, label) {
    if (is.null(gc) || is.null(nongc)) return(NULL)
    if (nrow(gc) < MIN_N || nrow(nongc) < MIN_N) return(NULL)
    pv <- tryCatch(wilcox.test(gc$mean_r2, nongc$mean_r2)$p.value,
                   error = function(e) NA)
    data.frame(
        gc_mean       = mean(gc$mean_r2, na.rm = TRUE),
        gc_n          = nrow(gc),
        nongc_mean    = mean(nongc$mean_r2, na.rm = TRUE),
        nongc_n       = nrow(nongc),
        diff          = mean(nongc$mean_r2, na.rm = TRUE) - mean(gc$mean_r2, na.rm = TRUE),
        p_value       = pv,
        gc_recomb     = mean(gc$recomb_rate, na.rm = TRUE),
        nongc_recomb  = mean(nongc$recomb_rate, na.rm = TRUE),
        recomb_bin    = label,
        stringsAsFactors = FALSE
    )
}

#===============================================================================
# 3. Main analysis loop
#===============================================================================

cat(rep("=", 70), "\n", sep = "")
cat("LD vs Recombination Rate Analysis (1:1 Paralogs)\n")
cat(rep("=", 70), "\n\n")

results <- list()

for (mapp in MAPP_CATEGORIES) {
    cat(sprintf("\n=== Mappability: %s ===\n", mapp))

    for (pop in POPULATIONS) {
        cat(sprintf("\n  Population: %s\n", pop))

        for (win in WINDOW_SIZES) {
            cat(sprintf("    Window: %s\n", win))
            win_dir <- file.path(LD_BASE, mapp, pop, win)
            if (!dir.exists(win_dir)) { cat("      [not found]\n"); next }

            for (k in K_VALUES) {
                cat(sprintf("      k=%d: ", k))

                for (maf in MAF_BINS) {
                    gc_path    <- file.path(win_dir, sprintf("ld_k%d_gc_%s_pervar.txt", k, maf))
                    nongc_path <- file.path(win_dir, sprintf("ld_k%d_nongc_%s_pervar.txt", k, maf))

                    gc_all    <- load_pervar(gc_path)
                    nongc_all <- load_pervar(nongc_path)

                    base_row <- data.frame(mappability = mapp, population = pop,
                                           window = win, k = k, maf_bin = maf,
                                           stringsAsFactors = FALSE)

                    # --- overall ---
                    res <- compare_ld(gc_all, nongc_all, "all")
                    if (!is.null(res)) results[[length(results) + 1]] <- cbind(base_row, res)

                    # --- per recomb bin ---
                    if (!is.null(gc_all) && !is.null(nongc_all)) {
                        gc_all$rbin    <- cut(gc_all$recomb_rate,    RECOMB_BREAKS,
                                              labels = RECOMB_LABELS, right = FALSE,
                                              include.lowest = TRUE)
                        nongc_all$rbin <- cut(nongc_all$recomb_rate, RECOMB_BREAKS,
                                              labels = RECOMB_LABELS, right = FALSE,
                                              include.lowest = TRUE)

                        for (rb in RECOMB_LABELS) {
                            gc_sub    <- gc_all[!is.na(gc_all$rbin) & gc_all$rbin == rb, ]
                            nongc_sub <- nongc_all[!is.na(nongc_all$rbin) & nongc_all$rbin == rb, ]
                            res <- compare_ld(gc_sub, nongc_sub, rb)
                            if (!is.null(res)) results[[length(results) + 1]] <- cbind(base_row, res)
                        }
                    }
                    cat(".")
                }
                cat("\n")
            }
        }
    }
}

df <- bind_rows(results)
df$sig <- ""
df$sig[df$p_value < 0.05]  <- "*"
df$sig[df$p_value < 0.01]  <- "**"
df$sig[df$p_value < 0.001] <- "***"
df$maf_bin    <- factor(df$maf_bin, levels = MAF_BINS)
df$recomb_bin <- factor(df$recomb_bin, levels = c("all", RECOMB_LABELS))

write.csv(df, file.path(LD_BASE, "ld_vs_recombrate_1to1.csv"), row.names = FALSE)
cat(sprintf("\nSaved %d rows to ld_vs_recombrate_1to1.csv\n", nrow(df)))

#===============================================================================
# 4. Heatmaps per recomb bin
#===============================================================================

cat("\nGenerating heatmaps...\n")

for (mapp in MAPP_CATEGORIES) {
    for (pop in POPULATIONS) {
        for (win in WINDOW_SIZES) {
            for (rb in c("all", RECOMB_LABELS)) {
                sub <- df %>% filter(mappability == mapp, population == pop,
                                     window == win, recomb_bin == rb)
                if (nrow(sub) == 0) next

                rb_file <- gsub(">=", "gte", rb)

                p <- ggplot(sub, aes(x = maf_bin, y = factor(k), fill = diff * 1000)) +
                    geom_tile(color = "white") +
                    geom_text(aes(label = sprintf("%.1f%s", diff * 1000, sig)), size = 2.5) +
                    scale_fill_gradient2(low = "#E41A1C", mid = "white", high = "#4DAF4A",
                                        midpoint = 0, limits = c(-40, 40),
                                        oob = scales::squish, name = "LD diff\n(\u00d71000)") +
                    labs(title = sprintf("1:1 LD GC\u2013nonGC:  %s | %s | %s | recomb %s cM/Mb",
                                        mapp, pop, win, rb),
                         subtitle = "Green = GC has lower LD  |  Red = GC has higher LD",
                         x = "AF bin", y = "k") +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))

                ggsave(file.path(LD_BASE, sprintf("ld_recomb_heatmap_1to1_%s_%s_%s_%s.png",
                                                   mapp, pop, win, rb_file)),
                       p, width = 12, height = 8, dpi = 150)
            }
        }
    }
}

#===============================================================================
# 5. Trend plots
#===============================================================================

cat("Generating trend plots...\n")

trend_k   <- c(17, 21, 41, 61, 81, 91)
trend_maf <- c("0.001_0.005", "0.005_0.01", "0.01_0.05", "0.05_0.1", "0.1_0.5")

for (mapp in MAPP_CATEGORIES) {
    for (pop in POPULATIONS) {
        for (win in WINDOW_SIZES) {
            td <- df %>%
                filter(mappability == mapp, population == pop, window == win,
                       recomb_bin != "all",
                       k %in% trend_k,
                       maf_bin %in% trend_maf) %>%
                mutate(recomb_bin = factor(recomb_bin, levels = RECOMB_LABELS),
                       k_label = paste0("k=", k))

            if (nrow(td) == 0) next

            p <- ggplot(td, aes(x = recomb_bin, y = diff * 1000,
                                color = factor(k), group = factor(k))) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
                geom_point(size = 2) +
                geom_line() +
                facet_wrap(~ maf_bin, ncol = 5, scales = "free_y") +
                scale_color_brewer(palette = "Set1", name = "k") +
                labs(title = sprintf("1:1 LD diff (nonGC \u2013 GC) vs recomb rate  [%s | %s | %s]",
                                     mapp, pop, win),
                     subtitle = "Positive = GC has lower LD  |  Negative = GC has higher LD",
                     x = "Recombination rate (cM/Mb)", y = "LD diff (\u00d71000)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                      strip.text = element_text(size = 8))

            ggsave(file.path(LD_BASE, sprintf("ld_recomb_trend_1to1_%s_%s_%s.png",
                                               mapp, pop, win)),
                   p, width = 16, height = 6, dpi = 150)
        }
    }
}

#===============================================================================
# 6. Confounding plot
#===============================================================================

cat("Generating confounding plots...\n")

for (mapp in MAPP_CATEGORIES) {
    for (pop in POPULATIONS) {
        for (win in WINDOW_SIZES) {
            cd <- df %>%
                filter(mappability == mapp, population == pop,
                       window == win, recomb_bin == "all") %>%
                select(k, maf_bin, gc_recomb, nongc_recomb) %>%
                pivot_longer(cols = c(gc_recomb, nongc_recomb),
                             names_to = "group", values_to = "mean_recomb") %>%
                mutate(group = ifelse(group == "gc_recomb", "GC", "nonGC"))

            if (nrow(cd) == 0) next

            p <- ggplot(cd, aes(x = maf_bin, y = mean_recomb, fill = group)) +
                geom_col(position = "dodge") +
                facet_wrap(~ paste0("k=", k), ncol = 5, scales = "free_y") +
                scale_fill_manual(values = c("GC" = "#E41A1C", "nonGC" = "#4DAF4A")) +
                labs(title = sprintf("1:1 Mean recomb rate: GC vs nonGC  [%s | %s | %s]",
                                     mapp, pop, win),
                     x = "AF bin", y = "Mean recomb rate (cM/Mb)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                      strip.text = element_text(size = 8))

            ggsave(file.path(LD_BASE, sprintf("ld_recomb_confounding_1to1_%s_%s_%s.png",
                                               mapp, pop, win)),
                   p, width = 16, height = 10, dpi = 150)
        }
    }
}

#===============================================================================
# 7. Matrix tables per recomb bin
#===============================================================================

cat("Generating matrix tables...\n")

for (mapp in MAPP_CATEGORIES) {
    for (pop in POPULATIONS) {
        for (win in WINDOW_SIZES) {
            for (rb in c("all", RECOMB_LABELS)) {
                sub <- df %>% filter(mappability == mapp, population == pop,
                                     window == win, recomb_bin == rb)
                if (nrow(sub) == 0) next

                rb_file <- gsub(">=", "gte", rb)
                out_path <- file.path(LD_BASE,
                    sprintf("ld_recomb_matrix_1to1_%s_%s_%s_%s.tsv", mapp, pop, win, rb_file))
                sink(out_path)

                cat(sprintf("# 1:1 LD GC-nonGC | %s | %s | %s | recomb %s cM/Mb\n",
                            mapp, pop, win, rb))
                cat("# diff = (nonGC_mean - GC_mean) x 1000\n")
                cat("# Positive = GC has lower LD\n\n")

                cat("k\t", paste(MAF_BINS, collapse = "\t"), "\n", sep = "")
                for (kv in K_VALUES) {
                    row <- sapply(MAF_BINS, function(m) {
                        v <- sub %>% filter(k == kv, maf_bin == m) %>% pull(diff)
                        if (length(v) == 0) "." else sprintf("%.1f", v * 1000)
                    })
                    cat(kv, "\t", paste(row, collapse = "\t"), "\n", sep = "")
                }

                cat("\n# Significance (* p<0.05  ** p<0.01  *** p<0.001)\n")
                cat("k\t", paste(MAF_BINS, collapse = "\t"), "\n", sep = "")
                for (kv in K_VALUES) {
                    row <- sapply(MAF_BINS, function(m) {
                        v <- sub %>% filter(k == kv, maf_bin == m) %>% pull(sig)
                        if (length(v) == 0) "." else ifelse(v == "", "-", v)
                    })
                    cat(kv, "\t", paste(row, collapse = "\t"), "\n", sep = "")
                }

                cat("\n# N (GC / nonGC)\n")
                cat("k\t", paste(MAF_BINS, collapse = "\t"), "\n", sep = "")
                for (kv in K_VALUES) {
                    row <- sapply(MAF_BINS, function(m) {
                        r <- sub %>% filter(k == kv, maf_bin == m)
                        if (nrow(r) == 0) "." else sprintf("%d/%d", r$gc_n, r$nongc_n)
                    })
                    cat(kv, "\t", paste(row, collapse = "\t"), "\n", sep = "")
                }

                sink()
                cat(sprintf("  %s\n", basename(out_path)))
            }
        }
    }
}

cat("\nDone. Output in:", LD_BASE, "\n")
