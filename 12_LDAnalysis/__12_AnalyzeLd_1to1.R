#!/usr/bin/env Rscript
#===============================================================================
# LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 12
# Analyse mean R² per variant (1:1 paralogs only)
#
# Output: LD_analysis_resampling_1to1/
#===============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

#===============================================================================
# Configuration
#===============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot("PROJECT_ROOT env var must be set" = nzchar(project_root))
BASE    <- file.path(project_root, "geneconv_complete")
ld_base <- file.path(BASE, "LD_analysis_resampling_1to1")

mapp_categories <- c("allmapp", "nosegdupmapp", "segdupmapp")
populations <- c("EUR", "ALL")
k_values <- c(17, 19, 21, 31, 41, 51, 61, 71, 81, 91)
window_sizes <- c("win10kb", "win25kb", "win50kb", "win100kb")

maf_bins <- c("0.001_0.005", "0.005_0.01", "0.01_0.05", "0.05_0.1", "0.1_0.5", "0.5_2")

#===============================================================================
# Helper functions
#===============================================================================

load_pervar_data <- function(file_path) {
    if (!file.exists(file_path)) return(NULL)
    if (file.info(file_path)$size < 50) return(NULL)

    tryCatch({
        df <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE,
                        col.names = c("chr", "pos", "mean_r2", "n_pairs"))
        if (nrow(df) == 0) return(NULL)
        return(df)
    }, error = function(e) NULL)
}

#===============================================================================
# Main Analysis
#===============================================================================

cat(rep("=", 70), "\n")
cat("LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 12\n")
cat(rep("=", 70), "\n\n")

all_results <- list()

for (mapp in mapp_categories) {
    cat(sprintf("\n=== Mappability: %s ===\n", mapp))

    for (pop in populations) {
        cat(sprintf("\n  Population: %s\n", pop))

        for (win in window_sizes) {
            cat(sprintf("    Window: %s\n", win))
            win_dir <- file.path(ld_base, mapp, pop, win)

            if (!dir.exists(win_dir)) {
                cat("      [Directory not found]\n")
                next
            }

            for (k in k_values) {
                cat(sprintf("      k = %d: ", k))

                for (maf_bin in maf_bins) {
                    gc_file <- file.path(win_dir, sprintf("ld_k%d_gc_%s_pervar.txt", k, maf_bin))
                    nongc_file <- file.path(win_dir, sprintf("ld_k%d_nongc_%s_pervar.txt", k, maf_bin))

                    gc_data <- load_pervar_data(gc_file)
                    nongc_data <- load_pervar_data(nongc_file)

                    if (is.null(gc_data) || is.null(nongc_data)) next
                    if (nrow(gc_data) < 30 || nrow(nongc_data) < 30) next

                    gc_mean <- mean(gc_data$mean_r2, na.rm = TRUE)
                    nongc_mean <- mean(nongc_data$mean_r2, na.rm = TRUE)

                    test_result <- tryCatch({
                        wilcox.test(gc_data$mean_r2, nongc_data$mean_r2)$p.value
                    }, error = function(e) NA)

                    all_results[[length(all_results) + 1]] <- data.frame(
                        mappability = mapp,
                        population = pop,
                        window = win,
                        k = k,
                        maf_bin = maf_bin,
                        gc_mean = gc_mean,
                        gc_n = nrow(gc_data),
                        nongc_mean = nongc_mean,
                        nongc_n = nrow(nongc_data),
                        diff = nongc_mean - gc_mean,
                        p_value = test_result,
                        direction = ifelse(gc_mean < nongc_mean, "GC < Non-GC", "GC >= Non-GC")
                    )
                }
                cat(".")
            }
            cat("\n")
        }
    }
}

if (length(all_results) == 0) {
    stop("Geen resultaten gevonden. Check of de LD berekening voltooid is.")
}

results_df <- bind_rows(all_results)

results_df$sig <- ""
results_df$sig[results_df$p_value < 0.05] <- "*"
results_df$sig[results_df$p_value < 0.01] <- "**"
results_df$sig[results_df$p_value < 0.001] <- "***"

results_df$maf_bin <- factor(results_df$maf_bin, levels = maf_bins)
results_df$window <- factor(results_df$window, levels = window_sizes)

write.csv(results_df, file.path(ld_base, "ld_analysis_results_1to1.csv"), row.names = FALSE)
cat(sprintf("\nTotaal %d resultaten opgeslagen\n", nrow(results_df)))

#===============================================================================
# Create Matrix Tables
#===============================================================================

cat("\n\nMatrix tabellen...\n")

for (win in window_sizes) {
    for (mapp in mapp_categories) {
        for (pop in populations) {
            subset_data <- results_df %>%
                filter(mappability == mapp, population == pop, window == win)

            if (nrow(subset_data) == 0) next

            output_file <- file.path(ld_base, sprintf("ld_matrix_1to1_%s_%s_%s.tsv", mapp, pop, win))
            sink(output_file)

            cat(sprintf("# LD Analysis 1:1 Paralogs (RESAMPLING): %s, %s, %s\n", mapp, pop, win))
            cat("# LD Difference (Non-GC - GC) x 1000\n")
            cat("# Positive = GC has LOWER LD\n\n")
            cat("k\t", paste(maf_bins, collapse = "\t"), "\n", sep = "")

            for (k in k_values) {
                row_data <- sapply(maf_bins, function(bin) {
                    val <- subset_data %>% filter(k == !!k, maf_bin == bin) %>% pull(diff)
                    if (length(val) == 0) return("NA")
                    return(sprintf("%.2f", val * 1000))
                })
                cat(k, "\t", paste(row_data, collapse = "\t"), "\n", sep = "")
            }

            cat("\n# Significance\n")
            cat("k\t", paste(maf_bins, collapse = "\t"), "\n", sep = "")
            for (k in k_values) {
                row_data <- sapply(maf_bins, function(bin) {
                    val <- subset_data %>% filter(k == !!k, maf_bin == bin) %>% pull(sig)
                    if (length(val) == 0) return("NA")
                    return(ifelse(val == "", "-", val))
                })
                cat(k, "\t", paste(row_data, collapse = "\t"), "\n", sep = "")
            }

            cat("\n# N variants GC / Non-GC\n")
            cat("k\t", paste(maf_bins, collapse = "\t"), "\n", sep = "")
            for (k in k_values) {
                row_data <- sapply(maf_bins, function(bin) {
                    row <- subset_data %>% filter(k == !!k, maf_bin == bin)
                    if (nrow(row) == 0) return("NA")
                    return(sprintf("%d/%d", row$gc_n, row$nongc_n))
                })
                cat(k, "\t", paste(row_data, collapse = "\t"), "\n", sep = "")
            }

            sink()
            cat(sprintf("  Saved: %s\n", basename(output_file)))
        }
    }
}

#===============================================================================
# Visualizations
#===============================================================================

cat("\n\nVisualisaties...\n")

for (win in window_sizes) {
    for (mapp in mapp_categories) {
        for (pop in populations) {
            subset_data <- results_df %>%
                filter(mappability == mapp, population == pop, window == win)

            if (nrow(subset_data) == 0) next

            p1 <- ggplot(subset_data, aes(x = maf_bin, y = factor(k), fill = diff * 1000)) +
                geom_tile(color = "white") +
                geom_text(aes(label = sprintf("%.1f%s", diff * 1000, sig)), size = 2.5) +
                scale_fill_gradient2(low = "#E41A1C", mid = "white", high = "#4DAF4A",
                                    midpoint = 0, limits = c(-20, 20), oob = scales::squish) +
                labs(
                    title = sprintf("1:1 PARALOGS RESAMPLING: %s, %s, %s", mapp, pop, win),
                    x = "AF bin", y = "k"
                ) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))

            ggsave(file.path(ld_base, sprintf("ld_heatmap_1to1_%s_%s_%s.png", mapp, pop, win)),
                   p1, width = 12, height = 8, dpi = 150)
        }
    }
}

cat("\nResultaten opgeslagen in: ", ld_base, "\n")
