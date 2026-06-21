#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=0:30:00

# =============================================================================
# Step 3: Build "excluding peaks" supplementary table (Supplementary Table 12).
#
# Thin wrapper around 3_excl_peaks_table.R, which, for each k x AF bin over the
# non-overlapping 1-Mb windows, reports mean log2(OR) and % positive windows
# (all windows and after peak exclusion) AND a genome-wide chi-square enrichment
# test (Bonferroni-corrected log10 p, m = 90) before and after exclusion, so the
# "enrichment survives peak exclusion" statement is backed by a stored statistic.
#
# Input:  $PROJECT_ROOT/geneconv_complete/sliding_window_corrected/
#             enrichment_corrected_combined.tsv, peaks_corrected.tsv
# Output: same dir / SupTable_EnrichmentExclPeaks_corrected.tsv
#
# Requires R (conda env 'homolo' provides R 4.2/4.3.2).
# =============================================================================

: "${PROJECT_ROOT:?PROJECT_ROOT must be set (see 00_Configuration/README.md)}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Rscript "${SCRIPT_DIR}/3_excl_peaks_table.R"
