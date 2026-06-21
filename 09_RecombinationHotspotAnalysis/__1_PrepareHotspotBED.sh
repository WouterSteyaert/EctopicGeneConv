#!/bin/bash
#===============================================================================
# __1_PrepareHotspotBED.sh
#
# Derive the recombination hotspot BED from the deCODE sex-averaged
# recombination map (UCSC recombAvg bigWig).  Hotspots are defined as regions
# with rate >= 10 cM/Mb (Halldorsson et al. 2019; ~10x the genome-wide
# background of ~1 cM/Mb).
#
# Procedure (matches the working pipeline in FetchRegions/__4_PrepareBedFiles.pl):
#   1. bigWigToBedGraph $RECOMB_BW -> recombAvg.bedGraph
#   2. awk '$4 >= 10' (strip 'chr' prefix)
#   3. bedtools sort
#   4. bedtools merge   (collapse adjacent/overlapping hotspot intervals)
#   5. bgzip + tabix
#
# Requires: bigWigToBedGraph (UCSC utility), bedtools, bgzip, tabix.
#
# Inputs (env):
#   PROJECT_ROOT
#   RECOMB_BW             (optional)  override path to recombAvg.bw; defaults
#                                     to $PROJECT_ROOT/geneconv_complete/bed/deCODE/recombAvg.bw
#   HOTSPOT_THRESHOLD     (optional)  override cM/Mb threshold; default 10
#
# Output:
#   $PROJECT_ROOT/geneconv_complete/bed/deCODE/RecombinationHotspots.sorted.bed.gz (+ .tbi)
#===============================================================================
set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT must be set (see 00_Configuration/README.md)}"

DEC_DIR="${PROJECT_ROOT}/geneconv_complete/bed/deCODE"
RECOMB_BW="${RECOMB_BW:-${DEC_DIR}/recombAvg.bw}"
RECOMB_BG="${DEC_DIR}/recombAvg.bedGraph"
THRESHOLD="${HOTSPOT_THRESHOLD:-10}"

TEMP_HOT="${DEC_DIR}/hotspots_temp.bed"
TEMP_SRT="${DEC_DIR}/hotspots_sorted_temp.bed"
HOT_BED="${DEC_DIR}/RecombinationHotspots.bed"
HOT_SORTED="${DEC_DIR}/RecombinationHotspots.sorted.bed"

[ -f "$RECOMB_BW" ] || { echo "ERROR: $RECOMB_BW not found"; exit 1; }
mkdir -p "$DEC_DIR"

echo "Step 1/4: bigWigToBedGraph..."
bigWigToBedGraph "$RECOMB_BW" "$RECOMB_BG"

echo "Step 2/4: thresholding at >= ${THRESHOLD} cM/Mb..."
awk -v t="$THRESHOLD" 'BEGIN{OFS="\t"} $4 >= t {
    sub(/^chr/, "", $1)
    print $1, $2, $3
}' "$RECOMB_BG" > "$TEMP_HOT"

echo "Step 3/4: bedtools sort + merge..."
bedtools sort -i "$TEMP_HOT" > "$TEMP_SRT"
bedtools merge -i "$TEMP_SRT" > "$HOT_BED"
bedtools sort -i "$HOT_BED" > "$HOT_SORTED"

echo "Step 4/4: bgzip + tabix..."
bgzip -f "$HOT_SORTED"
tabix -p bed "${HOT_SORTED}.gz"

rm -f "$TEMP_HOT" "$TEMP_SRT"

n=$(zcat "${HOT_SORTED}.gz" | wc -l)
total=$(zcat "${HOT_SORTED}.gz" | awk '{s+=$3-$2} END{print s}')
echo ""
echo "Done."
echo "  Output:    ${HOT_SORTED}.gz"
echo "  Intervals: $n  (expected ~70,000)"
echo "  Total bp:  $total  (~2% of genome)"
