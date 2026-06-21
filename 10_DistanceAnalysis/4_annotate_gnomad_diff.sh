#!/bin/bash
# =============================================================================
# Step 4: Annotate diff pairs with gnomAD (AF + ref/alt)
# =============================================================================
#
# For each pair in distances.diff.txt.gz, look up both positions in gnomAD
# and check for the EXACT gene-conversion-consistent variant:
#   - pos1 (base1) paired with pos2 (base2):
#     look for variant at pos1 with ref=base1, alt=base2
#   - pos2 (base2) paired with pos1 (base1):
#     look for variant at pos2 with ref=base2, alt=base1
#
# Input:  distance_analysis/{REPLEN}/distances.diff.txt.gz
# Output: distance_analysis/{REPLEN}/distances.diff.gnomad.txt.gz
#
# Output columns:
#   chr1 pos1 base1 chr2 pos2 base2 distance type AF1 ref1 alt1 AF2 ref2 alt2
#
# AF1/ref1/alt1: gnomAD variant at pos1 matching base1→base2 change
# AF2/ref2/alt2: gnomAD variant at pos2 matching base2→base1 change
# No matching variant → AF=0, ref=".", alt="."
# =============================================================================

: "${PROJECT_ROOT:?PROJECT_ROOT must be set (see 00_Configuration/README.md)}"
BASEDIR="${PROJECT_ROOT}/geneconv_complete"
GNOMAD="${BASEDIR}/gnomad/genome"

# Parse arguments
for arg in "$@"; do
    case $arg in
        --REPLEN=*) REPLEN="${arg#*=}" ;;
    esac
done

if [ -z "$REPLEN" ]; then
    echo "Usage: $0 --REPLEN=<value>"
    exit 1
fi

INFILE="${BASEDIR}/distance_analysis/${REPLEN}/distances.diff.txt.gz"
OUTFILE="${BASEDIR}/distance_analysis/${REPLEN}/distances.diff.gnomad.txt.gz"
TMPDIR="${BASEDIR}/distance_analysis/${REPLEN}/tmp_gnomad_diff"

if [ ! -f "$INFILE" ]; then
    echo "ERROR: ${INFILE} does not exist"
    exit 1
fi

mkdir -p "$TMPDIR"

echo "=== Step 4: gnomAD annotation for k=${REPLEN} diff pairs ==="

# --- Extract unique positions ---
echo "Phase 1: Extracting unique positions..."

zcat "$INFILE" | awk -F'\t' 'BEGIN{OFS="\t"}{
    print $1, $2
    print $4, $5
}' | sort -u -k1,1V -k2,2n > "${TMPDIR}/all_positions.txt"

n_pos=$(wc -l < "${TMPDIR}/all_positions.txt")
echo "  Unique positions: ${n_pos}"

# --- Query gnomAD per chromosome ---
echo "Phase 2: Looking up positions in gnomAD..."

> "${TMPDIR}/gnomad_hits.txt"

for CHR in $(cut -f1 "${TMPDIR}/all_positions.txt" | sort -uV); do
    GNOMAD_FILE="${GNOMAD}/${CHR}.sorted.bed.gz"

    if [ ! -f "$GNOMAD_FILE" ]; then
        echo "  WARNING: No gnomAD file for chr ${CHR}"
        continue
    fi

    # Create BED regions for tabix (0-based start, 1-based end)
    awk -v chr="$CHR" '$1 == chr {print chr "\t" ($2-1) "\t" $2}' \
        "${TMPDIR}/all_positions.txt" > "${TMPDIR}/regions_${CHR}.bed"

    n_regions=$(wc -l < "${TMPDIR}/regions_${CHR}.bed")

    if [ "$n_regions" -gt 0 ]; then
        tabix "$GNOMAD_FILE" -R "${TMPDIR}/regions_${CHR}.bed" \
            >> "${TMPDIR}/gnomad_hits.txt" 2>/dev/null
        echo "  chr${CHR}: queried ${n_regions} positions"
    fi
done

n_hits=$(wc -l < "${TMPDIR}/gnomad_hits.txt")
echo "  Total gnomAD hits: ${n_hits}"

# --- Build AF lookup keyed by chr:pos:ref:alt ---
echo "Phase 3: Building allele-specific AF lookup..."

# gnomAD BED: chr  start  end  variant_id  AF  AC  AN
# variant_id format: CHR_POS_POS_REF/ALT
# Key: chr:pos:ref:alt → AF (keeps ALL allele combinations)
awk -F'\t' '{
    # Parse ref/alt from variant_id (e.g., 1_10111_10111_C/A)
    n = split($4, vid, "_")
    split(vid[n], alleles, "/")
    ref = alleles[1]
    alt = alleles[2]
    af = $5 + 0

    key = $1 ":" $3 ":" ref ":" alt
    if (!(key in maxaf) || af > maxaf[key]) {
        maxaf[key] = af
    }
}
END {
    for (key in maxaf)
        print key "\t" maxaf[key]
}' "${TMPDIR}/gnomad_hits.txt" > "${TMPDIR}/af_lookup.txt"

n_lookup=$(wc -l < "${TMPDIR}/af_lookup.txt")
echo "  Unique chr:pos:ref:alt entries: ${n_lookup}"

# --- Annotate pairs with gene-conversion-consistent variants ---
echo "Phase 4: Annotating pairs (allele-matched)..."

# For each pair (chr1,pos1,base1 -- chr2,pos2,base2):
#   pos1: GC-consistent variant = ref=base1, alt=base2 → look up chr1:pos1:base1:base2
#   pos2: GC-consistent variant = ref=base2, alt=base1 → look up chr2:pos2:base2:base1
awk -F'\t' '
NR == FNR {
    af[$1] = $2
    next
}
{
    key1 = $1 ":" $2 ":" $3 ":" $6    # chr1:pos1:base1:base2
    key2 = $4 ":" $5 ":" $6 ":" $3    # chr2:pos2:base2:base1

    af1 = (key1 in af) ? af[key1] : 0
    r1  = (af1 > 0) ? $3 : "."
    a1  = (af1 > 0) ? $6 : "."

    af2 = (key2 in af) ? af[key2] : 0
    r2  = (af2 > 0) ? $6 : "."
    a2  = (af2 > 0) ? $3 : "."

    print $0 "\t" af1 "\t" r1 "\t" a1 "\t" af2 "\t" r2 "\t" a2
}' "${TMPDIR}/af_lookup.txt" <(zcat "$INFILE") | gzip > "$OUTFILE"

# --- Cleanup ---
rm -rf "$TMPDIR"

# --- Summary ---
total=$(zcat "$OUTFILE" | wc -l)
both_found=$(zcat "$OUTFILE" | awk -F'\t' '$9+0 > 0 && $12+0 > 0' | wc -l)
one_found=$(zcat "$OUTFILE" | awk -F'\t' '($9+0 > 0 && $12+0 == 0) || ($9+0 == 0 && $12+0 > 0)' | wc -l)
neither=$(zcat "$OUTFILE" | awk -F'\t' '$9+0 == 0 && $12+0 == 0' | wc -l)

echo ""
echo "=== Summary for k=${REPLEN} ==="
echo "Total diff pairs:             ${total}"
echo "Both positions in gnomAD:     ${both_found}"
echo "Only one position in gnomAD:  ${one_found}"
echo "Neither in gnomAD:            ${neither}"

# --- Transition counts ---
echo ""
echo "=== Nucleotide transitions (positions with variant) ==="
echo "--- Position 1 ---"
zcat "$OUTFILE" | awk -F'\t' '$9+0 > 0 && $10 != "." && $11 != "." {
    print $10 ">" $11
}' | sort | uniq -c | sort -rn | head -12

echo "--- Position 2 ---"
zcat "$OUTFILE" | awk -F'\t' '$12+0 > 0 && $13 != "." && $14 != "." {
    print $13 ">" $14
}' | sort | uniq -c | sort -rn | head -12

# --- GC-bias summary ---
echo ""
echo "=== GC-bias direction ==="
zcat "$OUTFILE" | awk -F'\t' '
function gc_class(ref, alt) {
    if ((ref == "A" || ref == "T") && (alt == "G" || alt == "C")) return "WS"
    if ((ref == "G" || ref == "C") && (alt == "A" || alt == "T")) return "SW"
    if ((ref == "A" || ref == "T") && (alt == "A" || alt == "T")) return "WW"
    if ((ref == "G" || ref == "C") && (alt == "G" || alt == "C")) return "SS"
    return "other"
}
{
    if ($9+0 > 0 && $10 != ".") {
        cl = gc_class($10, $11)
        counts[cl]++
        total++
    }
    if ($12+0 > 0 && $13 != ".") {
        cl = gc_class($13, $14)
        counts[cl]++
        total++
    }
}
END {
    printf "  W→S (GC-compatible): %d (%.1f%%)\n", counts["WS"]+0, 100*(counts["WS"]+0)/total
    printf "  S→W (anti-GC):       %d (%.1f%%)\n", counts["SW"]+0, 100*(counts["SW"]+0)/total
    printf "  S→S:                 %d (%.1f%%)\n", counts["SS"]+0, 100*(counts["SS"]+0)/total
    printf "  W→W:                 %d (%.1f%%)\n", counts["WW"]+0, 100*(counts["WW"]+0)/total
}'

echo ""
echo "Done. Output: ${OUTFILE}"
