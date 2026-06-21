#!/bin/bash
# =============================================================================
# Step 3: Compute pair-level gcdiffs for a specific k (v2 — simplified)
# =============================================================================
#
# For a given k, find pairs specific to this identity level:
#   1. Build set of ALL positions appearing at ANY higher k-value
#   2. Filter: keep only pairs where NEITHER position is in that set
#
# The position filter automatically handles pair-level diffs:
#   if a pair exists at higher k, both its positions are at higher k.
#
# No separate maxk build step needed. No per-chromosome splitting.
# One awk pass for building the set, one awk pass for filtering.
#
# Memory: dominated by higher-k positions set
#   k=91: 0 (just copy)
#   k=81-21: <1.5GB
#   k=19: ~1.3GB
#   k=17: ~10GB (k=19 has 120M pairs → ~129M unique positions)
#
# Input:  distance_analysis/{k}/distances.txt.gz
#         distance_analysis/{higher_k}/distances.txt.gz for all higher k
# Output: distance_analysis/{k}/distances.diff.txt.gz
#
# =============================================================================

: "${PROJECT_ROOT:?PROJECT_ROOT must be set (see 00_Configuration/README.md)}"
BASEDIR="${PROJECT_ROOT}/geneconv_complete"

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

INFILE="${BASEDIR}/distance_analysis/${REPLEN}/distances.txt.gz"
OUTFILE="${BASEDIR}/distance_analysis/${REPLEN}/distances.diff.txt.gz"

if [ ! -f "$INFILE" ]; then
    echo "ERROR: $INFILE not found"
    exit 1
fi

total_in=$(zcat "$INFILE" | wc -l)
echo "k=${REPLEN}: ${total_in} total pairs"

# --- Determine all higher k-values ---
ALL_K="91 81 71 61 51 41 31 21 19 17"
higher_ks=""
for kk in $ALL_K; do
    [ "$kk" -eq "$REPLEN" ] && break
    higher_ks="$higher_ks $kk"
done

if [ -z "$higher_ks" ]; then
    # k=91: highest level, no filtering needed
    echo "k=91: highest level, copying input directly."
    cp "$INFILE" "$OUTFILE"
else
    echo "Higher k-values: $higher_ks"

    # --- Step 1: Build higher-k positions set ---
    echo "Building higher-k positions set..."
    HIGHER_TMP=$(mktemp "${BASEDIR}/distance_analysis/higher_pos_${REPLEN}_XXXXXX")
    trap "rm -f $HIGHER_TMP" EXIT

    {
        for hk in $higher_ks; do
            HK_FILE="${BASEDIR}/distance_analysis/${hk}/distances.txt.gz"
            if [ ! -f "$HK_FILE" ]; then
                echo "  WARNING: $HK_FILE not found, skipping" >&2
                continue
            fi
            echo "  Reading k=${hk}..." >&2
            zcat "$HK_FILE" | awk -F'\t' '{print $1":"$2; print $4":"$5}'
        done
    } | awk '!seen[$0]++' > "$HIGHER_TMP"

    n_higher=$(wc -l < "$HIGHER_TMP")
    echo "  Unique higher-k positions: ${n_higher}"

    # --- Step 2: Filter pairs ---
    echo "Filtering pairs..."
    zcat "$INFILE" | awk -F'\t' '
    NR == FNR { higher[$0] = 1; next }
    {
        if (!($1 ":" $2 in higher) && !($4 ":" $5 in higher)) print
    }' "$HIGHER_TMP" - | gzip > "$OUTFILE"
fi

# --- Summary ---
total_out=$(zcat "$OUTFILE" | wc -l)
intra_out=$(zcat "$OUTFILE" | awk -F'\t' '$8=="intra"' | wc -l)
inter_out=$(zcat "$OUTFILE" | awk -F'\t' '$8=="inter"' | wc -l)

pct=$(echo "scale=1; 100 * $total_out / $total_in" | bc)
echo ""
echo "k=${REPLEN}: ${total_in} total → ${total_out} diff pairs (${pct}%)"
echo "  intra-chromosomal: ${intra_out}"
echo "  inter-chromosomal: ${inter_out}"

if [ "$intra_out" -gt 0 ]; then
    echo ""
    echo "Diff pair distance distribution (intra):"
    zcat "$OUTFILE" | awk -F'\t' '$8=="intra" {
        d = $7
        if      (d < 1000)        bin = "<1kb"
        else if (d < 10000)       bin = "1-10kb"
        else if (d < 100000)      bin = "10-100kb"
        else if (d < 1000000)     bin = "100kb-1Mb"
        else if (d < 10000000)    bin = "1-10Mb"
        else if (d < 100000000)   bin = "10-100Mb"
        else                      bin = ">100Mb"
        counts[bin]++
    }
    END {
        split("<1kb,1-10kb,10-100kb,100kb-1Mb,1-10Mb,10-100Mb,>100Mb", order, ",")
        for (i = 1; i <= 7; i++) {
            if (counts[order[i]]) printf "  %-14s %d\n", order[i], counts[order[i]]
        }
    }'
fi

echo ""
echo "Done. Output: ${OUTFILE}"
