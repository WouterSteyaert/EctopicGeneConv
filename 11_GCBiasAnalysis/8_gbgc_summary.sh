#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=4:00:00

# =============================================================================
# Step 7: gBGC (GC-biased gene conversion) summary tables
# =============================================================================
#
# The existing WS/SW COUNT ratio (plot 09) reflects mutation bias, not gBGC.
# S→W mutations are ~2x more frequent than W→S, so the count ratio is always <1.
#
# The proper gBGC signal is in ALLELE FREQUENCY: WS variants reach 2-3x higher
# AF than SW variants, because gBGC favors fixation of GC alleles.
#
# This script produces 4 summary tables from existing gnomAD-annotated files.
#
# Input:  distance_analysis/{k}/distances.diff.gnomad.txt.gz
# Output: distance_analysis/summary/gbgc_*.tsv
#
# Columns in input:
#   1:chr1  2:pos1  3:base1  4:chr2  5:pos2  6:base2
#   7:distance  8:type  9:AF1  10:ref1  11:alt1  12:AF2  13:ref2  14:alt2
# =============================================================================

source ~/.bashrc
conda activate homolo

if [[ -z "$PROJECT_ROOT" ]]; then
    echo "ERROR: PROJECT_ROOT env var must be set" >&2
    exit 1
fi
BASEDIR="${PROJECT_ROOT}/geneconv_complete"
OUTDIR="${BASEDIR}/distance_analysis/summary"
mkdir -p "$OUTDIR"

K_VALUES="17 19 21 31 41 51 61 71 81 91"

# =============================================================================
# Table 1: Mean and median AF of WS vs SW per k
# =============================================================================
echo "=== Table 1: gBGC AF by k ==="

OUTFILE1="${OUTDIR}/gbgc_af_by_k.tsv"
echo -e "k\tdirection\tn\tmean_AF\tmedian_AF\tsum_AF" > "$OUTFILE1"

for k in $K_VALUES; do
    INFILE="${BASEDIR}/distance_analysis/${k}/distances.diff.gnomad.txt.gz"
    [ ! -f "$INFILE" ] && echo "  WARNING: $INFILE missing, skipping" && continue
    echo "  Processing k=$k..."

    # Extract all WS and SW variants with their AF
    TMPFILE=$(mktemp)
    zcat "$INFILE" | awk -F'\t' '
    function gc_class(ref, alt) {
        if ((ref == "A" || ref == "T") && (alt == "G" || alt == "C")) return "WS"
        if ((ref == "G" || ref == "C") && (alt == "A" || alt == "T")) return "SW"
        return ""
    }
    {
        if ($9+0 > 0 && $10 != "." && $11 != ".") {
            cl = gc_class($10, $11)
            if (cl != "") print cl "\t" $9+0
        }
        if ($12+0 > 0 && $13 != "." && $14 != ".") {
            cl = gc_class($13, $14)
            if (cl != "") print cl "\t" $12+0
        }
    }' > "$TMPFILE"

    # Compute mean + sum per direction (single awk pass)
    # Then compute median per direction (requires sorted values)
    for dir in WS SW; do
        DIRFILE=$(mktemp)
        awk -F'\t' -v d="$dir" '$1==d {print $2+0}' "$TMPFILE" | sort -g > "$DIRFILE"
        n_dir=$(wc -l < "$DIRFILE")
        if [ "$n_dir" -gt 0 ]; then
            read mean_af sum_af <<< $(awk '{s+=$1} END {printf "%.8f %.4f", s/NR, s}' "$DIRFILE")
            mid=$(( (n_dir + 1) / 2 ))
            median_af=$(sed -n "${mid}p" "$DIRFILE")
            printf "%s\t%s\t%d\t%s\t%s\t%s\n" "$k" "$dir" "$n_dir" "$mean_af" "$median_af" "$sum_af" >> "$OUTFILE1"
        fi
        rm -f "$DIRFILE"
    done

    rm -f "$TMPFILE"
done

echo "  Written: $OUTFILE1"

# =============================================================================
# Table 2: Mean AF WS vs SW per k × distance bin
# =============================================================================
echo ""
echo "=== Table 2: gBGC AF by k x distance ==="

OUTFILE2="${OUTDIR}/gbgc_af_by_k_distance.tsv"
echo -e "k\tdist_bin\tdist_order\tdirection\tn\tmean_AF" > "$OUTFILE2"

for k in $K_VALUES; do
    INFILE="${BASEDIR}/distance_analysis/${k}/distances.diff.gnomad.txt.gz"
    [ ! -f "$INFILE" ] && continue
    echo "  Processing k=$k..."

    zcat "$INFILE" | awk -F'\t' -v k="$k" '
    function dist_bin(d) {
        if (d == "NA")            return "inter"
        else if (d+0 < 1000)     return "<1kb"
        else if (d+0 < 10000)    return "1-10kb"
        else if (d+0 < 100000)   return "10-100kb"
        else if (d+0 < 1000000)  return "100kb-1Mb"
        else if (d+0 < 10000000) return "1-10Mb"
        else if (d+0 < 1e8)      return "10-100Mb"
        else                     return ">100Mb"
    }
    function dist_order(b) {
        if (b == "<1kb") return 1; if (b == "1-10kb") return 2
        if (b == "10-100kb") return 3; if (b == "100kb-1Mb") return 4
        if (b == "1-10Mb") return 5; if (b == "10-100Mb") return 6
        if (b == ">100Mb") return 7; if (b == "inter") return 8
        return 9
    }
    function gc_class(ref, alt) {
        if ((ref == "A" || ref == "T") && (alt == "G" || alt == "C")) return "WS"
        if ((ref == "G" || ref == "C") && (alt == "A" || alt == "T")) return "SW"
        return ""
    }
    {
        b = ($8 == "inter") ? "inter" : dist_bin($7)
        bo = dist_order(b)

        if ($9+0 > 0 && $10 != "." && $11 != ".") {
            cl = gc_class($10, $11)
            if (cl != "") {
                key = b "\t" bo "\t" cl
                n[key]++
                sum_af[key] += $9+0
            }
        }
        if ($12+0 > 0 && $13 != "." && $14 != ".") {
            cl = gc_class($13, $14)
            if (cl != "") {
                key = b "\t" bo "\t" cl
                n[key]++
                sum_af[key] += $12+0
            }
        }
    }
    END {
        for (key in n) {
            printf "%s\t%s\t%d\t%.8f\n", k, key, n[key], sum_af[key] / n[key]
        }
    }' >> "$OUTFILE2"
done

{ head -1 "$OUTFILE2"; tail -n +2 "$OUTFILE2" | sort -t$'\t' -k1,1n -k3,3n -k4,4; } > "${OUTFILE2}.tmp"
mv "${OUTFILE2}.tmp" "$OUTFILE2"
echo "  Written: $OUTFILE2"

# =============================================================================
# Table 3: %WS among WS+SW at AF thresholds
# =============================================================================
echo ""
echo "=== Table 3: gBGC proportion by AF threshold ==="

OUTFILE3="${OUTDIR}/gbgc_proportion_by_af.tsv"
echo -e "k\taf_threshold\tn_WS\tn_SW\tn_total\tpct_WS" > "$OUTFILE3"

for k in $K_VALUES; do
    INFILE="${BASEDIR}/distance_analysis/${k}/distances.diff.gnomad.txt.gz"
    [ ! -f "$INFILE" ] && continue
    echo "  Processing k=$k..."

    zcat "$INFILE" | awk -F'\t' -v k="$k" '
    function gc_class(ref, alt) {
        if ((ref == "A" || ref == "T") && (alt == "G" || alt == "C")) return "WS"
        if ((ref == "G" || ref == "C") && (alt == "A" || alt == "T")) return "SW"
        return ""
    }
    BEGIN {
        # AF thresholds
        split("0,0.0001,0.001,0.005,0.01,0.05,0.1,0.5", thresholds, ",")
        n_thresh = 8
    }
    {
        # Position 1
        if ($9+0 > 0 && $10 != "." && $11 != ".") {
            cl = gc_class($10, $11)
            af = $9+0
            if (cl == "WS") {
                for (t = 1; t <= n_thresh; t++) {
                    if (af >= thresholds[t]+0) ws[t]++
                }
            } else if (cl == "SW") {
                for (t = 1; t <= n_thresh; t++) {
                    if (af >= thresholds[t]+0) sw[t]++
                }
            }
        }
        # Position 2
        if ($12+0 > 0 && $13 != "." && $14 != ".") {
            cl = gc_class($13, $14)
            af = $12+0
            if (cl == "WS") {
                for (t = 1; t <= n_thresh; t++) {
                    if (af >= thresholds[t]+0) ws[t]++
                }
            } else if (cl == "SW") {
                for (t = 1; t <= n_thresh; t++) {
                    if (af >= thresholds[t]+0) sw[t]++
                }
            }
        }
    }
    END {
        for (t = 1; t <= n_thresh; t++) {
            total = ws[t] + sw[t] + 0
            pct = (total > 0) ? 100 * ws[t] / total : 0
            printf "%s\t%s\t%d\t%d\t%d\t%.4f\n", k, thresholds[t], ws[t]+0, sw[t]+0, total, pct
        }
    }' >> "$OUTFILE3"
done

echo "  Written: $OUTFILE3"

# =============================================================================
# Table 4: AF spectrum in log-spaced bins, WS vs SW
# =============================================================================
echo ""
echo "=== Table 4: gBGC AF spectrum ==="

OUTFILE4="${OUTDIR}/gbgc_af_spectrum.tsv"
echo -e "k\tdirection\taf_bin\taf_order\tn\tpct_of_direction" > "$OUTFILE4"

for k in $K_VALUES; do
    INFILE="${BASEDIR}/distance_analysis/${k}/distances.diff.gnomad.txt.gz"
    [ ! -f "$INFILE" ] && continue
    echo "  Processing k=$k..."

    zcat "$INFILE" | awk -F'\t' -v k="$k" '
    function gc_class(ref, alt) {
        if ((ref == "A" || ref == "T") && (alt == "G" || alt == "C")) return "WS"
        if ((ref == "G" || ref == "C") && (alt == "A" || alt == "T")) return "SW"
        return ""
    }
    function af_bin(af) {
        if (af < 0.0001)     return "<1e-4"
        else if (af < 0.001) return "1e-4_1e-3"
        else if (af < 0.01)  return "1e-3_0.01"
        else if (af < 0.05)  return "0.01_0.05"
        else if (af < 0.1)   return "0.05_0.1"
        else if (af < 0.5)   return "0.1_0.5"
        else                 return ">=0.5"
    }
    function af_order(b) {
        if (b == "<1e-4")     return 1
        if (b == "1e-4_1e-3") return 2
        if (b == "1e-3_0.01") return 3
        if (b == "0.01_0.05") return 4
        if (b == "0.05_0.1")  return 5
        if (b == "0.1_0.5")   return 6
        if (b == ">=0.5")     return 7
        return 9
    }
    {
        if ($9+0 > 0 && $10 != "." && $11 != ".") {
            cl = gc_class($10, $11)
            if (cl != "") {
                b = af_bin($9+0)
                key = cl "\t" b "\t" af_order(b)
                n[key]++
                total[cl]++
            }
        }
        if ($12+0 > 0 && $13 != "." && $14 != ".") {
            cl = gc_class($13, $14)
            if (cl != "") {
                b = af_bin($12+0)
                key = cl "\t" b "\t" af_order(b)
                n[key]++
                total[cl]++
            }
        }
    }
    END {
        for (key in n) {
            split(key, parts, "\t")
            dir = parts[1]
            pct = 100 * n[key] / total[dir]
            printf "%s\t%s\t%d\t%.4f\n", k, key, n[key], pct
        }
    }' >> "$OUTFILE4"
done

{ head -1 "$OUTFILE4"; tail -n +2 "$OUTFILE4" | sort -t$'\t' -k1,1n -k2,2 -k4,4n; } > "${OUTFILE4}.tmp"
mv "${OUTFILE4}.tmp" "$OUTFILE4"
echo "  Written: $OUTFILE4"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=== All gBGC summary tables written ==="
ls -lh "$OUTDIR"/gbgc_*.tsv
