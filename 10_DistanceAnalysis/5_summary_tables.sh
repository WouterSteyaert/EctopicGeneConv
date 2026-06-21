#!/bin/bash
# =============================================================================
# Step 5: Generate summary tables from annotated diff pairs
# =============================================================================
#
# Reads distances.diff.gnomad.txt.gz for all k values and produces
# summary tables for downstream analysis and plotting.
#
# Input:  distance_analysis/{k}/distances.diff.gnomad.txt.gz
# Output: distance_analysis/summary/
#           - distance_af_by_k.tsv
#           - transitions_by_k.tsv
#           - transitions_by_k_distance.tsv
#           - complementary_pairs.tsv
#           - symmetry_by_k.tsv
#
# Columns in input:
#   1:chr1  2:pos1  3:base1  4:chr2  5:pos2  6:base2
#   7:distance  8:type  9:AF1  10:ref1  11:alt1  12:AF2  13:ref2  14:alt2
# =============================================================================

: "${PROJECT_ROOT:?PROJECT_ROOT must be set (see 00_Configuration/README.md)}"
BASEDIR="${PROJECT_ROOT}/geneconv_complete"
OUTDIR="${BASEDIR}/distance_analysis/summary"
mkdir -p "$OUTDIR"

K_VALUES="17 19 21 31 41 51 61 71 81 91"

# =============================================================================
# Table 1: Distance vs AF by k
# =============================================================================
echo "=== Table 1: Distance vs AF by k ==="

OUTFILE1="${OUTDIR}/distance_af_by_k.tsv"
echo -e "k\tdist_bin\tdist_order\tn_total\tn_intra\tn_inter\tn_any_var\tn_both_var\tmean_AF\tpct_WS\tpct_SW" \
    > "$OUTFILE1"

for k in $K_VALUES; do
    INFILE="${BASEDIR}/distance_analysis/${k}/distances.diff.gnomad.txt.gz"
    [ ! -f "$INFILE" ] && echo "  WARNING: $INFILE missing, skipping" && continue
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
        if (b == "<1kb")       return 1
        if (b == "1-10kb")     return 2
        if (b == "10-100kb")   return 3
        if (b == "100kb-1Mb")  return 4
        if (b == "1-10Mb")     return 5
        if (b == "10-100Mb")   return 6
        if (b == ">100Mb")     return 7
        if (b == "inter")      return 8
        return 9
    }
    function gc_class(ref, alt) {
        if ((ref == "A" || ref == "T") && (alt == "G" || alt == "C")) return "WS"
        if ((ref == "G" || ref == "C") && (alt == "A" || alt == "T")) return "SW"
        return "other"
    }
    {
        b = ($8 == "inter") ? "inter" : dist_bin($7)
        n_total[b]++
        if ($8 == "intra") n_intra[b]++
        else               n_inter[b]++

        has_var1 = ($9+0 > 0 && $10 != ".")
        has_var2 = ($12+0 > 0 && $13 != ".")

        if (has_var1 || has_var2) n_any_var[b]++
        if (has_var1 && has_var2) n_both_var[b]++

        # Mean AF: use max(AF1, AF2) per pair when at least one has variant
        if (has_var1 || has_var2) {
            maxaf = ($9+0 > $12+0) ? $9+0 : $12+0
            sum_af[b] += maxaf
            n_af[b]++
        }

        # GC-bias counts
        if (has_var1) {
            cl = gc_class($10, $11)
            if (cl == "WS") n_ws[b]++
            if (cl == "SW") n_sw[b]++
            n_gc_total[b]++
        }
        if (has_var2) {
            cl = gc_class($13, $14)
            if (cl == "WS") n_ws[b]++
            if (cl == "SW") n_sw[b]++
            n_gc_total[b]++
        }
    }
    END {
        for (b in n_total) {
            mean_af = (n_af[b] > 0) ? sum_af[b] / n_af[b] : 0
            pct_ws = (n_gc_total[b] > 0) ? 100 * n_ws[b] / n_gc_total[b] : 0
            pct_sw = (n_gc_total[b] > 0) ? 100 * n_sw[b] / n_gc_total[b] : 0
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%.1f\t%.1f\n",
                k, b, dist_order(b), n_total[b], n_intra[b]+0, n_inter[b]+0,
                n_any_var[b]+0, n_both_var[b]+0, mean_af, pct_ws, pct_sw
        }
    }' >> "$OUTFILE1"
done

{ head -1 "$OUTFILE1"; tail -n +2 "$OUTFILE1" | sort -t$'\t' -k1,1n -k3,3n; } > "${OUTFILE1}.tmp"
mv "${OUTFILE1}.tmp" "$OUTFILE1"
echo "  Written: $OUTFILE1"

# =============================================================================
# Table 2: Specific transitions by k
# =============================================================================
echo ""
echo "=== Table 2: Transitions by k ==="

OUTFILE2="${OUTDIR}/transitions_by_k.tsv"
echo -e "k\ttransition\tn\tmean_AF\tdirection" > "$OUTFILE2"

for k in $K_VALUES; do
    INFILE="${BASEDIR}/distance_analysis/${k}/distances.diff.gnomad.txt.gz"
    [ ! -f "$INFILE" ] && continue
    echo "  Processing k=$k..."

    zcat "$INFILE" | awk -F'\t' -v k="$k" '
    function gc_dir(ref, alt) {
        if ((ref == "A" || ref == "T") && (alt == "G" || alt == "C")) return "WS"
        if ((ref == "G" || ref == "C") && (alt == "A" || alt == "T")) return "SW"
        if ((ref == "A" || ref == "T") && (alt == "A" || alt == "T")) return "WW"
        if ((ref == "G" || ref == "C") && (alt == "G" || alt == "C")) return "SS"
        return "other"
    }
    {
        # Position 1 variant
        if ($9+0 > 0 && $10 != "." && $11 != ".") {
            tr = $10 ">" $11
            n[tr]++
            sum_af[tr] += $9+0
            dir[tr] = gc_dir($10, $11)
        }
        # Position 2 variant
        if ($12+0 > 0 && $13 != "." && $14 != ".") {
            tr = $13 ">" $14
            n[tr]++
            sum_af[tr] += $12+0
            dir[tr] = gc_dir($13, $14)
        }
    }
    END {
        for (tr in n) {
            mean = sum_af[tr] / n[tr]
            printf "%s\t%s\t%d\t%.6f\t%s\n", k, tr, n[tr], mean, dir[tr]
        }
    }' >> "$OUTFILE2"
done

echo "  Written: $OUTFILE2"

# =============================================================================
# Table 3: Transitions by k x distance
# =============================================================================
echo ""
echo "=== Table 3: Transitions by k x distance ==="

OUTFILE3="${OUTDIR}/transitions_by_k_distance.tsv"
echo -e "k\tdist_bin\tdist_order\ttransition\tn\tmean_AF\tdirection" > "$OUTFILE3"

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
    function gc_dir(ref, alt) {
        if ((ref == "A" || ref == "T") && (alt == "G" || alt == "C")) return "WS"
        if ((ref == "G" || ref == "C") && (alt == "A" || alt == "T")) return "SW"
        if ((ref == "A" || ref == "T") && (alt == "A" || alt == "T")) return "WW"
        if ((ref == "G" || ref == "C") && (alt == "G" || alt == "C")) return "SS"
        return "other"
    }
    {
        b = ($8 == "inter") ? "inter" : dist_bin($7)
        bo = dist_order(b)

        if ($9+0 > 0 && $10 != "." && $11 != ".") {
            tr = $10 ">" $11
            key = b "\t" bo "\t" tr
            n[key]++
            sum_af[key] += $9+0
            dir[key] = gc_dir($10, $11)
        }
        if ($12+0 > 0 && $13 != "." && $14 != ".") {
            tr = $13 ">" $14
            key = b "\t" bo "\t" tr
            n[key]++
            sum_af[key] += $12+0
            dir[key] = gc_dir($13, $14)
        }
    }
    END {
        for (key in n) {
            mean = sum_af[key] / n[key]
            printf "%s\t%s\t%d\t%.6f\t%s\n", k, key, n[key], mean, dir[key]
        }
    }' >> "$OUTFILE3"
done

echo "  Written: $OUTFILE3"

# =============================================================================
# Table 4: Complementary transition pairs
# =============================================================================
echo ""
echo "=== Table 4: Complementary transition pairs ==="

OUTFILE4="${OUTDIR}/complementary_pairs.tsv"
echo -e "k\tpair\tforward\tn_fwd\tmeanAF_fwd\treverse\tn_rev\tmeanAF_rev\tratio_fwd_rev" > "$OUTFILE4"

for k in $K_VALUES; do
    INFILE="${BASEDIR}/distance_analysis/${k}/distances.diff.gnomad.txt.gz"
    [ ! -f "$INFILE" ] && continue
    echo "  Processing k=$k..."

    zcat "$INFILE" | awk -F'\t' -v k="$k" '
    {
        if ($9+0 > 0 && $10 != "." && $11 != ".") {
            tr = $10 ">" $11
            n[tr]++
            sum_af[tr] += $9+0
        }
        if ($12+0 > 0 && $13 != "." && $14 != ".") {
            tr = $13 ">" $14
            n[tr]++
            sum_af[tr] += $12+0
        }
    }
    END {
        # Define complementary pairs
        split("A>G,G>A,A>C,C>A,A>T,T>A,C>T,T>C,C>G,G>C,G>T,T>G", pairs, ",")
        for (i = 1; i <= 12; i += 2) {
            fwd = pairs[i]
            rev = pairs[i+1]
            n_f = n[fwd]+0; n_r = n[rev]+0
            af_f = (n_f > 0) ? sum_af[fwd] / n_f : 0
            af_r = (n_r > 0) ? sum_af[rev] / n_r : 0
            ratio = (n_r > 0) ? n_f / n_r : 0
            printf "%s\t%s/%s\t%s\t%d\t%.6f\t%s\t%d\t%.6f\t%.3f\n",
                k, fwd, rev, fwd, n_f, af_f, rev, n_r, af_r, ratio
        }
    }' >> "$OUTFILE4"
done

echo "  Written: $OUTFILE4"

# =============================================================================
# Table 5: Symmetry analysis (AF1 vs AF2) by k and distance
# =============================================================================
echo ""
echo "=== Table 5: Symmetry by k ==="

OUTFILE5="${OUTDIR}/symmetry_by_k.tsv"
echo -e "k\tdist_bin\tdist_order\tn_both\tmean_AF1\tmean_AF2\tmean_abs_diff\tmean_min_max_ratio" \
    > "$OUTFILE5"

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
    $9+0 > 0 && $12+0 > 0 {
        b = ($8 == "inter") ? "inter" : dist_bin($7)
        n[b]++
        sum_af1[b] += $9+0
        sum_af2[b] += $12+0
        diff = $9 - $12; if (diff < 0) diff = -diff
        sum_diff[b] += diff
        if ($9+0 > $12+0)
            sum_ratio[b] += $12 / $9
        else
            sum_ratio[b] += $9 / $12
    }
    END {
        for (b in n) {
            printf "%s\t%s\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.3f\n",
                k, b, dist_order(b), n[b],
                sum_af1[b]/n[b], sum_af2[b]/n[b],
                sum_diff[b]/n[b], sum_ratio[b]/n[b]
        }
    }' >> "$OUTFILE5"
done

{ head -1 "$OUTFILE5"; tail -n +2 "$OUTFILE5" | sort -t$'\t' -k1,1n -k3,3n; } > "${OUTFILE5}.tmp"
mv "${OUTFILE5}.tmp" "$OUTFILE5"
echo "  Written: $OUTFILE5"

echo ""
echo "=== All summary tables written to: ${OUTDIR}/ ==="
ls -lh "$OUTDIR"/*.tsv
