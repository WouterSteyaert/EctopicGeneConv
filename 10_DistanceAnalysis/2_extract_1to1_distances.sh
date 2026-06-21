#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=24:00:00

# Extract coordinate pairs and distances from 1:1 filtered gene conversion data.
# For each 1:1 entry, extracts the two genomic coordinates and computes the
# distance for intra-chromosomal pairs.
#
# Input:  distance_analysis/{REPLEN}/*.1to1.txt.gz  (from filter_1to1_coords.sh)
# Output: distance_analysis/{REPLEN}/distances.txt.gz
#
# Each pair appears twice in the input (donor->acceptor and reverse).
# We deduplicate by always putting the smaller coordinate first.

: "${PROJECT_ROOT:?PROJECT_ROOT must be set (see 00_Configuration/README.md)}"
BASEDIR="${PROJECT_ROOT}/geneconv_complete"

# Parse --REPLEN= argument
for arg in "$@"; do
    case $arg in
        --REPLEN=*) REPLEN="${arg#*=}" ;;
    esac
done

if [ -z "$REPLEN" ]; then
    echo "ERROR: --REPLEN=<value> required"
    exit 1
fi

INPATH="${BASEDIR}/distance_analysis/${REPLEN}"
OUTFILE="${INPATH}/distances.txt.gz"

if [ ! -d "$INPATH" ]; then
    echo "ERROR: ${INPATH} does not exist"
    exit 1
fi

echo "Processing repeat length ${REPLEN}..."

# Process all 1:1 files:
# 1. Extract the two coordinates from each line
# 2. Deduplicate pairs (each appears twice with swapped columns)
# 3. Compute distance for intra-chromosomal, flag inter-chromosomal
#
# Output columns: chr1  pos1  base1  chr2  pos2  base2  distance  type

zcat "${INPATH}"/*.1to1.txt.gz | awk '
BEGIN {
    FS = "\t"
    OFS = "\t"
}
{
    # Extract coord from col2 or col3 (whichever is not "/")
    if ($2 != "/") { coord1 = $2 } else { coord1 = $3 }
    # Extract coord from col4 or col5 (whichever is not "/")
    if ($4 != "/") { coord2 = $4 } else { coord2 = $5 }

    # Parse Base=Chr:Pos
    split(coord1, c1, /[=:]/)
    base1 = c1[1]; chr1 = c1[2]; pos1 = c1[3] + 0

    split(coord2, c2, /[=:]/)
    base2 = c2[1]; chr2 = c2[2]; pos2 = c2[3] + 0

    # Deduplicate: canonical order (smaller chr first, or smaller pos if same chr)
    if (chr1 > chr2 || (chr1 == chr2 && pos1 > pos2)) {
        # swap
        tmp = chr1; chr1 = chr2; chr2 = tmp
        tmp = pos1; pos1 = pos2; pos2 = tmp
        tmp = base1; base1 = base2; base2 = tmp
    }

    key = chr1 ":" pos1 ":" chr2 ":" pos2
    if (seen[key]) next
    seen[key] = 1

    # Compute distance
    if (chr1 == chr2) {
        dist = pos2 - pos1
        type = "intra"
    } else {
        dist = "NA"
        type = "inter"
    }

    print chr1, pos1, base1, chr2, pos2, base2, dist, type
}
' | sort -k1,1V -k2,2n | gzip > "$OUTFILE"

# Summary statistics
total=$(zcat "$OUTFILE" | wc -l)
intra=$(zcat "$OUTFILE" | awk '$8=="intra"' | wc -l)
inter=$(zcat "$OUTFILE" | awk '$8=="inter"' | wc -l)

echo "RepLen ${REPLEN}: ${total} unique pairs (${intra} intra-chromosomal, ${inter} inter-chromosomal)"

# Distance distribution for intra-chromosomal (log10 bins)
echo ""
echo "Intra-chromosomal distance distribution:"
zcat "$OUTFILE" | awk '$8=="intra" {
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
    order[1]="<1kb"; order[2]="1-10kb"; order[3]="10-100kb"
    order[4]="100kb-1Mb"; order[5]="1-10Mb"; order[6]="10-100Mb"; order[7]=">100Mb"
    for (i=1; i<=7; i++) {
        if (counts[order[i]]) printf "  %-12s %d\n", order[i], counts[order[i]]
    }
}'
