#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=24:00:00

# Filter gene conversion coords files for 1:1 relationships
# for a single repeat length (passed via --REPLEN= or REPLEN env var).
#
# 1:1 = exactly 1 coordinate in columns 2+3 combined,
#        AND exactly 1 coordinate in columns 4+5 combined.
# Each conversion event appears twice (donor->acceptor and acceptor->donor).

: "${PROJECT_ROOT:?PROJECT_ROOT must be set (see 00_Configuration/README.md)}"
BASEDIR="${PROJECT_ROOT}/geneconv_complete"
INDIR="${BASEDIR}/conv"
OUTDIR="${BASEDIR}/distance_analysis"

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

INPATH="${INDIR}/${REPLEN}"
OUTPATH="${OUTDIR}/${REPLEN}"

if [ ! -d "$INPATH" ]; then
    echo "ERROR: ${INPATH} does not exist"
    exit 1
fi

mkdir -p "$OUTPATH"

# AWK filter for 1:1 entries
AWK_FILTER='
BEGIN { FS="\t" }
{
    n23 = 0
    if ($2 != "/") n23 += split($2, a, ",")
    if ($3 != "/") n23 += split($3, a, ",")
    n45 = 0
    if ($4 != "/") n45 += split($4, a, ",")
    if ($5 != "/") n45 += split($5, a, ",")
    if (n23 == 1 && n45 == 1) print
}
'

total=0
kept=0
files=0

for f in "${INPATH}"/*.repsum.geneconv.coords.txt.gz; do
    [ -f "$f" ] || continue
    base=$(basename "$f" .txt.gz)
    outfile="${OUTPATH}/${base}.1to1.txt.gz"

    before=$(zcat "$f" | wc -l)
    zcat "$f" | awk "$AWK_FILTER" | gzip > "$outfile"
    after=$(zcat "$outfile" | wc -l)

    total=$((total + before))
    kept=$((kept + after))
    files=$((files + 1))

    # Remove empty output files
    if [ "$after" -eq 0 ]; then
        rm "$outfile"
    fi
done

if [ "$total" -gt 0 ]; then
    pct=$(echo "scale=1; ${kept}*100/${total}" | bc)
else
    pct="0"
fi

echo "RepLen ${REPLEN}: ${files} files, ${kept}/${total} lines kept (${pct}%)"
