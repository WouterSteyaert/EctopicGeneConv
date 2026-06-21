#!/usr/bin/env perl
#===============================================================================
# LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 8
# Extract 1:1 paralog GC positions from DistanceAnalysis output
#
# Input:  distance_analysis/{k}/distances.diff.txt.gz
#         (already k-specific, already filtered for strict 1:1 pairs)
#
# Output: LD_analysis_resampling_1to1/1to1_positions/{k}/{chr}.bed.gz
#         (0-based BED, sorted, deduplicated, chr1-22 only)
#
# Generates one SLURM job per k-value.
#===============================================================================

use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

#===============================================================================
# Configuration
#===============================================================================

my %cfg = load_config();

my $BASE        = $cfg{paths}{data_root};
my $DIST_BASE   = "$BASE/distance_analysis";
my $OUTPUT_BASE = "$BASE/LD_analysis_resampling_1to1";
my $JOBS_DIR    = "$OUTPUT_BASE/jobs";

my @k_values    = split /,/, $cfg{repeat_lengths}{values};
my @chromosomes = (1..22);

#===============================================================================
# Setup
#===============================================================================

make_path($JOBS_DIR) unless -d $JOBS_DIR;
for my $k (@k_values) {
    make_path("$OUTPUT_BASE/1to1_positions/$k") unless -d "$OUTPUT_BASE/1to1_positions/$k";
}

#===============================================================================
# Generate jobs
#===============================================================================

print "="x70 . "\n";
print "LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 8: EXTRACT POSITIONS\n";
print "="x70 . "\n\n";

my $job_count = 0;

for my $k (@k_values) {
    my $input = "$DIST_BASE/$k/distances.diff.txt.gz";
    unless (-f $input) {
        print "WARNING: $input not found, skipping k=$k\n";
        next;
    }

    $job_count++;
    my $job_file = "$JOBS_DIR/extract_k${k}.sh";
    my $out_dir = "$OUTPUT_BASE/1to1_positions/$k";

    # Estimate resources: k=17/19 are 100M+ lines
    my $mem = ($k <= 19) ? "16G" : "4G";
    my $time = ($k <= 19) ? "2:00:00" : "0:30:00";

    open(my $jf, ">", $job_file) or die "Cannot write $job_file: $!";
    print $jf <<"EOF";
#!/bin/bash
#SBATCH --job-name=ex1to1_k${k}
#SBATCH --output=$JOBS_DIR/extract_k${k}.out
#SBATCH --error=$JOBS_DIR/extract_k${k}.err
#SBATCH --time=$time
#SBATCH --mem=$mem
#SBATCH --cpus-per-task=1

echo "========================================"
echo "Extract 1:1 positions: k=$k"
echo "========================================"

INPUT="$input"
OUT_DIR="$out_dir"
TMP="/tmp/ex1to1_k${k}_\$\$"
mkdir -p \$TMP

echo "Reading \$INPUT ..."
N_PAIRS=\$(zcat \$INPUT | wc -l)
echo "Total pairs: \$N_PAIRS"

# Extract both pos1 (chr1:pos1) and pos2 (chr2:pos2) into per-chromosome files
# Convert to 0-based BED: chr  pos-1  pos
# Only keep chr 1-22
zcat \$INPUT | awk -F'\\t' '{
    # pos1 side
    if (\$1 >= 1 && \$1 <= 22) {
        print \$1"\\t"(\$2-1)"\\t"\$2 >> "'\$TMP'/chr_"\$1".bed"
    }
    # pos2 side
    if (\$4 >= 1 && \$4 <= 22) {
        print \$4"\\t"(\$5-1)"\\t"\$5 >> "'\$TMP'/chr_"\$4".bed"
    }
}'

echo "Sorting and deduplicating per chromosome..."
for chr in {1..22}; do
    if [[ -f \$TMP/chr_\${chr}.bed ]]; then
        sort -k1,1 -k2,2n \$TMP/chr_\${chr}.bed | uniq | gzip > \$OUT_DIR/\${chr}.bed.gz
        N=\$(zcat \$OUT_DIR/\${chr}.bed.gz | wc -l)
        echo "  chr\$chr: \$N unique positions"
    else
        # Empty file for missing chromosomes
        echo -n | gzip > \$OUT_DIR/\${chr}.bed.gz
        echo "  chr\$chr: 0 positions"
    fi
done

TOTAL=0
for chr in {1..22}; do
    N=\$(zcat \$OUT_DIR/\${chr}.bed.gz | wc -l)
    TOTAL=\$((TOTAL + N))
done
echo ""
echo "Total unique 1:1 positions (chr1-22): \$TOTAL"

rm -rf \$TMP
echo "Done!"
EOF
    close($jf);
    chmod 0755, $job_file;
    print "  k=$k: $job_file\n";
}

print "\n" . "="x70 . "\n";
print "TOTAAL: $job_count jobs in $JOBS_DIR\n";
print "="x70 . "\n";
print "\nSubmit:\n  for f in $JOBS_DIR/extract_*.sh; do sbatch -M shinx \$f; done\n\n";
