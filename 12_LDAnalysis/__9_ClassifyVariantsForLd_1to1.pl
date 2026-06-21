#!/usr/bin/env perl
#===============================================================================
# LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 9
# Classify gnomAD variants using 1:1 GC positions
#
# Same logic as __1 but uses 1:1 positions from __8 instead of all GC positions.
# NonGC exclusion: excludes 1:1 GC positions from ALL k-values.
#
# Input:  LD_analysis_resampling_1to1/1to1_positions/{k}/{chr}.bed.gz
# Output: LD_analysis_resampling_1to1/{mapp}/k{k}_chr{chr}_{gc|nongc}_{maf}.bed
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
my $SCRIPTS_DIR = $FindBin::Bin;

my %mapp_beds = (
    "allmapp"      => "$cfg{paths}{regions_dir}/GRCh38_notinlowmappabilityall.bedsort.nochr.bed.gz",
    "nosegdupmapp" => "$cfg{paths}{regions_dir}/GRCh38_notinlowmappabilitynosegdups.bedsort.nochr.bed.gz",
    "segdupmapp"   => "$cfg{paths}{regions_dir}/GRCh38_notinlowmappabilitysegdups.bedsort.nochr.bed.gz",
);

my @k_values    = split /,/, $cfg{repeat_lengths}{values};
my @chromosomes = (1..22);

my @maf_bins = (
    { name => "0.001_0.005",   min => 0.001,  max => 0.005 },
    { name => "0.005_0.01",    min => 0.005,  max => 0.01 },
    { name => "0.01_0.05",     min => 0.01,   max => 0.05 },
    { name => "0.05_0.1",      min => 0.05,   max => 0.1 },
    { name => "0.1_0.5",       min => 0.1,    max => 0.5 },
    { name => "0.5_2",         min => 0.5,    max => 2 },
);

my $gnomad_dir     = $cfg{paths}{gnomad_genome_dir};
my $positions_base = "$BASE/LD_analysis_resampling_1to1/1to1_positions";
my $output_base    = "$BASE/LD_analysis_resampling_1to1";
my $jobs_dir       = "$output_base/jobs";

my $PERL_MOD = $cfg{slurm}{perl_module};
my $BED_MOD  = $cfg{modules}{bedtools};
my $HTS_MOD  = $cfg{modules}{htslib};

#===============================================================================
# Helper functions
#===============================================================================

sub get_higher_ks {
    my $k = shift;
    my @higher = ();
    for my $kv (@k_values) {
        push @higher, $kv if $kv > $k;
    }
    return @higher;
}

#===============================================================================
# Setup
#===============================================================================

make_path($jobs_dir) unless -d $jobs_dir;
for my $mapp (keys %mapp_beds) {
    make_path("$output_base/$mapp") unless -d "$output_base/$mapp";
}

#===============================================================================
# Generate jobs
#===============================================================================

print "="x70 . "\n";
print "LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 9: CLASSIFY VARIANTS\n";
print "Output: $output_base\n";
print "="x70 . "\n\n";

my $job_count = 0;

for my $mapp (sort keys %mapp_beds) {
    my $mapp_bed = $mapp_beds{$mapp};
    next unless -e $mapp_bed;

    print "=== $mapp ===\n";

    for my $k (@k_values) {
        # Check that 1:1 position files exist
        my $test_file = "$positions_base/$k/1.bed.gz";
        unless (-e $test_file) {
            print "  WARNING: $test_file not found, skipping k=$k\n";
            next;
        }

        for my $chr (@chromosomes) {
            $job_count++;
            my $job_file = "$jobs_dir/classify_${mapp}_k${k}_chr${chr}.sh";

            my $gc_file_chr = "$positions_base/$k/$chr.bed.gz";
            my @higher_ks = get_higher_ks($k);

            open(my $jf, ">", $job_file) or die "Cannot write $job_file: $!";
            print $jf <<"HEADER";
#!/bin/bash
#SBATCH --job-name=c1_${mapp}_k${k}_c${chr}
#SBATCH --output=$jobs_dir/classify_${mapp}_k${k}_chr${chr}.out
#SBATCH --error=$jobs_dir/classify_${mapp}_k${k}_chr${chr}.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module load $PERL_MOD
module load $BED_MOD
module load $HTS_MOD

echo "Classify 1:1: $mapp k=$k chr=$chr"

MAPP_BED="$mapp_bed"
GC_FILE="$gc_file_chr"
GNOMAD="$gnomad_dir/${chr}.sorted.bed.gz"
OUTPUT_DIR="$output_base/$mapp"
TMP="/tmp/c1_${mapp}_k${k}_chr${chr}_\$\$"
mkdir -p \$TMP

# Step 1: Extract gnomAD SNVs with MAF
zcat \$GNOMAD | awk -F'\\t' '
{
    split(\$4, a, "_"); split(a[4], b, "/")
    ref=b[1]; alt=b[2]
    if (ref ~ /^[ACGT]\$/ && alt ~ /^[ACGT]\$/) {
        print \$1"\\t"\$2"\\t"\$3"\\t"\$5
    }
}' > \$TMP/snvs.bed

# Step 2: Filter to mappable regions
zcat \$MAPP_BED | awk '\$1 == "$chr"' > \$TMP/mapp.bed
bedtools intersect -a \$TMP/snvs.bed -b \$TMP/mapp.bed -u -sorted > \$TMP/snvs_mapp.bed

# Step 3: Extract 1:1 GC regions for this k
zcat \$GC_FILE > \$TMP/gc_regions.bed

# Step 4: Classify
bedtools intersect -a \$TMP/snvs_mapp.bed -b \$TMP/gc_regions.bed -u -sorted > \$TMP/gc.bed
bedtools intersect -a \$TMP/snvs_mapp.bed -b \$TMP/gc_regions.bed -v -sorted > \$TMP/nongc_raw.bed

HEADER

            # Step 5: Exclude 1:1 GC positions from higher k-values from non-GC
            if (@higher_ks) {
                print $jf "# Step 5: Exclude 1:1 GC positions from higher k's from non-GC\n";
                print $jf "awk '{print \$1\"\\t\"\$2\"\\t\"\$3}' \$TMP/gc_regions.bed > \$TMP/all_gc_to_exclude.bed\n";
                for my $higher_k (@higher_ks) {
                    my $higher_file = "$positions_base/$higher_k/$chr.bed.gz";
                    print $jf "zcat $higher_file >> \$TMP/all_gc_to_exclude.bed\n";
                }
                print $jf "sort -k1,1 -k2,2n \$TMP/all_gc_to_exclude.bed > \$TMP/all_gc_sorted.bed\n";
                print $jf "bedtools intersect -a \$TMP/nongc_raw.bed -b \$TMP/all_gc_sorted.bed -v -sorted > \$TMP/nongc.bed\n";
            } else {
                print $jf "mv \$TMP/nongc_raw.bed \$TMP/nongc.bed\n";
            }

            print $jf <<"STATS";

echo "Total SNVs mappable: \$(wc -l < \$TMP/snvs_mapp.bed)"
echo "GC variants 1:1 (k=$k): \$(wc -l < \$TMP/gc.bed)"
echo "Non-GC variants: \$(wc -l < \$TMP/nongc.bed)"

# Step 6: Split by MAF bin
STATS

            for my $bin (@maf_bins) {
                my $name = $bin->{name};
                my $min = $bin->{min};
                my $max = $bin->{max};
                print $jf "awk '\$4 >= $min && \$4 < $max {print \$1\"\\t\"\$2\"\\t\"\$3}' \$TMP/gc.bed > \$OUTPUT_DIR/k${k}_chr${chr}_gc_${name}.bed\n";
                print $jf "awk '\$4 >= $min && \$4 < $max {print \$1\"\\t\"\$2\"\\t\"\$3}' \$TMP/nongc.bed > \$OUTPUT_DIR/k${k}_chr${chr}_nongc_${name}.bed\n";
            }

            print $jf <<"FOOTER";

rm -rf \$TMP
echo "Done: $mapp k=$k chr=$chr"
FOOTER
            close($jf);
            chmod(0755, $job_file);
        }
        print "  k=$k: 22 jobs\n";
    }
}

print "\n" . "="x70 . "\n";
print "TOTAAL: $job_count jobs in $jobs_dir\n";
print "="x70 . "\n\n";
print "Submit:\n  for f in $jobs_dir/classify_*.sh; do sbatch -M shinx \$f; done\n\n";
