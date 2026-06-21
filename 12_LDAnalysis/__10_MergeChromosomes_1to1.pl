#!/usr/bin/env perl
#===============================================================================
# LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 10
# Merge chromosome files per k/maf combination
#
# Output: LD_analysis_resampling_1to1/
#===============================================================================

use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $BASE        = $cfg{paths}{data_root};
my $output_base = "$BASE/LD_analysis_resampling_1to1";
my $TMP_DIR     = $ENV{TMPDIR} // "$BASE/temp";
my $PARALLEL_CORES = $ENV{SLURM_CPUS_PER_TASK} // 4;

my @mapp_categories = ("allmapp", "nosegdupmapp", "segdupmapp");
my @k_values        = split /,/, $cfg{repeat_lengths}{values};
my @maf_bins = ("0.001_0.005", "0.005_0.01", "0.01_0.05", "0.05_0.1", "0.1_0.5", "0.5_2");

# Ensure tmp dir exists
make_path($TMP_DIR) unless -d $TMP_DIR;

print "="x70 . "\n";
print "LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 10: MERGE CHROMOSOMES\n";
print "Using $PARALLEL_CORES cores for sorting, tmp: $TMP_DIR\n";
print "="x70 . "\n\n";

for my $mapp (@mapp_categories) {
    print "=== $mapp ===\n";
    my $dir = "$output_base/$mapp";

    for my $k (@k_values) {
        for my $bin (@maf_bins) {
            my $gc_out = "$dir/k${k}_gc_${bin}.bed";
            my $nongc_out = "$dir/k${k}_nongc_${bin}.bed";

            # Merge GC files
            my @gc_files = glob("$dir/k${k}_chr*_gc_${bin}.bed");
            if (@gc_files) {
                system("cat @gc_files | sort --parallel=$PARALLEL_CORES -T $TMP_DIR -k1,1n -k2,2n > $gc_out");
            }

            # Merge non-GC files
            my @nongc_files = glob("$dir/k${k}_chr*_nongc_${bin}.bed");
            if (@nongc_files) {
                system("cat @nongc_files | sort --parallel=$PARALLEL_CORES -T $TMP_DIR -k1,1n -k2,2n > $nongc_out");
            }
        }
        print "  k=$k merged\n";
    }
}

print "\nDone! Merged files in $output_base/*/\n";
