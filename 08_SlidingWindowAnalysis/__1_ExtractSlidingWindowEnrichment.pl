#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use POSIX qw(log10);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

# ============================================================================
# Sliding Window Enrichment Analysis (uncorrected; diagnostic only).
#
# Reads per-window contingency counts from step 06 and emits per-window
# enrichment (gc_rate, nongc_rate, odds_ratio, log2_or) plus region annotation
# (pericentromeric / subtelomeric / interstitial) using GRCh38 cytoBand acen
# bands.
#
# The CORRECTED sliding-window analysis used in the paper is in
# 1_correct_enrichment.pl (control group is restricted to positions not GC-
# compatible at ANY k).  This script is retained for diagnostic comparison.
# ============================================================================

my %cfg = load_config();

# --- Configuration ---
my $BASE_STATS = $cfg{paths}{stats_dir};
my $OUTPUT_DIR = "$cfg{paths}{data_root}/sliding_window_analysis";

# Parameters with defaults
my $window_size = 1000000;   # 1 Mb
my $step_size   = 50000;     # 50 kb
my $source      = "gnomad~genome";
my $stratum     = "allmapp";
my @k_values    = (21, 31, 41, 51, 61, 71, 81, 91);
my @maf_bins    = ("01_05", "005_01", "001_005");  # Focus on MAF >= 0.01 by default
my $help        = 0;

GetOptions(
    "window=i"    => \$window_size,
    "step=i"      => \$step_size,
    "source=s"    => \$source,
    "stratum=s"   => \$stratum,
    "k=s"         => sub { @k_values = split(/,/, $_[1]) },
    "maf=s"       => sub { @maf_bins = split(/,/, $_[1]) },
    "output=s"    => \$OUTPUT_DIR,
    "help"        => \$help,
) or die "Error in command line arguments\n";

if ($help) {
    print_usage();
    exit 0;
}

# --- GRCh38 Centromere boundaries (acen bands from UCSC cytoBand) ---
# p-arm edge (start of centromeric region)
my %acen_start = (
    1  => 121700000, 2  => 91800000,  3  => 87800000,  4  => 48200000,
    5  => 46100000,  6  => 58500000,  7  => 58100000,  8  => 43200000,
    9  => 42200000,  10 => 38000000,  11 => 51000000,  12 => 33200000,
    13 => 16500000,  14 => 16100000,  15 => 17500000,  16 => 35300000,
    17 => 22700000,  18 => 15400000,  19 => 24200000,  20 => 25700000,
    21 => 10900000,  22 => 13700000,
    X  => 58100000,  Y  => 10300000,
);
# q-arm edge (end of centromeric region)
my %acen_end = (
    1  => 125100000, 2  => 96000000,  3  => 94000000,  4  => 51800000,
    5  => 51400000,  6  => 62600000,  7  => 62100000,  8  => 47200000,
    9  => 45500000,  10 => 41600000,  11 => 55800000,  12 => 37800000,
    13 => 18900000,  14 => 18200000,  15 => 20500000,  16 => 38400000,
    17 => 27400000,  18 => 21500000,  19 => 28100000,  20 => 30400000,
    21 => 13000000,  22 => 17400000,
    X  => 63800000,  Y  => 10600000,
);
# Midpoints for centromere_dist column and plotting
my %centromeres = map { $_ => int(($acen_start{$_} + $acen_end{$_}) / 2) } keys %acen_start;

# --- GRCh38 Chromosome sizes ---
my %chrom_sizes = (
    1  => 248956422,
    2  => 242193529,
    3  => 198295559,
    4  => 190214555,
    5  => 181538259,
    6  => 170805979,
    7  => 159345973,
    8  => 145138636,
    9  => 138394717,
    10 => 133797422,
    11 => 135086622,
    12 => 133275309,
    13 => 114364328,
    14 => 107043718,
    15 => 101991189,
    16 => 90338345,
    17 => 83257441,
    18 => 80373285,
    19 => 58617616,
    20 => 64444167,
    21 => 46709983,
    22 => 50818468,
    X  => 156040895,
    Y  => 57227415,
);

# --- Create output directory ---
system("mkdir -p $OUTPUT_DIR") == 0 or die "Cannot create output directory: $OUTPUT_DIR\n";

# --- Main processing ---
print STDERR "Sliding Window Enrichment Analysis\n";
print STDERR "===================================\n";
print STDERR "Stratum: $stratum\n";
print STDERR "K values: " . join(", ", @k_values) . "\n";
print STDERR "MAF bins: " . join(", ", @maf_bins) . "\n";
print STDERR "Window: ${window_size}bp, Step: ${step_size}bp\n\n";

my $stats_dir = "$BASE_STATS/stats_$stratum";

# Process each k value
foreach my $k (@k_values) {
    print STDERR "Processing k=$k...\n";

    my $k_dir = "$stats_dir/$k";
    unless (-d $k_dir) {
        print STDERR "  WARNING: Directory not found: $k_dir\n";
        next;
    }

    # Process each MAF bin
    foreach my $maf (@maf_bins) {
        print STDERR "  MAF bin: $maf\n";

        my $output_file = "$OUTPUT_DIR/enrichment_${stratum}_k${k}_maf${maf}.tsv";
        open(my $out_fh, ">", $output_file) or die "Cannot open output: $output_file\n";

        # Header
        print $out_fh join("\t", qw(
            chr window_start window_end window_mid
            gc_var gc_novar nongc_var nongc_novar
            gc_rate nongc_rate odds_ratio log2_or
            centromere_dist telomere_dist region_type
        )) . "\n";

        my $total_windows = 0;

        # Process each chromosome
        foreach my $chr (1..22, 'X', 'Y') {
            my $chrom_size = $chrom_sizes{$chr};
            my $pattern = "${chr}_${window_size}_${step_size}_${source}_${maf}.*.txt";

            # Find the file (exclude annotation-specific files)
            my @files = glob("$k_dir/${chr}_${window_size}_${step_size}_${source}_${maf}.1.*.txt");
            @files = grep { $_ !~ /\.(Exons|Genes|RecombHotspots|SegmentalDup|cpgIsland|Enhancer|CorePromoter|OpenChromatine|TfBinding|CtcfBinding|EarlyReplicating|LateReplicating|cpgIslandUnmasked)\.txt$/ } @files;

            unless (@files) {
                print STDERR "    Chr $chr: No file found\n";
                next;
            }

            my $file = $files[0];

            open(my $fh, "<", $file) or do {
                print STDERR "    Chr $chr: Cannot open $file\n";
                next;
            };

            my $header = <$fh>;  # Skip header
            my $window_count = 0;

            while (<$fh>) {
                chomp;
                my @cols = split(/\t/);

                # Extract columns
                my $chrom        = $cols[0];
                my $win_start    = $cols[1];
                my $win_end      = $cols[2];
                my $gc_var       = $cols[3];   # NrOfVarPosConvPos
                my $nongc_var    = $cols[4];   # NrOfVarPosNoConvPos
                my $gc_novar     = $cols[5];   # NrOfNoVarPosConvPos
                my $nongc_novar  = $cols[6];   # NrOfNoVarPosNoConvPos

                # Skip windows with no data
                next if ($gc_var + $gc_novar == 0 || $nongc_var + $nongc_novar == 0);

                # Calculate rates
                my $gc_rate    = $gc_var / ($gc_var + $gc_novar);
                my $nongc_rate = $nongc_var / ($nongc_var + $nongc_novar);

                # Calculate odds ratio (with pseudocount to avoid division by zero)
                my $pseudo = 0.5;
                my $or = (($gc_var + $pseudo) * ($nongc_novar + $pseudo)) /
                         (($gc_novar + $pseudo) * ($nongc_var + $pseudo));
                my $log2_or = log($or) / log(2);

                # Calculate distances
                my $win_mid = ($win_start + $win_end) / 2;
                # Distance to nearest centromere EDGE (not midpoint)
                my $centro_dist;
                if ($win_mid >= $acen_start{$chr} && $win_mid <= $acen_end{$chr}) {
                    $centro_dist = 0;  # Inside centromere
                } else {
                    my $dist_p = abs($win_mid - $acen_start{$chr});
                    my $dist_q = abs($win_mid - $acen_end{$chr});
                    $centro_dist = ($dist_p < $dist_q) ? $dist_p : $dist_q;
                }
                my $telo_dist = min($win_mid, $chrom_size - $win_mid);

                # Classify region
                my $region_type = classify_region($win_mid, $chr, $chrom_size);

                # Output
                print $out_fh join("\t",
                    $chr, $win_start, $win_end, int($win_mid),
                    $gc_var, $gc_novar, $nongc_var, $nongc_novar,
                    sprintf("%.6f", $gc_rate),
                    sprintf("%.6f", $nongc_rate),
                    sprintf("%.4f", $or),
                    sprintf("%.4f", $log2_or),
                    int($centro_dist),
                    int($telo_dist),
                    $region_type
                ) . "\n";

                $window_count++;
            }
            close($fh);
            $total_windows += $window_count;
        }

        close($out_fh);
        print STDERR "    Output: $output_file ($total_windows windows)\n";
    }
}

# --- Generate combined file for all MAF bins ---
print STDERR "\nGenerating combined output...\n";
my $combined_file = "$OUTPUT_DIR/enrichment_${stratum}_combined.tsv";
open(my $comb_fh, ">", $combined_file) or die "Cannot open: $combined_file\n";
print $comb_fh join("\t", qw(k maf_bin chr window_start window_end window_mid gc_var gc_novar nongc_var nongc_novar gc_rate nongc_rate odds_ratio log2_or centromere_dist telomere_dist region_type)) . "\n";

foreach my $k (@k_values) {
    foreach my $maf (@maf_bins) {
        my $file = "$OUTPUT_DIR/enrichment_${stratum}_k${k}_maf${maf}.tsv";
        next unless -f $file;

        open(my $fh, "<", $file) or next;
        <$fh>;  # Skip header
        while (<$fh>) {
            chomp;
            print $comb_fh "$k\t$maf\t$_\n";
        }
        close($fh);
    }
}
close($comb_fh);
print STDERR "Combined output: $combined_file\n";

print STDERR "\nDone!\n";

# ============================================================================
# Subroutines
# ============================================================================

sub min {
    my ($a, $b) = @_;
    return $a < $b ? $a : $b;
}

sub classify_region {
    my ($pos, $chr, $chrom_size) = @_;

    # Distance from nearest centromere EDGE (not midpoint)
    my $dist_to_p_edge = abs($pos - $acen_start{$chr});
    my $dist_to_q_edge = abs($pos - $acen_end{$chr});
    my $edge_dist = ($dist_to_p_edge < $dist_to_q_edge) ? $dist_to_p_edge : $dist_to_q_edge;

    # Inside centromere region itself = pericentromeric
    if ($pos >= $acen_start{$chr} && $pos <= $acen_end{$chr}) {
        return "pericentromeric";
    }

    # 5 Mb from centromere edge = pericentromeric, 5 Mb from telomere = subtelomeric
    # (Methods §Sliding window analysis; matches the thresholds used in
    # 1_correct_enrichment.pl and 2_identify_peaks.pl.)
    my $pericen_threshold = 5_000_000;
    my $subtelo_threshold = 5_000_000;

    if ($edge_dist < $pericen_threshold) {
        return "pericentromeric";
    } elsif ($pos < $subtelo_threshold || $pos > $chrom_size - $subtelo_threshold) {
        return "subtelomeric";
    } else {
        return "interstitial";
    }
}

sub print_usage {
    print <<EOF;
Usage: $0 [options]

Options:
    --stratum STR    Mappability stratum (allmapp, nosegdupmapp, segdupmapp)
                     Default: allmapp
    --k LIST         Comma-separated list of k values
                     Default: 21,31,41,51,61,71,81,91
    --maf LIST       Comma-separated list of MAF bins
                     Default: 01_05,005_01,001_005
    --window INT     Window size in bp (default: 1000000)
    --step INT       Step size in bp (default: 50000)
    --source STR     Data source (default: gnomad~genome)
    --output DIR     Output directory
    --help           Show this help message

Example:
    $0 --stratum allmapp --k 21,31,41 --maf 01_05,005_01

EOF
}
