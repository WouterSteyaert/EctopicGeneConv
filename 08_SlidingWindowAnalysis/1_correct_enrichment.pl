#!/usr/bin/env perl
use strict;
use warnings;
use POSIX qw(log10);
use File::Basename;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

# =============================================================================
# Step 1 (corrected): Sliding Window Enrichment with higher-k control group
# correction (cf. __10_SummarizeContingencyStatistics.pl in step 06).
#
# For each window at k, subtract GC-compatible counts from all higher k from
# the non-GC columns so the control group contains ONLY positions that are
# not GC-compatible at ANY k value (consistent with the genome-wide analysis).
#
# This is the variant of the sliding-window analysis used in the paper.
# =============================================================================

my %cfg = load_config();

my $BASE       = $cfg{paths}{stats_dir};
my $STATS_DIR  = "$BASE/stats_allmapp";
my $OUTPUT_DIR = "$cfg{paths}{data_root}/sliding_window_corrected";
make_path($OUTPUT_DIR) unless -d $OUTPUT_DIR;

my $window_size = $cfg{enrichment}{window_size} || 1_000_000;
my $step_size   = 50_000;     # sliding step (smaller than enrichment.step_size)
my $source      = "gnomad~genome";

my @k_values = split /,/, $cfg{repeat_lengths}{values};
my @maf_bins;
{
    my @e = split /,/, $cfg{enrichment}{af_edges};
    for (my $i = 1; $i < @e; $i++) {
        my $lo = $e[$i-1]; $lo =~ s|/|~|; $lo =~ s/\.//;
        my $hi = $e[$i];   $hi =~ s|/|~|; $hi =~ s/\.//;
        push @maf_bins, "${lo}_${hi}";
    }
}

my @chroms = (1..22, 'X', 'Y');

# --- GRCh38 Centromere boundaries ---
my %acen_start = (
    1=>121700000, 2=>91800000,  3=>87800000,  4=>48200000,
    5=>46100000,  6=>58500000,  7=>58100000,  8=>43200000,
    9=>42200000,  10=>38000000, 11=>51000000, 12=>33200000,
    13=>16500000, 14=>16100000, 15=>17500000, 16=>35300000,
    17=>22700000, 18=>15400000, 19=>24200000, 20=>25700000,
    21=>10900000, 22=>13700000, X=>58100000,  Y=>10300000,
);
my %acen_end = (
    1=>125100000, 2=>96000000,  3=>94000000,  4=>51800000,
    5=>51400000,  6=>62600000,  7=>62100000,  8=>47200000,
    9=>45500000,  10=>41600000, 11=>55800000, 12=>37800000,
    13=>18900000, 14=>18200000, 15=>20500000, 16=>38400000,
    17=>27400000, 18=>21500000, 19=>28100000, 20=>30400000,
    21=>13000000, 22=>17400000, X=>63800000,  Y=>10600000,
);
my %chrom_sizes = (
    1=>248956422, 2=>242193529, 3=>198295559, 4=>190214555,
    5=>181538259, 6=>170805979, 7=>159345973, 8=>145138636,
    9=>138394717, 10=>133797422, 11=>135086622, 12=>133275309,
    13=>114364328, 14=>107043718, 15=>101991189, 16=>90338345,
    17=>83257441, 18=>80373285, 19=>58617616, 20=>64444167,
    21=>46709983, 22=>50818468, X=>156040895, Y=>57227415,
);

sub classify_region {
    my ($mid, $chr, $chrsize) = @_;
    my $centro_mid = ($acen_start{$chr} + $acen_end{$chr}) / 2;
    my $telo_dist = $mid < $centro_mid ? $mid : $chrsize - $mid;
    my $centro_dist;
    if ($mid >= $acen_start{$chr} && $mid <= $acen_end{$chr}) {
        $centro_dist = 0;
    } else {
        my $d1 = abs($mid - $acen_start{$chr});
        my $d2 = abs($mid - $acen_end{$chr});
        $centro_dist = $d1 < $d2 ? $d1 : $d2;
    }
    return "pericentromeric" if $centro_dist < 5_000_000;
    return "subtelomeric"    if $telo_dist < 5_000_000;
    return "interstitial";
}

sub min { return $_[0] < $_[1] ? $_[0] : $_[1]; }

# =============================================================================
# Main: process each MAF bin x chromosome
# =============================================================================

# Combined output
my $combined_file = "$OUTPUT_DIR/enrichment_corrected_combined.tsv";
open(my $combo_fh, ">", $combined_file) or die "Cannot open $combined_file\n";
print $combo_fh join("\t", "k", "maf_bin", "chr", "window_start", "window_end", "window_mid",
    "gc_var", "gc_novar", "nongc_var", "nongc_novar",
    "gc_rate", "nongc_rate", "odds_ratio", "log2_or",
    "centromere_dist", "telomere_dist", "region_type") . "\n";

foreach my $maf (@maf_bins) {
    print STDERR "=== MAF: $maf ===\n";

    foreach my $chr (@chroms) {
        my $chrsize = $chrom_sizes{$chr};

        # --- Load per-window data for ALL k values ---
        # Key: window_start -> { k -> [gc_var, gc_novar, nongc_var, nongc_novar] }
        my %windows;
        my @window_order;
        my %seen_window;

        foreach my $k (@k_values) {
            my $file = "$STATS_DIR/$k/${chr}_${window_size}_${step_size}_${source}_${maf}.1.${chrsize}.txt";
            unless (-f $file) {
                print STDERR "  WARNING: $file not found\n";
                next;
            }

            open(my $fh, "<", $file) or die "Cannot open $file\n";
            my $header = <$fh>;  # skip header
            while (<$fh>) {
                chomp;
                my @c = split(/\t/);
                my $ws = $c[1];
                $windows{$ws}{$k} = [$c[3], $c[5], $c[4], $c[6]];
                # gc_var=$c[3], gc_novar=$c[5], nongc_var=$c[4], nongc_novar=$c[6]
                unless ($seen_window{$ws}) {
                    push @window_order, $ws;
                    $seen_window{$ws} = 1;
                }
            }
            close($fh);
        }

        # --- Apply correction and output ---
        foreach my $k (@k_values) {
            foreach my $ws (@window_order) {
                next unless exists $windows{$ws}{$k};

                my ($gc_var, $gc_novar, $nongc_var, $nongc_novar) = @{$windows{$ws}{$k}};

                # Correction: subtract higher-k GC counts from non-GC columns
                foreach my $k2 (@k_values) {
                    if ($k < $k2 && exists $windows{$ws}{$k2}) {
                        $nongc_var    -= $windows{$ws}{$k2}[0];  # gc_var of higher k
                        $nongc_novar  -= $windows{$ws}{$k2}[1];  # gc_novar of higher k
                    }
                }

                # Safety: clamp to 0 (rounding/edge effects)
                $nongc_var    = 0 if $nongc_var < 0;
                $nongc_novar  = 0 if $nongc_novar < 0;

                # Skip if no data in either group
                next if ($gc_var + $gc_novar == 0 || $nongc_var + $nongc_novar == 0);

                # Calculate rates and OR
                my $gc_rate    = $gc_var / ($gc_var + $gc_novar);
                my $nongc_rate = $nongc_var / ($nongc_var + $nongc_novar);

                my $pseudo = 0.5;
                my $or = (($gc_var + $pseudo) * ($nongc_novar + $pseudo)) /
                         (($gc_novar + $pseudo) * ($nongc_var + $pseudo));
                my $log2_or = log($or) / log(2);

                # Window metadata
                my $we = $ws + $window_size - 1;
                my $mid = ($ws + $we) / 2;
                my $centro_dist;
                if ($mid >= $acen_start{$chr} && $mid <= $acen_end{$chr}) {
                    $centro_dist = 0;
                } else {
                    my $d1 = abs($mid - $acen_start{$chr});
                    my $d2 = abs($mid - $acen_end{$chr});
                    $centro_dist = $d1 < $d2 ? $d1 : $d2;
                }
                my $telo_dist = min($mid, $chrsize - $mid);
                my $region = classify_region($mid, $chr, $chrsize);

                print $combo_fh join("\t",
                    $k, $maf, $chr, $ws, $we, int($mid),
                    $gc_var, $gc_novar, $nongc_var, $nongc_novar,
                    sprintf("%.6f", $gc_rate),
                    sprintf("%.6f", $nongc_rate),
                    sprintf("%.4f", $or),
                    sprintf("%.4f", $log2_or),
                    int($centro_dist),
                    int($telo_dist),
                    $region
                ) . "\n";
            }
        }

        print STDERR "  chr$chr done\n";
    }
}

close($combo_fh);
print STDERR "\nWritten: $combined_file\n";

# Also write per-k per-maf files for convenience
print STDERR "Splitting into per-k per-maf files...\n";
open(my $in_fh, "<", $combined_file) or die;
my $header_line = <$in_fh>;
my %out_handles;
while (<$in_fh>) {
    my @f = split(/\t/);
    my $key = "k$f[0]_maf$f[1]";
    unless ($out_handles{$key}) {
        my $fn = "$OUTPUT_DIR/enrichment_corrected_${key}.tsv";
        open($out_handles{$key}, ">", $fn) or die;
        # Write header without k and maf_bin columns
        my @hdr = split(/\t/, $header_line);
        print {$out_handles{$key}} join("\t", @hdr[2..$#hdr]);
    }
    print {$out_handles{$key}} join("\t", @f[2..$#f]);
}
close($in_fh);
close($_) for values %out_handles;

print STDERR "Done.\n";
