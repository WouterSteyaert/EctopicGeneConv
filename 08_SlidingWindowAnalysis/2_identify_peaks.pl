#!/usr/bin/env perl
use strict;
use warnings;
use POSIX qw(log10);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

# =============================================================================
# Step 2: Peak identification on the corrected sliding-window enrichment
# (produced by 1_correct_enrichment.pl).
#
# Outlier 1Mb windows (non-overlapping) are those where log2(OR) exceeds the
# genome-wide mean + Z_SCORE * SD; adjacent outliers on the same chromosome
# are merged into peaks; peaks with fewer than MIN_GC_VAR conversion-compatible
# variants are dropped (MIN_GC_VAR matches what the paper reports).
# =============================================================================

my %cfg = load_config();

my $BASE   = $cfg{paths}{data_root};
my $INPUT  = "$BASE/sliding_window_corrected/enrichment_corrected_combined.tsv";
my $PEAKS  = "$BASE/sliding_window_corrected/peaks_corrected.tsv";
my $THRESH = "$BASE/sliding_window_corrected/peaks_thresholds_corrected.tsv";

my $Z_SCORE    = 2;
my $MIN_GC_VAR = 10;   # minimum gc_var per peak (matches paper Methods)

# --- Segmental duplication BED for overlap annotation ---
my $SEGDUP_BED = "$cfg{paths}{regions_dir}/GRCh38_segdups.sorted.merged.nochr.bed.gz";

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
    my $centro_dist;
    if ($mid >= $acen_start{$chr} && $mid <= $acen_end{$chr}) {
        $centro_dist = 0;
    } else {
        my $d1 = abs($mid - $acen_start{$chr});
        my $d2 = abs($mid - $acen_end{$chr});
        $centro_dist = $d1 < $d2 ? $d1 : $d2;
    }
    my $telo_dist = $mid < $chrsize/2 ? $mid : $chrsize - $mid;
    return "pericentromeric" if $centro_dist < 5_000_000;
    return "subtelomeric"    if $telo_dist < 5_000_000;
    return "interstitial";
}

sub min { return $_[0] < $_[1] ? $_[0] : $_[1]; }

# --- Load segmental duplication intervals ---
my %segdup;
if (-f $SEGDUP_BED) {
    open(my $sd_fh, "zcat $SEGDUP_BED |") or warn "Cannot read $SEGDUP_BED\n";
    if ($sd_fh) {
        while (<$sd_fh>) {
            chomp;
            my @f = split(/\t/);
            push @{$segdup{$f[0]}}, [$f[1], $f[2]];
        }
        close($sd_fh);
    }
}

sub segdup_overlap_fraction {
    my ($chr, $start, $end) = @_;
    return (0, "none") unless exists $segdup{$chr};
    my $overlap = 0;
    foreach my $seg (@{$segdup{$chr}}) {
        my $os = $seg->[0] > $start ? $seg->[0] : $start;
        my $oe = $seg->[1] < $end   ? $seg->[1] : $end;
        $overlap += ($oe - $os) if $oe > $os;
    }
    my $frac = $overlap / ($end - $start);
    my $label = $frac > 0.5 ? "segdup" : $frac > 0 ? "partial" : "none";
    return ($frac, $label);
}

# =============================================================================
# Read corrected enrichment data — non-overlapping 1Mb windows only
# =============================================================================
print STDERR "Reading corrected enrichment data...\n";

open(my $in, "<", $INPUT) or die "Cannot open $INPUT\n";
my $header = <$in>;

# Structure: {k}{maf} -> [ {chr, start, end, gc_var, gc_novar, nongc_var, nongc_novar, log2_or}, ... ]
my %data;
my %all_k;
my %all_maf;

while (<$in>) {
    chomp;
    my @f = split(/\t/);
    my ($k, $maf, $chr, $ws, $we, $mid, $gc_var, $gc_novar, $nongc_var, $nongc_novar,
        $gc_rate, $nongc_rate, $or, $log2_or, $cdist, $tdist, $region) = @f;

    # Non-overlapping windows only (start aligns to 1Mb grid)
    next unless (($ws - 1) % 1000000 == 0);

    $all_k{$k} = 1;
    $all_maf{$maf} = 1;

    push @{$data{$k}{$maf}}, {
        chr => $chr, start => $ws, end => $we,
        gc_var => $gc_var, gc_novar => $gc_novar,
        nongc_var => $nongc_var, nongc_novar => $nongc_novar,
        log2_or => $log2_or,
    };
}
close($in);

# =============================================================================
# Identify peaks per k x maf
# =============================================================================
print STDERR "Identifying peaks...\n";

open(my $peak_fh, ">", $PEAKS) or die "Cannot open $PEAKS\n";
print $peak_fh join("\t",
    "Chromosome", "Peak_start", "Peak_end", "Width_Mb",
    "Template_length_k", "MAF_bin", "Max_log2_OR", "Mean_log2_OR",
    "Chromosomal_region", "Centromere_dist_Mb", "Telomere_dist_Mb",
    "Segdup_overlap", "Segdup_fraction"
) . "\n";

open(my $thr_fh, ">", $THRESH) or die "Cannot open $THRESH\n";
print $thr_fh join("\t", "k", "maf", "n_windows", "mean_log2or", "sd_log2or",
    "threshold", "n_peaks", "n_peak_windows") . "\n";

my $total_peaks = 0;

foreach my $k (sort { $a <=> $b } keys %all_k) {
    foreach my $maf (sort keys %all_maf) {
        next unless exists $data{$k}{$maf};
        my @windows = @{$data{$k}{$maf}};
        next unless @windows > 10;

        # Compute mean and SD of log2_or
        my $sum = 0;
        my $sum_sq = 0;
        my $n = scalar @windows;
        foreach my $w (@windows) {
            $sum += $w->{log2_or};
            $sum_sq += $w->{log2_or} ** 2;
        }
        my $mean = $sum / $n;
        my $sd = sqrt($sum_sq / $n - $mean ** 2);
        my $threshold = $mean + $Z_SCORE * $sd;

        # Find peak windows
        my @peak_wins;
        foreach my $w (@windows) {
            if ($w->{log2_or} > $threshold && $w->{gc_var} >= $MIN_GC_VAR) {
                push @peak_wins, $w;
            }
        }

        # Merge adjacent peak windows per chromosome
        my %by_chr;
        foreach my $w (@peak_wins) {
            push @{$by_chr{$w->{chr}}}, $w;
        }

        my $n_peak_windows = scalar @peak_wins;
        my $n_peaks_km = 0;

        foreach my $chr (sort keys %by_chr) {
            my @sorted = sort { $a->{start} <=> $b->{start} } @{$by_chr{$chr}};

            my @merged;
            my $cur_start = $sorted[0]->{start};
            my $cur_end   = $sorted[0]->{end};
            my @cur_ors   = ($sorted[0]->{log2_or});

            for (my $i = 1; $i < @sorted; $i++) {
                if ($sorted[$i]->{start} <= $cur_end + 1) {
                    $cur_end = $sorted[$i]->{end} if $sorted[$i]->{end} > $cur_end;
                    push @cur_ors, $sorted[$i]->{log2_or};
                } else {
                    push @merged, [$cur_start, $cur_end, [@cur_ors]];
                    $cur_start = $sorted[$i]->{start};
                    $cur_end   = $sorted[$i]->{end};
                    @cur_ors   = ($sorted[$i]->{log2_or});
                }
            }
            push @merged, [$cur_start, $cur_end, [@cur_ors]];

            foreach my $peak (@merged) {
                my ($ps, $pe, $ors) = @$peak;
                my $width_mb = ($pe - $ps + 1) / 1e6;
                my $max_or = (sort { $b <=> $a } @$ors)[0];
                my $mean_or = 0;
                $mean_or += $_ for @$ors;
                $mean_or /= scalar @$ors;

                my $mid = ($ps + $pe) / 2;
                my $region = classify_region($mid, $chr, $chrom_sizes{$chr});
                my $cdist;
                if ($mid >= $acen_start{$chr} && $mid <= $acen_end{$chr}) {
                    $cdist = 0;
                } else {
                    my $d1 = abs($mid - $acen_start{$chr});
                    my $d2 = abs($mid - $acen_end{$chr});
                    $cdist = $d1 < $d2 ? $d1 : $d2;
                }
                my $tdist = min($mid, $chrom_sizes{$chr} - $mid);
                my ($sd_frac, $sd_label) = segdup_overlap_fraction($chr, $ps, $pe);

                print $peak_fh join("\t",
                    $chr, $ps, $pe, $width_mb,
                    $k, $maf,
                    sprintf("%.2f", $max_or), sprintf("%.2f", $mean_or),
                    $region,
                    sprintf("%.1f", $cdist / 1e6),
                    sprintf("%.1f", $tdist / 1e6),
                    $sd_label, sprintf("%.3f", $sd_frac)
                ) . "\n";

                $n_peaks_km++;
                $total_peaks++;
            }
        }

        print $thr_fh join("\t", $k, $maf, $n, sprintf("%.4f", $mean),
            sprintf("%.4f", $sd), sprintf("%.4f", $threshold),
            $n_peaks_km, $n_peak_windows) . "\n";
    }
}

close($peak_fh);
close($thr_fh);

print STDERR "\nTotal peaks: $total_peaks\n";
print STDERR "Written: $PEAKS\n";
print STDERR "Written: $THRESH\n";
