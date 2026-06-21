#!/usr/bin/env perl
#===============================================================================
# LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 11
# Calculate mean LD per GC/non-GC variant (1:1 paralogs only)
#
# MET RESAMPLING: max 5000 varianten per chromosoom
# Window sizes: 10kb, 25kb, 50kb, 100kb
#
# Output: LD_analysis_resampling_1to1/
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

my $ld_base  = "$BASE/LD_analysis_resampling_1to1";
my $bed_base = "$BASE/LD_analysis_resampling_1to1";
my $vcf_dir  = $cfg{paths}{thousand_genomes_dir} // "$BASE/1000G";
my $jobs_dir = "$ld_base/jobs";

my $PERL_MOD = $cfg{slurm}{perl_module};
my $PLINK_MOD = $cfg{modules}{plink} // "PLINK/1.9-x86_64";

my @mapp_categories = ("allmapp", "nosegdupmapp", "segdupmapp");
my @k_values        = split /,/, $cfg{repeat_lengths}{values};
my @populations     = ("EUR", "ALL");

my @maf_bins = (
    "0.001_0.005", "0.005_0.01", "0.01_0.05", "0.05_0.1", "0.1_0.5", "0.5_2"
);

my @ld_windows_kb = (10, 25, 50, 100);
my $max_variants_per_chr = 5000;

#===============================================================================
# Setup
#===============================================================================

make_path($jobs_dir) unless -d $jobs_dir;
for my $mapp (@mapp_categories) {
    for my $pop (@populations) {
        for my $win (@ld_windows_kb) {
            make_path("$ld_base/$mapp/$pop/win${win}kb") unless -d "$ld_base/$mapp/$pop/win${win}kb";
        }
    }
}

#===============================================================================
# Generate jobs
#===============================================================================

print "="x70 . "\n";
print "LD ANALYSE - 1:1 PARALOG PIPELINE - STAP 11: LD CALCULATION\n";
print "="x70 . "\n\n";
print "Window sizes: " . join(", ", map { "${_}kb" } @ld_windows_kb) . "\n";
print "Max variants per chr: $max_variants_per_chr (RESAMPLING)\n\n";

my $job_count = 0;

for my $mapp (@mapp_categories) {
    print "=== $mapp ===\n";

    for my $pop (@populations) {
        my $samples_file = "$vcf_dir/samples_${pop}.txt";

        for my $k (@k_values) {
            for my $bin (@maf_bins) {
                $job_count++;
                my $job_file = "$jobs_dir/ld_${mapp}_${pop}_k${k}_${bin}.sh";

                open(my $jf, ">", $job_file) or die "Cannot write $job_file: $!";

                my $windows_str = join(" ", @ld_windows_kb);

                print $jf <<"EOF";
#!/bin/bash
#SBATCH --job-name=l1_${mapp}_${pop}_k${k}
#SBATCH --output=$jobs_dir/ld_${mapp}_${pop}_k${k}_${bin}.out
#SBATCH --error=$jobs_dir/ld_${mapp}_${pop}_k${k}_${bin}.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load $PERL_MOD
module load $PLINK_MOD

echo "========================================"
echo "LD per variant 1:1 (RESAMPLING): mapp=$mapp pop=$pop k=$k maf=$bin"
echo "Window sizes: $windows_str kb"
echo "Max per chr: $max_variants_per_chr"
echo "========================================"

GC_VARIANTS="$bed_base/$mapp/k${k}_gc_$bin.bed"
NONGC_VARIANTS="$bed_base/$mapp/k${k}_nongc_$bin.bed"
SAMPLES="$samples_file"

WINDOWS=($windows_str)
MAX_VAR=$max_variants_per_chr

TMP="/tmp/l1_${mapp}_${pop}_k${k}_${bin}_\$\$"
mkdir -p \$TMP

# Initialize output files
for win in \${WINDOWS[\@]}; do
    GC_OUT="$ld_base/$mapp/$pop/win\${win}kb/ld_k${k}_gc_${bin}_pervar.txt"
    NONGC_OUT="$ld_base/$mapp/$pop/win\${win}kb/ld_k${k}_nongc_${bin}_pervar.txt"
    echo -e "chr\\tpos\\tmean_r2\\tn_pairs" > \$GC_OUT
    echo -e "chr\\tpos\\tmean_r2\\tn_pairs" > \$NONGC_OUT
done

for chr in {1..22}; do
    echo "Chr \$chr..."

    VCF="$vcf_dir/CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
    [[ ! -f "\$VCF" ]] && continue

    awk -v c=\$chr '\$1 == c' \$GC_VARIANTS > \$TMP/gc_chr.bed
    awk -v c=\$chr '\$1 == c' \$NONGC_VARIANTS > \$TMP/nongc_chr.bed

    GC_N=\$(wc -l < \$TMP/gc_chr.bed)
    NONGC_N=\$(wc -l < \$TMP/nongc_chr.bed)

    echo "  GC: \$GC_N, NonGC: \$NONGC_N"

    # RESAMPLING
    if [[ \$GC_N -gt \$MAX_VAR ]]; then
        shuf -n \$MAX_VAR \$TMP/gc_chr.bed | sort -k2,2n > \$TMP/gc_chr_s.bed
        mv \$TMP/gc_chr_s.bed \$TMP/gc_chr.bed
        echo "  Sampled GC to \$MAX_VAR"
    fi
    if [[ \$NONGC_N -gt \$MAX_VAR ]]; then
        shuf -n \$MAX_VAR \$TMP/nongc_chr.bed | sort -k2,2n > \$TMP/nongc_chr_s.bed
        mv \$TMP/nongc_chr_s.bed \$TMP/nongc_chr.bed
        echo "  Sampled NonGC to \$MAX_VAR"
    fi

    cat \$TMP/gc_chr.bed \$TMP/nongc_chr.bed | sort -k2,2n > \$TMP/all_chr.bed
    [[ ! -s \$TMP/all_chr.bed ]] && continue

    awk '{print \$2}' \$TMP/gc_chr.bed | sort -u > \$TMP/gc_pos.txt
    awk '{print \$2}' \$TMP/nongc_chr.bed | sort -u > \$TMP/nongc_pos.txt
    awk '{print "chr"\$1"\\t"\$2"\\t"\$3"\\tvar"NR}' \$TMP/all_chr.bed > \$TMP/range.bed

    for win in \${WINDOWS[\@]}; do
        echo "  Window: \${win}kb"

        GC_OUT="$ld_base/$mapp/$pop/win\${win}kb/ld_k${k}_gc_${bin}_pervar.txt"
        NONGC_OUT="$ld_base/$mapp/$pop/win\${win}kb/ld_k${k}_nongc_${bin}_pervar.txt"

        plink --vcf \$VCF \\
              --keep \$SAMPLES \\
              --extract range \$TMP/range.bed \\
              --snps-only just-acgt \\
              --r2 \\
              --ld-window-kb \$win \\
              --ld-window 99999 \\
              --ld-window-r2 0 \\
              --out \$TMP/ld_chr_\${win} \\
              --threads 4 \\
              --allow-extra-chr 2>/dev/null

        [[ ! -f \$TMP/ld_chr_\${win}.ld ]] && continue

        # GC: partners within +/- FLANK bp of focal are excluded
        # (k-flanking is identical between donor and acceptor by definition)
        # FLANK = (k - 1) / 2; here k=${k} so FLANK=@{[int(($k - 1) / 2)]}
        awk -v chr=\$chr -v flank=@{[int(($k - 1) / 2)]} '
        BEGIN {OFS="\\t"}
        NR == FNR {gc[\$1] = 1; next}
        FNR == 1 {next}
        {
            pos_a = \$2; pos_b = \$5; r2 = \$7
            d = pos_b - pos_a; if (d < 0) d = -d
            if (d <= flank) next
            if (pos_a in gc) { sum[pos_a] += r2; count[pos_a]++ }
            if (pos_b in gc) { sum[pos_b] += r2; count[pos_b]++ }
        }
        END {
            for (pos in sum) {
                if (count[pos] > 0) print chr, pos, sum[pos]/count[pos], count[pos]
            }
        }' \$TMP/gc_pos.txt \$TMP/ld_chr_\${win}.ld >> \$GC_OUT

        awk -v chr=\$chr '
        BEGIN {OFS="\\t"}
        NR == FNR {nongc[\$1] = 1; next}
        FNR == 1 {next}
        {
            pos_a = \$2; pos_b = \$5; r2 = \$7
            if (pos_a in nongc) { sum[pos_a] += r2; count[pos_a]++ }
            if (pos_b in nongc) { sum[pos_b] += r2; count[pos_b]++ }
        }
        END {
            for (pos in sum) {
                if (count[pos] > 0) print chr, pos, sum[pos]/count[pos], count[pos]
            }
        }' \$TMP/nongc_pos.txt \$TMP/ld_chr_\${win}.ld >> \$NONGC_OUT

        rm -f \$TMP/ld_chr_\${win}.*
    done
done

rm -rf \$TMP

echo ""
echo "========================================"
echo "Summary:"
for win in \${WINDOWS[\@]}; do
    GC_OUT="$ld_base/$mapp/$pop/win\${win}kb/ld_k${k}_gc_${bin}_pervar.txt"
    NONGC_OUT="$ld_base/$mapp/$pop/win\${win}kb/ld_k${k}_nongc_${bin}_pervar.txt"
    GC_TOTAL=\$(tail -n +2 \$GC_OUT | wc -l)
    NONGC_TOTAL=\$(tail -n +2 \$NONGC_OUT | wc -l)
    echo "Window \${win}kb: GC=\$GC_TOTAL NonGC=\$NONGC_TOTAL"
done
echo "Done!"
EOF
                close($jf);
                chmod 0755, $job_file;
            }
        }
        print "  $pop: " . scalar(@k_values) * scalar(@maf_bins) . " jobs\n";
    }
}

print "\n" . "="x70 . "\n";
print "TOTAAL: $job_count jobs\n";
print "="x70 . "\n";
print "\nSubmit:\n  for f in $jobs_dir/ld_*.sh; do sbatch -M shinx \$f; done\n\n";
