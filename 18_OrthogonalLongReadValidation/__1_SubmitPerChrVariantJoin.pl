#!/usr/bin/env perl
#===============================================================================
# __1_SubmitPerChrVariantJoin.pl
#
# Emit one SLURM job per chromosome (chr1-22, chrX, chrY) that runs
# __1_PerChrVariantJoin.py — the per-chromosome variant-level join between
# gnomAD genome SNVs and GA4K PacBio HiFi joint calls.
#
# Memory 32 GB per chr: chr1/chr2 pos_to_k dict (~15-18 GB) + GA4K preload (~2 GB).
#
# Usage:
#   perl __1_SubmitPerChrVariantJoin.pl --ConfigFile=../00_Configuration/config.GRCh38.ini
#   for f in $PROJECT_ROOT/geneconv_complete/jobs/orthogonal_val_chr*.job; do sbatch $f; done
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $JobDir       = $cfg{paths}{jobs_dir};
my $CodeDir      = $FindBin::Bin;
my $Hg38Dict     = $cfg{paths}{reference_dict};
my $ProjectRoot  = $cfg{_meta}{project_root};

make_path($JobDir) unless -d $JobDir;

# Read chromosomes from the reference dict (autosomes + X + Y)
my @Chroms;
open my $D, '<', $Hg38Dict or die "Can't open $Hg38Dict: $!\n";
while (<$D>) {
    chomp;
    my @F = split /\t/, $_;
    next unless defined $F[1];
    my $C = $F[1]; $C =~ s/SN://; $C =~ s/^chr//;
    push @Chroms, $C if ($C =~ /^\d+$/ || $C eq "X" || $C eq "Y");
}
close $D;

my $NJobs = 0;
foreach my $Chrom (@Chroms) {
    my $base = "$JobDir/orthogonal_val_chr$Chrom";
    open my $JOB, '>', "$base.job" or die "Can't open $base.job: $!\n";
    print $JOB <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=06:00:00

export PROJECT_ROOT=$ProjectRoot

# Worker uses Python stdlib only (gzip, re, os, bisect, math, collections, argparse)
# — no conda activation required.
/usr/bin/python3 $CodeDir/__1_PerChrVariantJoin.py --chr=$Chrom
JOB
    close $JOB;
    chmod 0755, "$base.job";
    $NJobs++;
}
print "Wrote $NJobs SLURM jobs to $JobDir/orthogonal_val_chr*.job\n";
print "Submit with: for f in \$PROJECT_ROOT/geneconv_complete/jobs/orthogonal_val_chr*.job; do sbatch \$f; done\n";
