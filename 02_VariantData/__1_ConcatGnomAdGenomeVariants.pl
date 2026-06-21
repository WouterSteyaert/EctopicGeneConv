#!/usr/bin/env perl
#===============================================================================
# Concat + mappability-restrict gnomAD v4.1 genome BEDs.
#
# Writes a single SLURM job that:
#   1. concatenates all per-chromosome *.sorted.bed files from __1_Fetch output
#   2. natural-sorts the combined BED
#   3. intersects with the GIAB all-mappable mask
#      (GRCh38_notinlowmappabilityall.bedsort.nochr.bed) -> VariantsToQuery.bed
#
# This produces the final 534 M SNV set used as the gnomAD input for downstream
# enrichment analyses (step 06 onwards).
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $JOBDIR     = $cfg{paths}{jobs_dir};
my $GENOMEDIR  = $cfg{paths}{gnomad_genome_dir};
my $REGIONS    = $cfg{paths}{regions_dir};
my $MAPP_BED   = "$REGIONS/GRCh38_notinlowmappabilityall.bedsort.nochr.bed";
my $PERL_MOD   = $cfg{slurm}{perl_module};
my $HTS_MOD    = $cfg{modules}{htslib};
my $BED_MOD    = $cfg{modules}{bedtools};

make_path($JOBDIR) unless -d $JOBDIR;

my $base = "$JOBDIR/ConcatGnomAdGenomeVariants";

open my $jf, '>', "$base.job" or die "Can't open $base.job: $!";
print $jf <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24400M
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=24:00:00

module load $PERL_MOD
module load $HTS_MOD
module load $BED_MOD

cd $GENOMEDIR

# 1. concatenate per-chr BEDs (sorted within chromosome from __1_Fetch)
cat *.sorted.bed > all.bed

# 2. natural-sort across chromosomes
sort -k1,1V -k2,2n all.bed > all.sorted.bed

# 3. restrict to GIAB high-mappability regions (final 534 M SNV set)
bedtools intersect -a all.sorted.bed -b $MAPP_BED > VariantsToQuery.bed
JOB
close $jf;
chmod 0755, "$base.job";

print "Wrote SLURM job: $base.job\n";
print "Submit:  sbatch $base.job\n";
