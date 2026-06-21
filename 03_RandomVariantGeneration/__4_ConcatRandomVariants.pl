#!/usr/bin/env perl
#===============================================================================
# __4_ConcatRandomVariants.pl
#
# Emit one SLURM job per random set that:
#   1. concatenates all per-chromosome <chr>.bedsort.bed inside the set;
#   2. natural-sorts the combined BED;
#   3. intersects with the GIAB v3.5 all-mappable mask
#      (regions_dir/GRCh38_notinlowmappabilityall.bedsort.nochr.bed)
#      to produce VariantsToQuery.bed (the analog of the gnomAD final set).
#
# Random set names are read from <random>set_names in the config.
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $DataRoot   = $cfg{paths}{data_root};
my $JobDir     = $cfg{paths}{jobs_dir};
my $Regions    = $cfg{paths}{regions_dir};
my $PerlMod    = $cfg{slurm}{perl_module};
my $HtsMod     = $cfg{modules}{htslib};
my $BedMod     = $cfg{modules}{bedtools};
my $MappBed    = "$Regions/GRCh38_notinlowmappabilityall.bedsort.nochr.bed";
my @RandomSets = split /,/, ($cfg{random}{set_names} // "random1,random2,random3");

make_path($JobDir) unless -d $JobDir;

foreach my $RandomSet (@RandomSets) {
    my $SetDir   = "$DataRoot/$RandomSet";
    my $base     = "$JobDir/ConcatRandomVariants.$RandomSet";
    my $JobFp    = "$base.job";

    open my $JOB, '>', $JobFp or die "Can't open $JobFp: $!\n";
    print $JOB <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24400M
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=24:00:00

module load $PerlMod
module load $HtsMod
module load $BedMod

cd $SetDir

cat *.bedsort.bed > all.bed
sort -k1,1V -k2,2n all.bed > all.sorted.bed
bedtools intersect -a all.sorted.bed -b $MappBed > VariantsToQuery.bed
JOB
    close $JOB;
    chmod 0755, $JobFp;
}

print "Wrote " . scalar(@RandomSets) . " concat jobs to $JobDir\n";
print "Submit: for f in $JobDir/ConcatRandomVariants.*.job; do sbatch \$f; done\n";
