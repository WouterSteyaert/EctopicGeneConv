#!/usr/bin/env perl
#===============================================================================
# __1_SubmitComputeTractLengthRandom.pl
#
# Sample N (default 5000) random control positions from the random-set
# VariantsToQuery.bed (step 03) and emit SLURM jobs that run
# __1_ComputeTractLengthRandom3 on each.  Random control positions are
# matched in number to one AF bin of the gnomAD analysis (5000), so they form
# a baseline of equal sampling depth.
#
# Inputs (config):
#   <paths>data_root/<random_set>/VariantsToQuery.bed   (default random3)
#   <paths>tractlengths_random_dir
#   <paths>jobs_dir
#   <tract>n_random_variants
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $RandomSet = "random3";
GetOptions("RandomSet=s" => \$RandomSet);

my $JobDir       = $cfg{paths}{jobs_dir};
my $TractDir     = $cfg{paths}{tractlengths_random_dir};
my $TractSumDir  = "$cfg{paths}{data_root}/tractlengthsSumRandom3";
my $QueryVarsFp  = "$cfg{paths}{data_root}/$RandomSet/VariantsToQuery.bed";
my $CodeDir      = $FindBin::Bin;
my $PerlMod      = $cfg{slurm}{perl_module};
my $HtsMod       = $cfg{modules}{htslib};
my $NrVariants   = $cfg{tract}{n_random_variants} || 5000;

my $RandomRunNr  = int(rand(9999999999999));
my $InputVarsFp  = "$TractSumDir/Variant.$RandomRunNr.input";

for my $d ($JobDir, $TractDir, $TractSumDir) { make_path($d) unless -d $d; }

# Skip variants already selected by previous runs
my %InputVariants;
if (-d $TractSumDir) {
    opendir(my $SUM, $TractSumDir);
    while (my $File = readdir($SUM)) {
        next unless $File =~ /^Variant\.\d+\.input$/;
        open my $F, '<', "$TractSumDir/$File" or next;
        while (<$F>) {
            chomp;
            s/\t.*//;
            $InputVariants{$_} = undef;
        }
        close $F;
    }
    closedir($SUM);
}

print "Select random variants\n";
my %SelectedVariants;
do {
    my @RawLines = `shuf -n 500000 $QueryVarsFp`;
    foreach my $Line (@RawLines) {
        chomp $Line;
        my @F = split(/\t/, $Line);
        my $SelectedVariant = $F[3];
        last if scalar keys %SelectedVariants >= $NrVariants;
        next if exists $InputVariants{$SelectedVariant};
        $SelectedVariants{$SelectedVariant} = undef;
    }
} while (scalar keys %SelectedVariants < $NrVariants);

# Batch 100 variants per job
my ($VarId, $JobId) = (0, 0);
open my $V, '>', $InputVarsFp or die "Can't open $InputVarsFp: $!\n";
foreach my $VARIANT (sort keys %SelectedVariants) {
    $JobId = int($VarId / 100);
    $VarId++;

    my $base  = "$JobDir/ComputeTractLengthRandom.$JobId";
    my $JobFp = "$base.job";

    unless (-f $JobFp) {
        open my $JOB, '>', $JobFp or die "Can't open $JobFp: $!\n";
        print $JOB <<"HEADER";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12400M
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=72:00:00

module load $PerlMod
module load $HtsMod

HEADER
        close $JOB;
        chmod 0755, $JobFp;
    }

    open my $JOB, '>>', $JobFp or die "Can't open $JobFp: $!\n";
    print $JOB "perl $CodeDir/__1_ComputeTractLengthRandom3.pl "
             . "--VARIANT=$VARIANT "
             . "--ConfigFile=$cfg{_meta}{config_file} "
             . "--ProjectRoot=$cfg{_meta}{project_root}\n";
    close $JOB;
    print $V "$VARIANT\n";
}
close $V;

print "Wrote variant list to $InputVarsFp\n";
print "Submit: for f in $JobDir/ComputeTractLengthRandom.*.job; do sbatch \$f; done\n";
