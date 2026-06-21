#!/usr/bin/env perl
#===============================================================================
# __2_SubmitComputeTractLengthWgs.pl
#
# Sample N (default 5000) gnomAD variants per AF bin from the post-mappability
# VariantsToQuery.bed (step 02), then emit SLURM jobs that run
# __2_ComputeTractLengthWgs on each.  Variants are batched 100 per job; jobs
# already containing previously-selected variants are skipped, allowing
# resubmission without redoing previous work.
#
# Inputs (config):
#   <paths>gnomad_genome_dir/VariantsToQuery.bed
#   <paths>tractlengths_sum_dir            (records previously-submitted variants)
#   <paths>jobs_dir
#   <enrichment>af_edges                   (AF bin edges; reused from step 06)
#   <tract>n_variants_per_bin              (5000 by default)
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $JobDir          = $cfg{paths}{jobs_dir};
my $TractDir        = $cfg{paths}{tractlengths_wgs_dir};
my $TractSumDir     = $cfg{paths}{tractlengths_sum_dir};
my $GnomadDir       = $cfg{paths}{gnomad_genome_dir};
my $QueryVarsFp     = "$GnomadDir/VariantsToQuery.bed";
my $CodeDir         = $FindBin::Bin;
my $PerlMod         = $cfg{slurm}{perl_module};
my $HtsMod          = $cfg{modules}{htslib};
my $NrPerFreq       = $cfg{tract}{n_variants_per_bin} || 5000;
my @AfEdges         = split /,/, $cfg{enrichment}{af_edges};

my $RandomRunNr     = int(rand(9999999999999));
my $InputVarsFp     = "$TractSumDir/Variant.$RandomRunNr.input";

for my $d ($JobDir, $TractDir, $TractSumDir) { make_path($d) unless -d $d; }

# Don't reselect variants that previous runs already covered
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

# Stratified-sampling loop: shuf 500k rows, bucket into AF bins, repeat until
# every bin has reached the target count.
print "Select random variants\n";
my %SelectedVariants;
my $SelectionComplete = "no";

do {
    my @RawLines = `shuf -n 500000 $QueryVarsFp`;
    foreach my $Line (@RawLines) {
        chomp $Line;
        my @F = split(/\t/, $Line);
        my ($SelectedVariant, $SelectedFrequency) = ($F[3], $F[4]);

        for (my $I = 1; $I < @AfEdges; $I++) {
            my $Min = $AfEdges[$I-1];
            my $Max = $AfEdges[$I];
            my $MinPr = $Min; $MinPr =~ s|/|~|; $MinPr =~ s/\.//;
            my $MaxPr = $Max; $MaxPr =~ s|/|~|; $MaxPr =~ s/\.//;
            my $Key   = "${MinPr}_${MaxPr}";

            if ($SelectedFrequency >= $Min && $SelectedFrequency < $Max) {
                next if exists $InputVariants{$SelectedVariant};
                next if exists $SelectedVariants{$Key}{$SelectedVariant};
                next if (scalar keys %{$SelectedVariants{$Key} // {}} >= $NrPerFreq);
                $SelectedVariants{$Key}{$SelectedVariant} = $SelectedFrequency;
            }
        }
    }

    $SelectionComplete = "yes";
    foreach my $Key (sort keys %SelectedVariants) {
        $SelectionComplete = "no"
            if scalar keys %{$SelectedVariants{$Key}} < $NrPerFreq;
    }
} while ($SelectionComplete eq "no");

# Emit jobs (100 variants per job)
open my $V, '>', $InputVarsFp or die "Can't open $InputVarsFp: $!\n";
foreach my $Key (sort keys %SelectedVariants) {
    my ($VarId, $JobId) = (0, 0);
    foreach my $VARIANT (sort keys %{$SelectedVariants{$Key}}) {
        $JobId = int($VarId / 100);
        $VarId++;

        my $base   = "$JobDir/ComputeTractLengthWgs.$Key.$JobId";
        my $JobFp  = "$base.job";
        my $FREQ   = $SelectedVariants{$Key}{$VARIANT};

        unless (-f $JobFp) {
            open my $JOB, '>', $JobFp or die "Can't open $JobFp: $!\n";
            print $JOB <<"HEADER";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16400M
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
        print $JOB "perl $CodeDir/__2_ComputeTractLengthWgs.pl "
                 . "--VARIANT=$VARIANT --FREQ=$FREQ "
                 . "--ConfigFile=$cfg{_meta}{config_file} "
                 . "--ProjectRoot=$cfg{_meta}{project_root}\n";
        close $JOB;
        print $V "$VARIANT\t$FREQ\n";
    }
}
close $V;

print "Wrote variant list to $InputVarsFp\n";
print "Submit: for f in $JobDir/ComputeTractLengthWgs.*.job; do sbatch \$f; done\n";
