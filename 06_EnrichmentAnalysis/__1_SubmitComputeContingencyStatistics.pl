#!/usr/bin/env perl
#===============================================================================
# __1_SubmitComputeContingencyStatistics.pl
#
# For each combination of (chromosome × template length k × AF bin × mappability
# context), emit one SLURM job that runs __1_ComputeContingencyStatistics on
# that combination.
#
# Defaults run the genome-wide gnomAD case (--QueryVarList=gnomad/genome,
# --AfColumn=4) across the three mappability contexts.  For the other variant
# sets, override:
#
#   --QueryVarList=gwas/clean/All   --AfColumn=9       # GWAS
#   --QueryVarList=dnm/DenovoVars/ALL_ORIGIN --AfColumn=7   # DNMs
#   --QueryVarList=random3          --VarFreqMode=unstratified
#
# Iteration parameters (chromosomes, k values, AF bin edges, contexts) are
# read from the config.
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $QueryVarList = "gnomad/genome";
my $AfColumn     = 4;
my $VarFreqMode  = "stratified";    # 'stratified' (default) or 'unstratified'
GetOptions(
    "QueryVarList=s" => \$QueryVarList,
    "AfColumn=i"     => \$AfColumn,
    "VarFreqMode=s"  => \$VarFreqMode,
);

my $JobDir     = $cfg{paths}{jobs_dir};
my $CodeDir    = $FindBin::Bin;
my $Hg38Dict   = $cfg{paths}{reference_dict};
my $PerlMod    = $cfg{slurm}{perl_module};
my $HtsMod     = $cfg{modules}{htslib};
my $WindowSize = $cfg{enrichment}{window_size};
my $StepSize   = $cfg{enrichment}{step_size};
my @AfEdges    = split /,/, $cfg{enrichment}{af_edges};
my @Contexts   = split /,/, $cfg{enrichment}{contexts};
my @KValues    = split /,/, $cfg{repeat_lengths}{values};
my @Chroms     = split /,/, $cfg{genome}{chromosomes};

make_path($JobDir) unless -d $JobDir;

# Validate chrom in dict (so we never submit jobs for a non-existent contig)
my %DictChrom;
open my $D, '<', $Hg38Dict or die "Can't open $Hg38Dict: $!\n";
while (<$D>) {
    chomp;
    my @F = split /\t/, $_;
    next unless defined $F[1];
    my $C = $F[1]; $C =~ s/SN://; $C =~ s/^chr//;
    $DictChrom{$C}++ if ($C =~ /^\d+$/ || $C eq "X" || $C eq "Y");
}
close $D;

my $QVarListPrint = $QueryVarList; $QVarListPrint =~ s|/|~|g;
my $NJobs = 0;

foreach my $Context (@Contexts) {
    foreach my $K (@KValues) {
        foreach my $Chrom (@Chroms) {
            next unless exists $DictChrom{$Chrom};

            # Build the AF bin list:
            #   stratified -> for each adjacent edge pair
            #   unstratified -> one job with Min=0, Max=all
            my @Bins;
            if ($VarFreqMode eq "unstratified") {
                push @Bins, ["0", "all"];
            } else {
                for (my $i = 1; $i < @AfEdges; $i++) {
                    push @Bins, [$AfEdges[$i-1], $AfEdges[$i]];
                }
            }

            foreach my $Bin (@Bins) {
                my ($MinF, $MaxF) = @$Bin;
                my $MinPr = $MinF; $MinPr =~ s|/|~|; $MinPr =~ s/\.//;
                my $MaxPr = $MaxF; $MaxPr =~ s|/|~|; $MaxPr =~ s/\.//;

                my $base = "$JobDir/ContingencyStats.$K.$Chrom.$QVarListPrint.$MinPr.$MaxPr.$Context";
                open my $JOB, '>', "$base.job" or die "Can't open $base.job: $!\n";
                print $JOB <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12400M
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=03:00:00

module load $PerlMod
module load $HtsMod

perl $CodeDir/__1_ComputeContingencyStatistics.pl \\
    --QueryChrom=$Chrom \\
    --QueryWindowSize=$WindowSize \\
    --QueryStepSize=$StepSize \\
    --QueryVarList=$QueryVarList \\
    --QueryRepLen=$K \\
    --QueryVarMinFreq=$MinF \\
    --QueryVarMaxFreq=$MaxF \\
    --AfColumn=$AfColumn \\
    --MappabilityContext=$Context \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
                close $JOB;
                chmod 0755, "$base.job";
                $NJobs++;
            }
        }
    }
}

print "Wrote $NJobs SLURM jobs to $JobDir\n";
print "Submit: for f in $JobDir/ContingencyStats.*.job; do sbatch \$f; done\n";
