#!/usr/bin/env perl
#===============================================================================
# __7_SubmitContingencyStatisticsHotspot.pl
#
# Emit one SLURM job per (chromosome × template length × AF bin × mappability
# context × chunk) that runs __7_ComputeContingencyStatisticsHotspot on the
# hotspot-restricted contingency stats.  Chunked because the per-chrom worker
# becomes slow on long chromosomes; chunks are merged with
# __10_MergeChunkedRecombHotspots.
#
# Default mappability context: allmapp (matches Methods text, which restricts
# the hotspot test to autosomes within the standard mappable mask).  Override
# with --MappabilityContexts=allmapp,segdupmapp,nosegdupmapp.
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $MappabilityContexts = "allmapp";
my $HotspotBED          = "$cfg{paths}{regions_dir}/recomb_hotspots.bed.gz";
my $ChunkSize           = 10_000_000;
my $RegionType          = "RecombHotspots";

GetOptions(
    "MappabilityContexts=s" => \$MappabilityContexts,
    "HotspotBED=s"          => \$HotspotBED,
    "ChunkSize=i"           => \$ChunkSize,
    "RegionType=s"          => \$RegionType,
);
die "Hotspot BED not found: $HotspotBED (run __1_PrepareHotspotBED.sh first)\n"
    unless -f $HotspotBED;

my $JobDir   = $cfg{paths}{jobs_dir};
my $CodeDir  = $FindBin::Bin;
my $Hg38Dict = $cfg{paths}{reference_dict};
my $PerlMod  = $cfg{slurm}{perl_module};
my $HtsMod   = $cfg{modules}{htslib};
my @AfEdges  = split /,/, $cfg{enrichment}{af_edges};
my @Contexts = split /,/, $MappabilityContexts;
my @KValues  = split /,/, $cfg{repeat_lengths}{values};
my @Chroms   = (1..22);     # autosomes only (Methods text)

make_path($JobDir) unless -d $JobDir;

# Chromosome lengths from dict
my %ChromLengths;
open my $D, '<', $Hg38Dict or die "Can't open $Hg38Dict: $!\n";
while (<$D>) {
    chomp;
    my @F = split /\t/, $_;
    next unless defined $F[1];
    my $C = $F[1]; $C =~ s/SN://; $C =~ s/^chr//;
    next unless $C =~ /^\d+$/;
    my $L = $F[2]; $L =~ s/LN://;
    $ChromLengths{$C} = $L;
}
close $D;

my $NJobs = 0;

foreach my $Context (@Contexts) {
    foreach my $K (@KValues) {
        foreach my $Chrom (@Chroms) {
            next unless exists $ChromLengths{$Chrom};
            my $ChromLen = $ChromLengths{$Chrom};

            for (my $i = 1; $i < @AfEdges; $i++) {
                my $MinF = $AfEdges[$i-1];
                my $MaxF = $AfEdges[$i];
                my $MinPr = $MinF; $MinPr =~ s|/|~|; $MinPr =~ s/\.//;
                my $MaxPr = $MaxF; $MaxPr =~ s|/|~|; $MaxPr =~ s/\.//;

                my $ChunkNumber = 0;
                for (my $ChunkStart = 1; $ChunkStart < $ChromLen; $ChunkStart += $ChunkSize) {
                    my $ChunkEnd = $ChunkStart + $ChunkSize - 1;
                    $ChunkEnd    = $ChromLen if $ChunkEnd > $ChromLen;
                    my $Pad      = sprintf("%03d", $ChunkNumber);

                    my $base = "$JobDir/ContStatHotspot.$K.chr$Chrom.chunk$Pad.${MinPr}-${MaxPr}.${RegionType}.$Context";
                    open my $JOB, '>', "$base.job" or die "Can't open $base.job: $!\n";
                    print $JOB <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=36000M
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=06:00:00
#SBATCH --job-name=ContStatHotspot.$K.chr$Chrom.chunk$Pad

module load $PerlMod
module load $HtsMod

perl $CodeDir/__7_ComputeContingencyStatisticsHotspot.pl \\
    --QueryChrom=$Chrom \\
    --QueryIntervalStart=$ChunkStart \\
    --QueryIntervalEnd=$ChunkEnd \\
    --QueryRepLen=$K \\
    --QueryVarMinFreq=$MinF \\
    --QueryVarMaxFreq=$MaxF \\
    --MappabilityContext=$Context \\
    --HotspotBED=$HotspotBED \\
    --RegionType=$RegionType \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
                    close $JOB;
                    chmod 0755, "$base.job";
                    $NJobs++;
                    $ChunkNumber++;
                }
            }
        }
    }
}

print "Wrote $NJobs SLURM jobs to $JobDir\n";
print "Submit: for f in $JobDir/ContStatHotspot.*.job; do sbatch \$f; done\n";
print "After all jobs finish, merge with __10_MergeChunkedRecombHotspots.pl\n";
