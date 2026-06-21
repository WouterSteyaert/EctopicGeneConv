#!/usr/bin/env perl
#===============================================================================
# __7_ComputeContingencyStatisticsHotspot.pl
#
# Recombination-hotspot variant of the step-06 contingency-stats worker.
# Identical counting logic, with the additional restriction that mappable
# positions must overlap the deCODE recombination-hotspot BED (intervals
# with sex-averaged rate >= 10 cM/Mb; produced by __1_PrepareHotspotBED.sh).
#
# Output naming mirrors step 06 but with a .<RegionType>.txt suffix so the
# step-06 summarizer (or the dedicated SupTable_RecombHotspots.R) can read
# both the all-positions and hotspot-only tables.
#
# Usage:
#   perl __7_ComputeContingencyStatisticsHotspot.pl \
#       --QueryChrom=1 --QueryIntervalStart=1 --QueryIntervalEnd=10000000 \
#       --QueryRepLen=21 --QueryVarMinFreq=0.001 --QueryVarMaxFreq=0.005 \
#       --HotspotBED=$PROJECT_ROOT/geneconv_complete/regions/recomb_hotspots.bed.gz \
#       --ConfigFile=../00_Configuration/config.GRCh38.ini
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use IO::Zlib;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $QueryChrom         = 1;
my $QueryIntervalStart = 0;
my $QueryIntervalEnd   = 0;
my $QueryWindowSize    = $cfg{enrichment}{window_size} || 1_000_000;
my $QueryStepSize      = $cfg{enrichment}{step_size}   || 1_000_000;
my $QueryVarList       = "gnomad/genome";
my $QueryRepLen        = 21;
my $QueryVarMinFreq    = 0;
my $QueryVarMaxFreq    = "all";
my $AfColumn           = 4;
my $MappabilityContext = "allmapp";
my $HotspotBED         = "";       # required
my $RegionType         = "RecombHotspots";

GetOptions(
    "QueryChrom=s"          => \$QueryChrom,
    "QueryIntervalStart=i"  => \$QueryIntervalStart,
    "QueryIntervalEnd=i"    => \$QueryIntervalEnd,
    "QueryWindowSize=i"     => \$QueryWindowSize,
    "QueryStepSize=i"       => \$QueryStepSize,
    "QueryVarList=s"        => \$QueryVarList,
    "QueryRepLen=i"         => \$QueryRepLen,
    "QueryVarMinFreq=s"     => \$QueryVarMinFreq,
    "QueryVarMaxFreq=s"     => \$QueryVarMaxFreq,
    "AfColumn=i"            => \$AfColumn,
    "MappabilityContext=s"  => \$MappabilityContext,
    "HotspotBED=s"          => \$HotspotBED,
    "RegionType=s"          => \$RegionType,
);
die "--HotspotBED=<bed.gz> required\n" unless $HotspotBED;
die "$HotspotBED not found\n"          unless -f $HotspotBED;
die "--MappabilityContext must be allmapp|segdupmapp|nosegdupmapp\n"
    unless $MappabilityContext =~ /^(allmapp|segdupmapp|nosegdupmapp)$/;

my %ContextToMask = (
    allmapp      => "GRCh38_notinlowmappabilityall.bedsort.nochr.bed.gz",
    segdupmapp   => "GRCh38_notinlowmappabilitysegdups.bedsort.nochr.bed.gz",
    nosegdupmapp => "GRCh38_notinlowmappabilitynosegdups.bedsort.nochr.bed.gz",
);
my $MappableBases = "$cfg{paths}{regions_dir}/$ContextToMask{$MappabilityContext}";
my $Hg38DictFp    = $cfg{paths}{reference_dict};
my $DataRoot      = $cfg{paths}{data_root};
my $ConvDir       = $cfg{paths}{conv_dir};
my $GcdiffsDir    = $cfg{paths}{gcdiffs_dir};
my $StatsDir      = $cfg{paths}{stats_dir};

my @Bases = ("A", "C", "G", "T");

my $QueryVarListPrint = $QueryVarList; $QueryVarListPrint =~ s|/|~|g;
my $MinPr             = $QueryVarMinFreq; $MinPr =~ s|/|~|; $MinPr =~ s/\.//;
my $MaxPr             = $QueryVarMaxFreq; $MaxPr =~ s|/|~|; $MaxPr =~ s/\.//;

my $VariantFilePath = "$DataRoot/$QueryVarList/$QueryChrom.sorted.bed.gz";
my $StatsSubdir     = "$StatsDir/stats_$MappabilityContext/$QueryRepLen";
my $OutputFilePath  = "$StatsSubdir/${QueryChrom}_${QueryWindowSize}_${QueryStepSize}_${QueryVarListPrint}_${MinPr}_${MaxPr}.${QueryIntervalStart}.${QueryIntervalEnd}.${RegionType}.txt";

my ($PositionsGeneConvFilePath, $ChangesGeneConvFilePath);
my $AllChangesGeneConvFilePath = "$ConvDir/$QueryRepLen/$QueryChrom.repsum.geneconv.sorted.bed.gz";

if ($QueryRepLen == 91) {
    $PositionsGeneConvFilePath = $AllChangesGeneConvFilePath;
    $ChangesGeneConvFilePath   = $AllChangesGeneConvFilePath;
} else {
    my $NextK = ($QueryRepLen < 21) ? $QueryRepLen + 2 : $QueryRepLen + 10;
    $PositionsGeneConvFilePath = "$GcdiffsDir/${QueryRepLen}_${NextK}/$QueryChrom.positions.gcdiffs.sorted.bed.gz";
    $ChangesGeneConvFilePath   = "$GcdiffsDir/${QueryRepLen}_${NextK}/$QueryChrom.changes.gcdiffs.sorted.bed.gz";
}

make_path($StatsSubdir) unless -d $StatsSubdir;

my %ChromLengths;
open my $D, '<', $Hg38DictFp or die "Can't open $Hg38DictFp: $!\n";
while (<$D>) {
    chomp;
    my @F = split(/\t/, $_);
    next unless defined $F[1];
    my $C = $F[1]; $C =~ s/SN://; $C =~ s/^chr//;
    next unless ($C =~ /^\d+$/ || $C eq "X" || $C eq "Y");
    my $L = $F[2]; $L =~ s/LN://;
    $ChromLengths{$C} = $L;
}
close $D;
die "Chromosome $QueryChrom not in dict\n" unless exists $ChromLengths{$QueryChrom};

open my $O, '>', $OutputFilePath or die "Can't open $OutputFilePath: $!\n";
print $O "QueryChrom\tWindowStart\tWindowEnd\t"
       . "NrOfVarPosConvPos\tNrOfVarPosNoConvPos\tNrOfNoVarPosConvPos\tNrOfNoVarPosNoConvPos\t"
       . "NrOfVarMatchConv\tNrOfVarMatchNoConv\tNrOfNoVarMatchConv\tNrOfNoVarMatchNoConv\t"
       . "NrOfConvPos\tNrOfConvMatchVar\tNrOfConvNoMatchVar\tNrOfNoConvPos\tNrOfNoConvVar\t"
       . "ConcorVar\tConcorNoVar\tDiscorVar\tDiscorNoVar\n";

for (my $WindowStart = 1; $WindowStart < $ChromLengths{$QueryChrom}; $WindowStart += $QueryStepSize) {
    print "$QueryChrom\t$WindowStart\n";
    my $WindowEnd = ($WindowStart + $QueryWindowSize - 1 < $ChromLengths{$QueryChrom})
                    ? $WindowStart + $QueryWindowSize - 1
                    : $ChromLengths{$QueryChrom};
    my $BedStart  = $WindowStart - 1;

    my ($NrOfVarPosConvPos, $NrOfVarPosNoConvPos, $NrOfNoVarPosConvPos, $NrOfNoVarPosNoConvPos)
        = (0, 0, 0, 0);
    my ($NrOfVarMatchConv, $NrOfVarMatchNoConv, $NrOfNoVarMatchConv, $NrOfNoVarMatchNoConv)
        = (0, 0, 0, 0);
    my ($NrOfConvPos, $NrOfConvMatchVar, $NrOfConvNoMatchVar, $NrOfNoConvPos, $NrOfNoConvVar)
        = (0, 0, 0, 0, 0);
    my ($ConcorVar, $ConcorNoVar, $DiscorVar, $DiscorNoVar) = (0, 0, 0, 0);

    my (%MappablePositions, %HotspotPositions, %Variants, %PotGcPositions, %PotGcChanges, %AllPotGcChanges);

    foreach my $line (`tabix $MappableBases $QueryChrom:$BedStart-$WindowEnd`) {
        chomp $line;
        my ($c, $s, $e) = split(/\t/, $line);
        for (my $p = $s + 1; $p <= $e; $p++) { $MappablePositions{$p} = undef; }
    }

    # Hotspot-only restriction: load hotspot intervals overlapping the window
    foreach my $line (`tabix $HotspotBED $QueryChrom:$BedStart-$WindowEnd`) {
        chomp $line;
        my ($c, $s, $e) = split(/\t/, $line);
        for (my $p = $s + 1; $p <= $e; $p++) { $HotspotPositions{$p} = undef; }
    }

    foreach my $line (`tabix $VariantFilePath $QueryChrom:$BedStart-$WindowEnd`) {
        chomp $line;
        my @F = split(/\t/, $line);
        my $Variant = $F[3];
        my @V = split(/_/, $Variant);
        my $VarPos = $V[1];
        my ($RefNuc, $AltNuc) = split(/\//, $V[3]);
        next unless ($RefNuc =~ /^[ACGT]$/ && $AltNuc =~ /^[ACGT]$/);

        if ($QueryVarMaxFreq eq "all") {
            $Variants{$VarPos}{$RefNuc}{$AltNuc} = undef;
        } else {
            my $AF = $F[$AfColumn];
            $Variants{$VarPos}{$RefNuc}{$AltNuc} = undef
                if ($AF >= $QueryVarMinFreq && $AF < $QueryVarMaxFreq);
        }
    }

    foreach my $line (`tabix $PositionsGeneConvFilePath $QueryChrom:$BedStart-$WindowEnd`) {
        chomp $line;
        my (undef, undef, $Position) = split(/\t/, $line);
        $PotGcPositions{$Position} = undef;
    }
    foreach my $line (`tabix $ChangesGeneConvFilePath $QueryChrom:$BedStart-$WindowEnd`) {
        chomp $line;
        my (undef, undef, $Position, $Info) = split(/\t/, $line);
        my ($Change) = split(/;/, $Info);
        my ($RefNuc, $AltNuc) = split(/>/, $Change);
        $PotGcChanges{$Position}{$RefNuc}{$AltNuc} = 1;
    }
    foreach my $line (`tabix $AllChangesGeneConvFilePath $QueryChrom:$BedStart-$WindowEnd`) {
        chomp $line;
        my (undef, undef, $Position, $Info) = split(/\t/, $line);
        my ($Change) = split(/;/, $Info);
        my ($RefNuc, $AltNuc) = split(/>/, $Change);
        $AllPotGcChanges{$Position}{$RefNuc}{$AltNuc} = 1;
    }

    for (my $Position = $WindowStart; $Position <= $WindowStart + $QueryWindowSize - 1; $Position++) {
        next unless exists $MappablePositions{$Position};
        next unless exists $HotspotPositions{$Position};   # HOTSPOT RESTRICTION

        if (exists $Variants{$Position}) {
            exists $PotGcPositions{$Position} ? $NrOfVarPosConvPos++   : $NrOfVarPosNoConvPos++;
        } else {
            exists $PotGcPositions{$Position} ? $NrOfNoVarPosConvPos++ : $NrOfNoVarPosNoConvPos++;
        }

        if (exists $Variants{$Position}) {
            foreach my $RefNuc (keys %{$Variants{$Position}}) {
                foreach my $Base (@Bases) {
                    next if $Base eq $RefNuc;
                    # 3-level chained exists check to AVOID Perl autovivification of
                    # PotGcChanges{$Position} when the position is not present.
                    my $is_gc = exists $PotGcChanges{$Position}
                             && exists $PotGcChanges{$Position}{$RefNuc}
                             && exists $PotGcChanges{$Position}{$RefNuc}{$Base};
                    if (exists $Variants{$Position}{$RefNuc}{$Base}) {
                        $is_gc ? $NrOfVarMatchConv++ : $NrOfVarMatchNoConv++;
                    } else {
                        $is_gc ? $NrOfNoVarMatchConv++ : $NrOfNoVarMatchNoConv++;
                    }
                }
            }
        } elsif (exists $PotGcChanges{$Position}) {
            foreach my $RefNuc (keys %{$PotGcChanges{$Position}}) {
                foreach my $Base (@Bases) {
                    next if $Base eq $RefNuc;
                    my $is_gc = exists $PotGcChanges{$Position}{$RefNuc}{$Base};
                    $is_gc ? $NrOfNoVarMatchConv++ : $NrOfNoVarMatchNoConv++;
                }
            }
        } else {
            $NrOfNoVarMatchNoConv += 3;
        }

        if (exists $PotGcChanges{$Position}) {
            foreach my $RefNuc (keys %{$PotGcChanges{$Position}}) {
                foreach my $Base (@Bases) {
                    next if $Base eq $RefNuc;
                    my $is_donor   = exists $PotGcChanges{$Position}{$RefNuc}{$Base};
                    my $is_in_full = exists $AllPotGcChanges{$Position}{$RefNuc}{$Base};
                    if ($is_donor) {
                        my $var_present = exists $Variants{$Position}
                                       && exists $Variants{$Position}{$RefNuc}
                                       && exists $Variants{$Position}{$RefNuc}{$Base};
                        $var_present ? $ConcorVar++ : $ConcorNoVar++;
                    } elsif (!$is_in_full) {
                        my $var_present = exists $Variants{$Position}
                                       && exists $Variants{$Position}{$RefNuc}
                                       && exists $Variants{$Position}{$RefNuc}{$Base};
                        $var_present ? $DiscorVar++ : $DiscorNoVar++;
                    }
                }
            }
        }

        if (exists $PotGcChanges{$Position}) {
            $NrOfConvPos++;
            if (exists $Variants{$Position}) {
                foreach my $RefNuc (keys %{$Variants{$Position}}) {
                    foreach my $AltBase (keys %{$Variants{$Position}{$RefNuc}}) {
                        if (exists $PotGcChanges{$Position}{$RefNuc}{$AltBase}) {
                            $NrOfConvMatchVar++;
                        } else {
                            $NrOfConvNoMatchVar++;
                        }
                    }
                }
            }
        } else {
            $NrOfNoConvPos++;
            if (exists $Variants{$Position}) {
                foreach my $RefNuc (keys %{$Variants{$Position}}) {
                    foreach my $AltBase (keys %{$Variants{$Position}{$RefNuc}}) {
                        $NrOfNoConvVar++;
                    }
                }
            }
        }
    }

    print $O join("\t",
        $QueryChrom, $WindowStart, $WindowEnd,
        $NrOfVarPosConvPos, $NrOfVarPosNoConvPos, $NrOfNoVarPosConvPos, $NrOfNoVarPosNoConvPos,
        $NrOfVarMatchConv, $NrOfVarMatchNoConv, $NrOfNoVarMatchConv, $NrOfNoVarMatchNoConv,
        $NrOfConvPos, $NrOfConvMatchVar, $NrOfConvNoMatchVar, $NrOfNoConvPos, $NrOfNoConvVar,
        $ConcorVar, $ConcorNoVar, $DiscorVar, $DiscorNoVar,
    ) . "\n";
}
close $O;

print "Done. Wrote $OutputFilePath\n";
