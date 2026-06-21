#!/usr/bin/env perl
#===============================================================================
# __1_ComputeContingencyStatistics.pl
#
# Worker: compute per-window contingency-table counts for one combination of
# chromosome × template length (k) × AF bin × mappability context.
#
# For each mappable position inside the window, the worker emits counts along
# two axes:
#   (1) position-based  — whether the position is gene conversion-compatible
#                         at this k, and whether it harbours a variant in the
#                         AF bin;
#   (2) alteration-based — for each nucleotide change at each position, whether
#                          the change matches a predicted donor allele;
# plus per-window concordant/discordant variant counts (donor-allele-aware).
#
# Mappability context is selected via --MappabilityContext:
#   allmapp        all-mappable (default; GRCh38_notinlowmappabilityall)
#   segdupmapp     segmental-duplications subset
#   nosegdupmapp   mappable outside segmental duplications
#
# Variant set is selected via --QueryVarList (relative path under data_root)
# and --AfColumn (0-based index of the gnomAD-AF column in the variant BED):
#   gnomAD genome variants  --QueryVarList=gnomad/genome      --AfColumn=4
#   GWAS clean variants     --QueryVarList=gwas/clean/<class> --AfColumn=9
#   DNM-with-gnomAD beds    --QueryVarList=dnm/DenovoVars/<o> --AfColumn=<idx>
#   Random control          --QueryVarList=random3            --AfColumn=4
#     (random has no AF; pass --QueryVarMaxFreq=all)
#
# Outputs one TSV row per window to:
#   <stats_dir>/stats_<context>/<k>/<chr>_<window>_<step>_<varlist>_<minF>_<maxF>.<start>.<end>.txt
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
);

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
my $OutputFilePath  = "$StatsSubdir/${QueryChrom}_${QueryWindowSize}_${QueryStepSize}_${QueryVarListPrint}_${MinPr}_${MaxPr}.${QueryIntervalStart}.${QueryIntervalEnd}.txt";

# Positions used in this k:
#   k == 91             -> conv/91 (no longer k to set-difference against)
#   k in {17,19}        -> gcdiffs/k_(k+2)
#   else                -> gcdiffs/k_(k+10)
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

# Read chromosome length for the requested chromosome
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

    my (%MappablePositions, %Variants, %PotGcPositions, %PotGcChanges, %AllPotGcChanges);

    foreach my $line (`tabix $MappableBases $QueryChrom:$BedStart-$WindowEnd`) {
        chomp $line;
        my ($c, $s, $e) = split(/\t/, $line);
        for (my $p = $s + 1; $p <= $e; $p++) { $MappablePositions{$p} = undef; }
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
        my ($Change, $Recurrence) = split(/;/, $Info);
        my ($RefNuc, $AltNuc)     = split(/>/, $Change);
        $PotGcChanges{$Position}{$RefNuc}{$AltNuc} = $Recurrence;
    }

    foreach my $line (`tabix $AllChangesGeneConvFilePath $QueryChrom:$BedStart-$WindowEnd`) {
        chomp $line;
        my (undef, undef, $Position, $Info) = split(/\t/, $line);
        my ($Change, $Recurrence) = split(/;/, $Info);
        my ($RefNuc, $AltNuc)     = split(/>/, $Change);
        $AllPotGcChanges{$Position}{$RefNuc}{$AltNuc} = $Recurrence;
    }

    for (my $Position = $WindowStart; $Position <= $WindowStart + $QueryWindowSize - 1; $Position++) {
        next unless exists $MappablePositions{$Position};

        # Position-based contingency
        if (exists $Variants{$Position}) {
            exists $PotGcPositions{$Position} ? $NrOfVarPosConvPos++   : $NrOfVarPosNoConvPos++;
        } else {
            exists $PotGcPositions{$Position} ? $NrOfNoVarPosConvPos++ : $NrOfNoVarPosNoConvPos++;
        }

        # Alteration-based contingency
        if (exists $Variants{$Position}) {
            foreach my $RefNuc (keys %{$Variants{$Position}}) {
                foreach my $Base (@Bases) {
                    next if $Base eq $RefNuc;
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

        # Concordant / Discordant
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

        # Conv-binomial
        if (exists $PotGcChanges{$Position}) {
            $NrOfConvPos++;
            if (exists $Variants{$Position}) {
                foreach my $RefNuc (keys %{$Variants{$Position}}) {
                    foreach my $AltBase (keys %{$Variants{$Position}{$RefNuc}}) {
                        if (exists $PotGcChanges{$Position}{$RefNuc}
                         && exists $PotGcChanges{$Position}{$RefNuc}{$AltBase}) {
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
