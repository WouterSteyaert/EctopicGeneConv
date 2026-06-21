#!/usr/bin/env perl
#===============================================================================
# 3_ProcessDNMs.pl
#
# Aggregate all per-sublist FilteredDnms_<N>.txt produced by step 2 into the
# final per-chromosome, per-sex DNM BED set used downstream.
#
# Filters applied here:
#   - offspring genotype must be heterozygous (0/1, 1/0, 0|1, 1|0);
#   - per-trio total DNM count must be between min_dnms and max_dnms (config)
#     - this is the cohort quality filter described in the paper.
#
# Outputs (under denovo_work_dir/):
#   ListsOfFilteredDnms/_DnmFrequencies.txt        Variant -> recurrence
#   ListsOfFilteredDnms/_DnmSamples.txt            TrioId / Proband mapping
#   ListsOfFilteredDnms/_DnmNumbersPerSample.txt   per-proband nr of DNMs
#   ListsOfFilteredDnms/_FilteredDnms.all.bed      flat BED of all DNMs
#   ListsOfFilteredDnms/_FilteredDnms.sorted.all.bed{,.gz,.gz.tbi}
#   DenovoVars/{ALL,MALE,FEMALE}/<chr>.bed         per-sex per-chr BED
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $LabKeyFileP        = $cfg{paths}{denovo_cohort_file};
my $WorkDir            = $cfg{paths}{denovo_work_dir};
my $DnmDirectory       = "$WorkDir/ListsOfFilteredDnms/";
my $DnmVarDirectory    = "$WorkDir/DenovoVars/";
my $FamiliesToExclude  = "$WorkDir/FamiliesToExclude.txt";
my $DnmFrequencies     = "${DnmDirectory}_DnmFrequencies.txt";
my $DnmSamples         = "${DnmDirectory}_DnmSamples.txt";
my $DnmNumbers         = "${DnmDirectory}_DnmNumbersPerSample.txt";
my $OutputBedFp        = "${DnmDirectory}_FilteredDnms.all.bed";
my $SortedOutputBedFp  = "${DnmDirectory}_FilteredDnms.sorted.all.bed";

my $MinNrOfDnms = $cfg{denovo}{min_dnms} || 30;
my $MaxNrOfDnms = $cfg{denovo}{max_dnms} || 150;

make_path($DnmVarDirectory) unless -d $DnmVarDirectory;

my %FamiliesToExcludeH    = ();
my %SampleMappings        = ();
my %ProbandToSex          = ();
my %DNMsToProband         = ();
my %DNMsPerChrom          = ();
my %ChildToDnmsRaw        = ();
my %ProbandsToDnmsClean   = ();
my %ProbandsPerSex        = ();
my %NrOfProbandsPerSex    = ();
my %Probands              = ();
my @Sexes                 = ("ALL", "MALE", "FEMALE");

open my $T, '<', $FamiliesToExclude or die "Can't open $FamiliesToExclude: $!\n";
while (<$T>) { chomp; $FamiliesToExcludeH{$_} = undef; }
close $T;

my $LineNr = 0;
open my $C, '<', $LabKeyFileP or die "Can't open $LabKeyFileP: $!\n";
while (<$C>) {
    if ($LineNr) {
        chomp;
        my @F = split(/\t/, $_);
        my ($TrioId, $FamilyId, $PlateKey, $Member, $Sex, $Assembly)
            = ($F[0], $F[1], $F[3], $F[5], $F[10], $F[11]);
        if ($Member eq "Offspring" && !exists $FamiliesToExcludeH{$FamilyId}) {
            $SampleMappings{$PlateKey} = $TrioId;
            $ProbandToSex{$PlateKey}   = $Sex;
        }
    }
    $LineNr++;
}
close $C;

# Pass 1: raw per-proband DNM counts
opendir(my $DH, $DnmDirectory) or die "Can't open $DnmDirectory: $!\n";
while (my $DnmFile = readdir($DH)) {
    next unless $DnmFile =~ /^FilteredDnms_\d+\.txt$/;
    my $fp = "$DnmDirectory$DnmFile";
    open my $D, '<', $fp or die "Can't open $fp: $!\n";
    while (<$D>) {
        chomp;
        my ($FamilyString, $Variant, $Filter, $Qual, $GenotypeString) = split(/\t/, $_);
        foreach my $gv (split(/~/, $GenotypeString)) {
            my ($PlateKey, $Member, $Genotype) = split(/:/, $gv);
            next unless $Member eq "Offspring";
            next unless ($Genotype eq "0/1" || $Genotype eq "1/0"
                      || $Genotype eq "1|0" || $Genotype eq "0|1");
            $ChildToDnmsRaw{$PlateKey}{$Variant} = undef;
        }
    }
    close $D;
}
closedir($DH);

foreach my $PlateKey (sort keys %ChildToDnmsRaw) {
    my $Nr = scalar keys %{$ChildToDnmsRaw{$PlateKey}};
    print "$PlateKey\t$Nr\n";
}

# Pass 2: keep only DNMs from probands inside the 30-150 window
opendir($DH, $DnmDirectory) or die "Can't open $DnmDirectory: $!\n";
while (my $DnmFile = readdir($DH)) {
    next unless $DnmFile =~ /^FilteredDnms_\d+\.txt$/;
    my $fp = "$DnmDirectory$DnmFile";
    open my $D, '<', $fp or die "Can't open $fp: $!\n";
    while (<$D>) {
        chomp;
        my ($FamilyString, $Variant, $Filter, $Qual, $GenotypeString) = split(/\t/, $_);
        my ($Chrom) = split(/_/, $Variant);
        foreach my $gv (split(/~/, $GenotypeString)) {
            my ($PlateKey, $Member, $Genotype) = split(/:/, $gv);
            next unless $Member eq "Offspring";
            next unless ($Genotype eq "0/1" || $Genotype eq "1/0"
                      || $Genotype eq "1|0" || $Genotype eq "0|1");
            my $N = scalar keys %{$ChildToDnmsRaw{$PlateKey}};
            next unless ($N >= $MinNrOfDnms && $N <= $MaxNrOfDnms);

            $DNMsToProband{$ProbandToSex{$PlateKey}}{$Variant}{$PlateKey} = undef;
            $DNMsToProband{"ALL"}{$Variant}{$PlateKey}                    = undef;
            $Probands{$PlateKey}                                          = undef;
            $DNMsPerChrom{$Chrom}{$Variant}                               = undef;
            $ProbandsToDnmsClean{$PlateKey}{$Variant}                     = undef;
            $ProbandsPerSex{$ProbandToSex{$PlateKey}}{$PlateKey}          = undef;
        }
    }
    close $D;
}
closedir($DH);

$NrOfProbandsPerSex{"MALE"}   = scalar keys %{$ProbandsPerSex{"MALE"}};
$NrOfProbandsPerSex{"FEMALE"} = scalar keys %{$ProbandsPerSex{"FEMALE"}};
$NrOfProbandsPerSex{"ALL"}    = $NrOfProbandsPerSex{"MALE"} + $NrOfProbandsPerSex{"FEMALE"};

# Write outputs
open my $B, '>', $OutputBedFp   or die "Can't open $OutputBedFp: $!\n";
open my $O, '>', $DnmFrequencies or die "Can't open $DnmFrequencies: $!\n";
print $O "Variant\tNumberOfObservations\n";
foreach my $Variant (sort keys %{$DNMsToProband{"ALL"}}) {
    my ($Chrom, $Position, undef, $Change) = split(/_/, $Variant);
    my $BedStart = $Position - 1;
    my $NrOfObs  = scalar keys %{$DNMsToProband{"ALL"}{$Variant}};
    my $ProbandList = join('~', sort keys %{$DNMsToProband{"ALL"}{$Variant}});
    print $O "$Variant\t$NrOfObs\n";
    print $B "$Chrom\t$BedStart\t$Position\tChange=$Change;Samples=$ProbandList\n";
}
close $O;
close $B;

system("bedtools sort -i $OutputBedFp > $SortedOutputBedFp");
system("bgzip -c $SortedOutputBedFp > $SortedOutputBedFp.gz");
system("tabix -p bed $SortedOutputBedFp.gz");

open my $S, '>', $DnmSamples or die "Can't open $DnmSamples: $!\n";
print $S "TrioId\tProband\n";
foreach my $Proband (sort keys %Probands) {
    print $S "$SampleMappings{$Proband}\t$Proband\n";
}
close $S;

open my $DN, '>', $DnmNumbers or die "Can't open $DnmNumbers: $!\n";
print $DN "Proband\tTrioId\tNrOfDnms\tSex\n";
foreach my $Proband (sort keys %ProbandsToDnmsClean) {
    my $NrOfDnms = scalar keys %{$ProbandsToDnmsClean{$Proband}};
    print $DN "$Proband\t$SampleMappings{$Proband}\t$NrOfDnms\t$ProbandToSex{$Proband}\n";
}
close $DN;

# Per-sex per-chrom BEDs (note: Freq column is only meaningful for ALL).
foreach my $Sex (@Sexes) {
    my $BedSubdir = "$DnmVarDirectory$Sex/";
    make_path($BedSubdir) unless -d $BedSubdir;

    foreach my $Chrom (sort keys %DNMsPerChrom) {
        my $BedFilePath = "$BedSubdir$Chrom.bed";
        open my $BF, '>', $BedFilePath or die "Can't open $BedFilePath: $!\n";
        foreach my $Variant (sort keys %{$DNMsPerChrom{$Chrom}}) {
            my (undef, $Position, undef) = split(/_/, $Variant);
            my $BedStart = $Position - 1;
            if (exists $DNMsToProband{$Sex}{$Variant}) {
                my $NrOfObs = scalar keys %{$DNMsToProband{$Sex}{$Variant}};
                my $Freq    = $NrOfObs / $NrOfProbandsPerSex{$Sex};
                print $BF "$Chrom\t$BedStart\t$Position\t$Variant\t$Freq\t$NrOfObs\t$NrOfProbandsPerSex{$Sex}\n";
            }
        }
        close $BF;
    }
}

print "Done. Per-sex per-chrom BEDs written under $DnmVarDirectory\n";
