#!/usr/bin/env perl
#===============================================================================
# 5_CheckAndSummarizeUnfazed.pl
#
# Aggregate Unfazed per-trio outputs into per-chromosome FATHER_ORIGIN /
# MOTHER_ORIGIN / ALL_ORIGIN BEDs.  For each surviving proband:
#   - check whether Unfazed reported "No phaseable variants" in any of the
#     (re)submitted .e logs; skip if so;
#   - parse the trio's unfazed.vcf;
#     - child GT "1|0" -> Father origin
#     - child GT "0|1" -> Mother origin
#   - aggregate to per-chrom variant tables.
#
# Inputs (paths from config):
#   denovo_cohort_file              Genomics England denovo cohort TSV
#   denovo_work_dir                 top-level DNM work directory
# Reads from <denovo_work_dir>/:
#   ListsOfFilteredDnms/_DnmNumbersPerSample.txt
#   _Jobs/<TrioId>.phasing[.re[.re]].e         (to test "No phaseable variants")
#   phasing/<TrioId>/<TrioId>.unfazed.vcf
#   DenovoVars/ALL/<chr>.bed                    (from 3_ProcessDNMs)
#
# Outputs:
#   DenovoVars/FATHER_ORIGIN/<chr>.bed
#   DenovoVars/MOTHER_ORIGIN/<chr>.bed
#   DenovoVars/ALL_ORIGIN/<chr>.bed             (ALL BED + father/mother fractions)
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $WorkDir          = $cfg{paths}{denovo_work_dir};
my $PhasingDir       = "$WorkDir/phasing/";
my $JobDirectory     = "$WorkDir/_Jobs/";
my $DnmVarsDir       = "$WorkDir/DenovoVars/";
my $NrOfDNMsFp       = "$WorkDir/ListsOfFilteredDnms/_DnmNumbersPerSample.txt";
my $TrioLabKeyFileP  = $cfg{paths}{denovo_cohort_file};

my $MinNrOfDnms = $cfg{denovo}{min_dnms} || 30;
my $MaxNrOfDnms = $cfg{denovo}{max_dnms} || 150;

my %TrioMembers = ();
my %PhasedDnms  = ();
my $NrOfTrios   = 0;

# Map TrioId -> Member -> PlateKey
my $LineNr = 0;
open my $C, '<', $TrioLabKeyFileP or die "Can't open $TrioLabKeyFileP: $!\n";
while (<$C>) {
    if ($LineNr) {
        chomp;
        my @F = split(/\t/, $_);
        if ($F[11] eq "GRCh38") {
            my ($TrioId, $PlateKey, $Member) = ($F[0], $F[3], $F[5]);
            $TrioMembers{$TrioId}{$Member} = $PlateKey;
        }
    }
    $LineNr++;
}
close $C;
$NrOfTrios = scalar keys %TrioMembers;

# Walk valid trios, accumulate phased DNMs
$LineNr = 0;
open my $N, '<', $NrOfDNMsFp or die "Can't open $NrOfDNMsFp: $!\n";
while (<$N>) {
    chomp;
    if ($LineNr) {
        my ($Proband, $TrioId, $NrOfDNMs) = split(/\t/, $_);
        next unless ($NrOfDNMs >= $MinNrOfDnms && $NrOfDNMs <= $MaxNrOfDnms);

        my $Phasable = 1;
        for my $tag ("phasing", "phasing.re", "phasing.re.re") {
            my $eFp = "$JobDirectory$TrioId.$tag.e";
            if (-f $eFp) {
                $Phasable = system("grep -q 'No phaseable variants' $eFp");
                last;
            }
        }
        next unless $Phasable;

        my $UnfazedVcf = "$PhasingDir$TrioId/$TrioId.unfazed.vcf";
        unless (-f $UnfazedVcf) {
            print "$TrioId\tmissing-unfazed-vcf\n";
            next;
        }

        my %SampleToCol = ();
        open my $V, '<', $UnfazedVcf or die "Can't open $UnfazedVcf: $!\n";
        while (<$V>) {
            chomp;
            next if /^##/;
            my @LV = split(/\t/, $_);
            if (/^#/) {
                for (my $I = 9; $I < @LV; $I++) {
                    $SampleToCol{$LV[$I]} = $I;
                }
                next;
            }

            my $ChildCol  = $SampleToCol{$TrioMembers{$TrioId}{"Offspring"}};
            my $ChildGeno = (split(/:/, $LV[$ChildCol]))[0];
            $LV[0] =~ s/^chr//;
            my $Variant = "$LV[0]_$LV[1]_$LV[1]_$LV[3]/$LV[4]";

            if    ($ChildGeno eq "1|0") {
                $PhasedDnms{"Father"}{$LV[0]}{$Variant}++;
            }
            elsif ($ChildGeno eq "0|1") {
                $PhasedDnms{"Mother"}{$LV[0]}{$Variant}++;
            }
        }
        close $V;
    }
    $LineNr++;
}
close $N;

make_path("$DnmVarsDir/FATHER_ORIGIN") unless -d "$DnmVarsDir/FATHER_ORIGIN";
make_path("$DnmVarsDir/MOTHER_ORIGIN") unless -d "$DnmVarsDir/MOTHER_ORIGIN";
make_path("$DnmVarsDir/ALL_ORIGIN")    unless -d "$DnmVarsDir/ALL_ORIGIN";

# Father/Mother per-chrom BEDs
for my $pair (["Father", "FATHER_ORIGIN"], ["Mother", "MOTHER_ORIGIN"]) {
    my ($Key, $Dir) = @$pair;
    foreach my $Chrom (sort keys %{$PhasedDnms{$Key}}) {
        my $OutputFilePath = "$DnmVarsDir$Dir/$Chrom.bed";
        open my $O, '>', $OutputFilePath or die "Can't open $OutputFilePath: $!\n";
        foreach my $Variant (sort keys %{$PhasedDnms{$Key}{$Chrom}}) {
            my (undef, $Position) = split(/_/, $Variant);
            my $BedStart   = $Position - 1;
            my $Recurrence = $PhasedDnms{$Key}{$Chrom}{$Variant};
            my $Frequency  = $Recurrence / $NrOfTrios;
            print $O "$Chrom\t$BedStart\t$Position\t$Variant\t$Frequency\t$Recurrence\t$NrOfTrios\n";
        }
        close $O;
    }
}

# ALL_ORIGIN BED: ALL per-chr BED augmented with father/mother counts+fractions
opendir(my $DH, "$DnmVarsDir/ALL") or die "Can't open $DnmVarsDir/ALL: $!\n";
while (my $File = readdir($DH)) {
    next unless ($File =~ /\.bed$/ && $File !~ /gnomad/);
    my $InFp  = "$DnmVarsDir/ALL/$File";
    my $OutFp = "$DnmVarsDir/ALL_ORIGIN/$File";
    open my $O, '>', $OutFp or die "Can't open $OutFp: $!\n";
    open my $F, '<', $InFp  or die "Can't open $InFp: $!\n";
    while (<$F>) {
        chomp;
        print $O $_;
        my @LV       = split(/\t/, $_);
        my $Chrom    = $LV[0];
        my $Variant  = $LV[3];
        my $TotCount = $LV[5];

        if (my $f = $PhasedDnms{"Father"}{$Chrom}{$Variant}) {
            print $O "\t$f\t" . ($f / $TotCount);
        } else {
            print $O "\t0\tNA";
        }
        if (my $m = $PhasedDnms{"Mother"}{$Chrom}{$Variant}) {
            print $O "\t$m\t" . ($m / $TotCount);
        } else {
            print $O "\t0\tNA";
        }
        print $O "\n";
    }
    close $F;
    close $O;
}
closedir($DH);

print "Done. Aggregated FATHER_ORIGIN, MOTHER_ORIGIN, ALL_ORIGIN under $DnmVarsDir\n";
