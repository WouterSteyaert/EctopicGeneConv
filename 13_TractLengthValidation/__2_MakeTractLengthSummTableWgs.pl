#!/usr/bin/env perl
#===============================================================================
# __2_MakeTractLengthSummTableWgs.pl
#
# Aggregate per-(window, alignment) outputs of __2_ComputeTractLengthWgs into
# a single summary table and a per-variant "max tract analysis length" scatter
# table (Variant Frequency TractLength RepeatLength).
#
# For each variant, the maximum tract analysis length (largest k for which the
# variant is centred within the homologous tract) is retained across ALL BLAST
# alignments AND all window sizes — this is the value used in the paper's
# tract-length validation figures and tables.
#
# Inputs (config):
#   <paths>tractlengths_wgs_dir   per-(variant, window) per-alignment .txt or .gz
#   <paths>tractlengths_sum_dir   contains Variant.<seed>.input listing all
#                                 variants targeted by the corresponding run
#===============================================================================
use strict;
use warnings;
use List::Util qw(min max);
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $TractLengthDir       = "$cfg{paths}{tractlengths_wgs_dir}/";
my $TractLengthSumDir    = "$cfg{paths}{tractlengths_sum_dir}/";
my $TractLengthSummTable = "${TractLengthSumDir}_TractLengthSummaryTable.sum";

make_path($TractLengthSumDir) unless -d $TractLengthSumDir;

my $UsedSeed   = 0;
my %InputVars  = ();
my %AllTract   = ();
my %AllRepeat  = ();

opendir(my $DIR, $TractLengthSumDir) or die "Can't open $TractLengthSumDir: $!\n";
while (my $File = readdir($DIR)) {
    next unless $File =~ /^Variant\.(\d+)\.input$/;
    $UsedSeed = $1;
    open my $F, '<', "$TractLengthSumDir$File" or die "Can't open $File: $!\n";
    while (<$F>) {
        chomp;
        my ($Variant, $Frequency) = split(/\t/, $_);
        $InputVars{$Variant} = $Frequency;
    }
    close $F;
}
closedir($DIR);

my $TractLengthScatter = "${TractLengthSumDir}_TractLengthsVsFreq.${UsedSeed}.scatter";

print "Cat Tract Length Tables\n";
chdir($TractLengthDir);
my $Header = "QUERY_VARIANT\tFREQ\tChrom\tAcceptorStart\tAcceptorEnd\tChromName\t"
           . "FirstSubjectCoord\tLastSubjectCoord\tNrOfMatches\tSUNs\tHomSUNs\t"
           . "TractLengthRange\tStrand\tLargestSmallerHom\tSmallestLargerHom\t"
           . "LargestSmallerSunHom\tSmallestLargerSunHom\tLargestSmallerSun\t"
           . "SmallestLargerSun\tLargestSmaller\tSmallestLarger\tTractLength\t"
           . "TractAnalysisLength\tLeftTractFlank\tRightTractFlank\tFlankRatio\t"
           . "MaxTractLength\tMaxRepeatLength";
system("echo -e \"$Header\" > $TractLengthSummTable");
system("find . -name '*.txt' -exec cat {} \\; >> $TractLengthSummTable");
system("find . -name '*.gz'  -exec zcat {} \\; >> $TractLengthSummTable");

my $LineNr = 0;
open my $T, '<', $TractLengthSummTable or die "Can't open $TractLengthSummTable: $!\n";
while (<$T>) {
    if ($LineNr) {
        chomp;
        my @F = split(/\t/, $_);
        my ($Variant, $Freq, $TractLen, $TractAnalysisLen) = ($F[0], $F[1], $F[21], $F[22]);
        $AllTract{"$Variant~$Freq"}{$TractLen} = undef;
        $AllRepeat{"$Variant~$Freq"}{$TractAnalysisLen} = undef;
    }
    $LineNr++;
}
close $T;

print "Make Table TractLength Vs Freq\n";
open my $O, '>', $TractLengthScatter or die "Can't open $TractLengthScatter: $!\n";
print $O "Variant\tFrequency\tTractLength\tRepeatLength\n";
foreach my $Variant (sort keys %InputVars) {
    my $Frequency = $InputVars{$Variant} // "/";
    if (exists $AllTract{"$Variant~$Frequency"}) {
        my $MaxT = max(keys %{$AllTract{"$Variant~$Frequency"}});
        my $MaxR = max(keys %{$AllRepeat{"$Variant~$Frequency"}});
        print $O "$Variant\t$Frequency\t$MaxT\t$MaxR\n";
    } else {
        print $O "$Variant\t$Frequency\t0\t0\n";
    }
}
close $O;

print "Done. Scatter table: $TractLengthScatter\n";
