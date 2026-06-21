#!/usr/bin/env perl
#===============================================================================
# 0_FindTriosToExclude.pl
#
# Build a list of family IDs whose trio members cannot be linked back to the
# rare-disease analysis LabKey table (i.e. families for which at least one
# member is missing from the current LabKey release).  These families are
# excluded from all downstream DNM analyses to avoid mixing trios with
# incomplete phenotypic linkage.
#
# Inputs  (paths read from config.GRCh38.ini, section <paths>):
#   labkey_rd_analysis    LabKey rare-disease analysis TSV (one row per
#                         participant; participant ID in column 1)
#   denovo_cohort_file    Genomics England denovo cohort information TSV
#                         (V9 fix file); FamilyId in col 2 (0-based: 1),
#                         ParticipantId in col 5 (0-based: 4), Assembly in
#                         col 12 (0-based: 11)
#
# Output:
#   denovo_work_dir/FamiliesToExclude.txt   one family ID per line
#
# Run on the Genomics England LSF cluster (login node is fine; this is a
# short read/write job).
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $LabKeyFilePath       = $cfg{paths}{labkey_rd_analysis};
my $DeNovoCohortFilePath = $cfg{paths}{denovo_cohort_file};
my $WorkDir              = $cfg{paths}{denovo_work_dir};
my $FamiliesToExclude    = "$WorkDir/FamiliesToExclude.txt";

make_path($WorkDir) unless -d $WorkDir;

my %RdParticipants = ();
my %DenovoCohort   = ();
my $LineNr         = 0;

# Read in rare-disease participants from current LabKey release
open my $L, '<', $LabKeyFilePath or die "Can't open $LabKeyFilePath: $!\n";
while (<$L>) {
    if ($LineNr) {
        chomp;
        my @F = split(/\t/, $_);
        $RdParticipants{$F[0]} = undef;
    }
    $LineNr++;
}
close $L;

# Read in denovo cohort (GRCh38 only)
$LineNr = 0;
open my $D, '<', $DeNovoCohortFilePath or die "Can't open $DeNovoCohortFilePath: $!\n";
while (<$D>) {
    if ($LineNr) {
        chomp;
        my @F             = split(/\t/, $_);
        my $FamilyId      = $F[1];
        my $ParticipantId = $F[4];
        my $Assembly      = $F[11];
        if ($Assembly eq "GRCh38") {
            $DenovoCohort{$FamilyId}{$ParticipantId} = undef;
        }
    }
    $LineNr++;
}
close $D;

# Emit families with any unlinked member
open my $F, '>', $FamiliesToExclude or die "Can't open $FamiliesToExclude: $!\n";
foreach my $FamilyId (sort keys %DenovoCohort) {
    foreach my $ParticipantId (keys %{$DenovoCohort{$FamilyId}}) {
        if (!exists $RdParticipants{$ParticipantId}) {
            print $F $FamilyId . "\n";
            last;
        }
    }
}
close $F;

print "Wrote families-to-exclude list: $FamiliesToExclude\n";
