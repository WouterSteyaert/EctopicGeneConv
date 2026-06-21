#!/usr/bin/env perl
#===============================================================================
# 2_FilterDNMs.pl  (worker; launched per sublist by the LSF jobs written by
#                   1_CreateListsOfDenovoVcfs.pl)
#
# For one VCF-path sublist:
#   - parse the per-trio denovo VCF;
#   - drop variants whose DE_NOVO_FLAG contains base_fail / altreadparent /
#     abratio (Genomics England's stage-5 DNM-call quality flags);
#   - drop indels and MNVs (require A/C/G/T single nucleotides for both REF
#     and ALT);
#   - emit one row per surviving variant with the per-member GT string.
#
# Inputs:
#   --SubListNr=<N>                            sublist index
#   --ConfigFile=.../config.GRCh38.ini         standard config (added by ProjectConfig)
#
# Outputs:
#   denovo_work_dir/ListsOfFilteredDnms/FilteredDnms_<N>.txt
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use IO::Zlib;
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $SubListNr = 1;
GetOptions("SubListNr=i" => \$SubListNr);

my $WorkDir              = $cfg{paths}{denovo_work_dir};
my $LIST_OF_VCF_FILEPATHs = "$WorkDir/ListsOfDenovoVcfs/DenovoVcfs_$SubListNr.txt";
my $LIST_OF_FILTERED_DNMs = "$WorkDir/ListsOfFilteredDnms/FilteredDnms_$SubListNr.txt";

open my $O, '>', $LIST_OF_FILTERED_DNMs or die "Can't open $LIST_OF_FILTERED_DNMs: $!\n";
open my $L, '<', $LIST_OF_VCF_FILEPATHs or die "Can't open $LIST_OF_VCF_FILEPATHs: $!\n";

while (<$L>) {
    chomp;
    my ($FlaggedDnmFilePath, $FamilyString) = split(/\t/, $_);
    my @FamilyMems = split(/~/, $FamilyString);
    my $TrioId = "";
    my %Family = ();

    foreach my $FamMem (@FamilyMems) {
        my ($T, $Member, $PlateKey) = split(/:/, $FamMem);
        $TrioId = $T;
        $Family{$PlateKey} = $Member;
    }

    my $FH           = new IO::Zlib;
    my %TrioColumns  = ();
    my %PlateKeys    = ();

    if ($FH->open($FlaggedDnmFilePath, "rb")) {
        while (<$FH>) {
            chomp;
            next if /^##/;
            my @LineValues = split(/\t/, $_);

            if (/^#/) {
                for (my $I = 9; $I < scalar @LineValues; $I++) {
                    if (exists $Family{$LineValues[$I]}) {
                        $TrioColumns{$I} = $Family{$LineValues[$I]};
                        $PlateKeys{$I}   = $LineValues[$I];
                    }
                }
                next;
            }

            my ($Chrom, $Position, undef, $RefNuc, $AltNuc,
                $Qual, $Filter, $InfoField, $FormatField, @FamMembers)
                = split(/\t/, $_);

            my %DE_NOVO_FLAGs = ();
            my $Selected      = 1;

            for (my $I = 9; $I < scalar @LineValues; $I++) {
                if (exists $PlateKeys{$I} && exists $Family{$PlateKeys{$I}}) {
                    $DE_NOVO_FLAGs{ (split(/:/, $LineValues[$I]))[-1] } = undef;
                }
            }

            foreach my $flag (keys %DE_NOVO_FLAGs) {
                if ($flag =~ /base_fail/ || $flag =~ /altreadparent/ || $flag =~ /abratio/) {
                    $Selected = 0;
                }
            }

            next unless $Selected;
            next unless ($RefNuc =~ /^[ACGT]$/ && $AltNuc =~ /^[ACGT]$/);

            my @Format = split(/:/, $FormatField);
            my %Format = ();
            for (my $I = 0; $I < scalar @Format; $I++) { $Format{$Format[$I]} = $I; }

            my $GenotypeString = "";
            for (my $I = 0; $I < scalar @LineValues; $I++) {
                if (exists $TrioColumns{$I}) {
                    my $Member       = $TrioColumns{$I};
                    my $PlateKey     = $PlateKeys{$I};
                    my @FormatValues = split(/:/, $LineValues[$I]);
                    $GenotypeString .= "$PlateKey:$Member:$FormatValues[$Format{GT}]~";
                }
            }
            $GenotypeString =~ s/~$//;

            print $O join("\t",
                          $FamilyString,
                          "${Chrom}_${Position}_${Position}_${RefNuc}/${AltNuc}",
                          $Filter, $Qual, $GenotypeString) . "\n";
        }
    }
    $FH->close();
}
close $L;
close $O;

print "Wrote $LIST_OF_FILTERED_DNMs\n";
