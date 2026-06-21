#!/usr/bin/env perl
#===============================================================================
# 4_UpdateDnmsWithGnomAd.pl
#
# Annotate DNM BEDs (per-chr) with gnomAD AF/AC/AN by exact variant match.
# Runs once per origin set (ALL / MALE / FEMALE / FATHER_ORIGIN / MOTHER_ORIGIN
# / ALL_ORIGIN).  Origin selects the input/output subdirectory inside
# denovo_work_dir/DenovoVars/.
#
# gnomAD input is the per-chromosome sorted BED produced by step 02:
#   gnomad_genome_dir/<chr>.sorted.bed.gz
#
# Usage:
#   perl 4_UpdateDnmsWithGnomAd.pl \
#       --Origin=MOTHER_ORIGIN \
#       --ConfigFile=../00_Configuration/config.GRCh38.ini
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use IO::Zlib;
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $Origin = "";
GetOptions("Origin=s" => \$Origin);
die "--Origin=<ALL|MALE|FEMALE|FATHER_ORIGIN|MOTHER_ORIGIN|ALL_ORIGIN> required\n"
    unless $Origin =~ /^(ALL|MALE|FEMALE|FATHER_ORIGIN|MOTHER_ORIGIN|ALL_ORIGIN)$/;

my $DnmVarDirectory  = "$cfg{paths}{denovo_work_dir}/DenovoVars/$Origin/";
my $GnomAdDirectory  = $cfg{paths}{gnomad_genome_dir};

for (my $I = 1; $I <= 24; $I++) {
    my $Chrom = $I;
    $Chrom = "X" if $I == 23;
    $Chrom = "Y" if $I == 24;

    print "$Chrom\n";

    my $DnmFilePath    = "${DnmVarDirectory}${Chrom}.bed";
    my $GnomAdFilePath = "${GnomAdDirectory}/${Chrom}.sorted.bed.gz";
    my $OutFilePath    = "${DnmVarDirectory}${Chrom}.gnomad.bed";

    unless (-f $DnmFilePath) {
        print "\tSkip: $DnmFilePath not present\n";
        next;
    }

    my %DNMsInfo   = ();
    my %GnomAdInfo = ();

    print "\tRead In DNMs\n";
    open my $D, '<', $DnmFilePath or die "Can't open $DnmFilePath: $!\n";
    while (<$D>) {
        chomp;
        my @F = split(/\t/, $_);
        $DNMsInfo{$F[3]} = $_;
    }
    close $D;

    print "\tIterate Over GnomAD\n";
    my $FH = new IO::Zlib;
    if ($FH->open($GnomAdFilePath, "rb")) {
        while (<$FH>) {
            chomp;
            my @F          = split(/\t/, $_);
            my $Variant    = $F[3];
            my $GnomAdFreq = "$F[4]\t$F[5]\t$F[6]";
            $GnomAdInfo{$Variant} = $GnomAdFreq if exists $DNMsInfo{$Variant};
        }
    } else {
        warn "\tCan't open $GnomAdFilePath\n";
    }
    $FH->close();

    print "\tWrite Out\n";
    open my $O, '>', $OutFilePath or die "Can't open $OutFilePath: $!\n";
    foreach my $Dnm (sort keys %DNMsInfo) {
        print $O $DNMsInfo{$Dnm};
        print $O "\t" . $GnomAdInfo{$Dnm} if exists $GnomAdInfo{$Dnm};
        print $O "\n";
    }
    close $O;
}
