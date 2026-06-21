#!/usr/bin/env perl
#===============================================================================
# __5_FetchGwasVariants.pl  (worker)
#
# For one GWAS Catalog sub-file (--GwasFilePath=*.NN.tsv produced by the
# submitter), resolve each "strongest SNP" rsID to GRCh38 coordinates and
# alleles via the Ensembl variation API, classify the catalog risk allele as
# Risk (risk allele != reference) or Protection (risk allele == reference,
# i.e. each non-risk allele in the variant becomes a protective variant), and
# emit a per-variant BED-like row:
#
#   chr  bedStart  end  variant(chr_start_end_REF/ALT)  Pvalue  OrBeta  Disease  MappedTrait  Type
#
# Coordinate match against the catalog's reported chr+pos serves as a sanity
# check (deviations are logged).
#
# Inputs (--GwasFilePath = one of the *.NN.tsv split-files written by the
# submitter):
#   each line is one catalog association; columns 12, 13, 21, 28, 31, 35 hold
#   chr-list, pos-list, strongest-SNP-list, p-value, OR/beta, mapped trait
#   (0-based: 11, 12, 20, 27, 30, 34).
#
# Output:
#   <GwasFilePath>.processed.tsv
#
# Requires the Ensembl Perl API to be in PERL5LIB (it ships with the public
# Ensembl release).  Connection details (release, host) are read from config.
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $GwasFilePath = "";
GetOptions("GwasFilePath=s" => \$GwasFilePath);
die "--GwasFilePath=<one sub-file> required\n" unless $GwasFilePath;

my $EnsemblRelease = $cfg{gwas}{ensembl_release} || 113;
my $OutputFilePath = $GwasFilePath;
$OutputFilePath =~ s/\.tsv$/.processed.tsv/;

# Connect to the public Ensembl MySQL server for the configured release.
require Bio::EnsEMBL::Registry;
my $Registry = "Bio::EnsEMBL::Registry";
$Registry->load_registry_from_db(
    -host    => "ensembldb.ensembl.org",
    -user    => "anonymous",
    -port    => 3306,
    -db_version => $EnsemblRelease,
);

my $SliceAdaptor     = $Registry->get_adaptor('Human', 'Core',      'Slice');
my $VariationAdaptor = $Registry->get_adaptor('Human', 'Variation', 'Variation');

open my $O, '>', $OutputFilePath or die "Can't open $OutputFilePath: $!\n";
open my $G, '<', $GwasFilePath   or die "Can't open $GwasFilePath: $!\n";

my $LineNr = 0;
while (<$G>) {
    if ($LineNr) {
        chomp;
        my @F = split(/\t/, $_);

        my $Disease     = $F[7];
        my $TableChroms = $F[11];   $TableChroms =~ s/ //g;
        my @TableChroms = split(/;/, $TableChroms);
        my $TablePos    = $F[12];   $TablePos    =~ s/ //g;
        my @TablePos    = split(/;/, $TablePos);
        my $StrongestSnps = $F[20];
        my @StrongestSnps = split(/;/, $StrongestSnps);
        my $Pvalue      = $F[27];
        my $OrBeta      = $F[30];
        my $MappedTrait = $F[34];

        next unless ($TableChroms && $TablePos);

        for (my $S = 0; $S < @StrongestSnps; $S++) {
            my $StrongestSnp  = $StrongestSnps[$S];
            my $TableChrom    = $TableChroms[$S];
            my $TablePosition = $TablePos[$S];
            next unless ($StrongestSnp && $TableChrom && $TablePosition);

            for ($StrongestSnp, $TableChrom, $TablePosition) { s/ //g; }

            if ($StrongestSnp =~ /^(rs\d+)-([ACGT])$/) {
                my ($RsId, $RiskAllele) = ($1, $2);
                my $VariationObject = $VariationAdaptor->fetch_by_name($RsId);
                next unless $VariationObject;

                my %Variation = ();
                foreach my $Feature (@{ $VariationObject->get_all_VariationFeatures() }) {
                    my $Chrom        = $Feature->seq_region_name;
                    my $StartPos     = $Feature->start;
                    my $EndPos       = $Feature->end;
                    my $AlleleString = $Feature->allele_string;
                    next unless ($Chrom =~ /^\d+$/ || $Chrom eq "X" || $Chrom eq "Y");

                    if ($Chrom eq $TableChrom && $StartPos eq $TablePosition && $EndPos eq $TablePosition) {
                        $Variation{"${Chrom}_${StartPos}_${AlleleString}"} = undef;
                    } else {
                        print "ERROR: non-matching variant ($StartPos / $EndPos vs $TablePosition)\n";
                    }
                }

                next unless (scalar keys %Variation == 1);

                foreach my $Variant (keys %Variation) {
                    my ($Chrom, $Position, $AlleleString) = split(/_/, $Variant);
                    my @AlleleString = split(/\//, $AlleleString);
                    my %AlleleStringH = map { $_ => 1 } @AlleleString;
                    my $RefAllele = $SliceAdaptor
                                    ->fetch_by_region('Chromosome', $Chrom, $Position, $Position)
                                    ->seq();

                    unless (exists $AlleleStringH{$RefAllele}) {
                        print "ERROR: ref allele $RefAllele not in catalog allele string $AlleleString ($RsId)\n";
                        next;
                    }

                    my $BedPos = $Position - 1;

                    if ($RefAllele ne $RiskAllele) {
                        my $Variation_id = "${Chrom}_${Position}_${Position}_${RefAllele}/${RiskAllele}";
                        print $O join("\t",
                            $Chrom, $BedPos, $Position, $Variation_id,
                            $Pvalue, $OrBeta, $Disease, $MappedTrait, "Risk") . "\n";
                    } else {
                        # ref allele equals risk allele -> each remaining catalog
                        # allele becomes a Protection variant
                        foreach my $Allele (@AlleleString) {
                            next if $Allele eq $RefAllele;
                            my $Variation_id = "${Chrom}_${Position}_${Position}_${RefAllele}/${Allele}";
                            print $O join("\t",
                                $Chrom, $BedPos, $Position, $Variation_id,
                                $Pvalue, $OrBeta, $Disease, $MappedTrait, "Protection") . "\n";
                        }
                    }
                }
            }
            elsif ($StrongestSnp !~ /^rs\d+-\?$/) {
                print "ERROR: format of strongest SNP ($StrongestSnp)\n";
            }
        }
    }
    $LineNr++;
}
close $G;
close $O;

print "Wrote $OutputFilePath\n";
