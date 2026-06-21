#!/usr/bin/env perl
#===============================================================================
# __4_GenerateRandomVariants.pl  (worker; one job per random set)
#
# Generate one random control set of single-nucleotide variants distributed
# uniformly across chromosomes 1..22, X, Y in proportion to chromosome length.
#
# For each variant:
#   - sample a position uniformly from a chromosome;
#   - fetch the reference base with twoBitToFa;
#   - skip if reference == 'N';
#   - pick an alternative allele uniformly from the three non-reference bases;
#   - emit a BED row identical in format to the per-chromosome gnomAD BEDs:
#       chr  bedStart  bedEnd  chr_start_end_REF/ALT
#
# Per-chrom output is bedtools-sorted, bgzipped, and tabix-indexed.
#
# Inputs (config <paths> / <random>):
#   reference_dict          GRCh38 Picard-style sequence dictionary (for chrom lengths)
#   reference_2bit          GRCh38 2bit file (for twoBitToFa)
#   data_root               base; random sets live under <data_root>/<RandomSet>/
#   n_variants_total        total target genome-wide
#   genome_size_bp          denominator used to scale per-chrom count (~3e9)
#
# Usage:
#   perl __4_GenerateRandomVariants.pl --RandomSet=random1
#     [ --Seed=<integer> ]      # optional; sets srand for reproducibility
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $RandomSet = "";
my $Seed;
GetOptions("RandomSet=s" => \$RandomSet, "Seed=i" => \$Seed);
die "--RandomSet=<name> required\n" unless $RandomSet;

srand($Seed) if defined $Seed;

my $Hg38DictFp    = $cfg{paths}{reference_dict};
my $Hg38TwoBitFp  = $cfg{paths}{reference_2bit};
my $VarDirectory  = "$cfg{paths}{data_root}/$RandomSet";
my $NrTotal       = $cfg{random}{n_variants_total} || 10_000_000;
my $GenomeSize    = $cfg{random}{genome_size_bp}   || 3_000_000_000;

make_path($VarDirectory) unless -d $VarDirectory;

# Read chromosome lengths for 1..22, X, Y from the .dict file
print "Read In Chrom Sizes\n";
my %ChromLengths;
open my $D, '<', $Hg38DictFp or die "Can't open $Hg38DictFp: $!\n";
while (<$D>) {
    chomp;
    my @F     = split(/\t/, $_);
    my $Chrom = $F[1] // next;
    $Chrom    =~ s/SN://;
    next unless ($Chrom =~ /^chr\d+$/ || $Chrom eq "chrX" || $Chrom eq "chrY");
    my $Length = $F[2] // next;
    $Length    =~ s/LN://;
    $ChromLengths{$Chrom} = $Length;
}
close $D;

foreach my $Chrom (sort keys %ChromLengths) {
    my $ChromLen   = $ChromLengths{$Chrom};
    my $NrOfVars   = int( ($NrTotal / $GenomeSize) * $ChromLen );
    my $ChromNoChr = $Chrom; $ChromNoChr =~ s/^chr//;
    print "$Chrom\n";

    my $RawBed = "$VarDirectory/$ChromNoChr.bed";
    open my $R, '>', $RawBed or die "Can't open $RawBed: $!\n";

    for (my $I = 1; $I <= $NrOfVars; $I++) {
        my $End   = int(rand($ChromLen)) + 1;
        my $Start = $End - 1;
        my $Ref   = `twoBitToFa $Hg38TwoBitFp stdout -seq=$Chrom -start=$Start -end=$End | tail -n 1`;
        chomp $Ref;
        $Ref = uc($Ref);
        next if $Ref eq "N";

        my @AltBases  = grep { $_ ne $Ref } ("A", "C", "G", "T");
        my $Alt       = $AltBases[ rand @AltBases ];
        my $Variant   = "${ChromNoChr}_${End}_${End}_${Ref}/${Alt}";

        print $R "$ChromNoChr\t$Start\t$End\t$Variant\n";
    }
    close $R;

    my $Sorted = "$VarDirectory/$ChromNoChr.bedsort.bed";
    system("bedtools sort -i $RawBed > $Sorted");
    system("bgzip -c $Sorted > $Sorted.gz");
    system("tabix -p bed $Sorted.gz");
}

print "Done. Random set '$RandomSet' written to $VarDirectory\n";
