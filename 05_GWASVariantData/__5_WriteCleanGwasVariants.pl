#!/usr/bin/env perl
#===============================================================================
# __5_WriteCleanGwasVariants.pl
#
# Take all per-shard `*.processed.tsv` files produced by __5_FetchGwasVariants,
# merge them, attach the gnomAD AF for each variant (from the per-chrom BED
# written by step 02), and emit clean per-chromosome BEDs stratified by
# phenotype class+type:
#
#   <gwas_clean_dir>/All/<chr>.bed
#   <gwas_clean_dir>/Disease_Risk/<chr>.bed
#   <gwas_clean_dir>/Disease_Protection/<chr>.bed
#   <gwas_clean_dir>/Trait_All/<chr>.bed
#
# Each is bedtools-sorted, bgzipped and tabix-indexed.
#
# A variant whose phenotype is not present in classified_phenotypes is dropped.
# A variant with no gnomAD match keeps a frequency of 0 in the BED; downstream
# enrichment scripts (06+) require a gnomAD match for inclusion.
#
# Inputs (config <paths>):
#   gwas_raw_dir              directory containing *.processed.tsv shards and
#                             the classified_phenotypes CSV
#   gwas_clean_dir            output directory
#   classified_phenotypes     CSV mapping phenotype label -> "Trait" or "Disease"
#   gnomad_genome_dir         per-chrom <chr>.bed from step 02 (pre-mappability)
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $GwasRawDir       = $cfg{paths}{gwas_raw_dir};
my $GwasCleanDir     = $cfg{paths}{gwas_clean_dir};
my $ClassifiedPheno  = $cfg{paths}{classified_phenotypes};
my $GnomAdDir        = $cfg{paths}{gnomad_genome_dir};

make_path($GwasCleanDir) unless -d $GwasCleanDir;

my %ClassifiedPhenos = ();
my %GwasVariants     = ();

# Phenotype -> Class (Trait | Disease)
print "Read In Classified Phenotypes\n";
my $LineNr = 0;
open my $C, '<', $ClassifiedPheno or die "Can't open $ClassifiedPheno: $!\n";
while (<$C>) {
    if ($LineNr) {
        s/[\n\r"]//g;
        if (/^(.*),(Trait|Disease)$/) {
            my $Phenotype = $1; $Phenotype =~ s/ $//;
            my $Class     = $2;
            $ClassifiedPhenos{$Phenotype} = $Class;
        }
    }
    $LineNr++;
}
close $C;

# Walk all .processed.tsv shards
print "Iterate Over Directory\n";
opendir(my $RAW, $GwasRawDir) or die "Can't open $GwasRawDir: $!\n";
while (my $File = readdir($RAW)) {
    next unless $File =~ /\.processed\.tsv$/;
    my $fp = "$GwasRawDir/$File";
    open my $F, '<', $fp or die "Can't open $fp: $!\n";
    while (<$F>) {
        s/[\n\r"]//g;
        next unless $_;
        my @F = split(/\t/, $_);
        my $Variant   = $F[3];
        my ($Chrom)   = split(/_/, $Variant);
        my $Phenotype = $F[7]; $Phenotype =~ s/ $//;
        my $Type      = $F[8];
        my $Class     = $ClassifiedPhenos{$Phenotype};

        next unless $Phenotype;
        next unless defined $Class;

        $Type = "All" if $Class eq "Trait";
        $GwasVariants{$Chrom}{"${Class}_${Type}"}{$Variant} = $_;
    }
    close $F;
}
closedir($RAW);

# Per-chromosome write-out, with gnomAD AF lookup
print "Write Out\n";
foreach my $Chrom (sort keys %GwasVariants) {
    print "\t$Chrom\n";

    my %GnomAdVars = ();
    my $GnomAdBed  = "$GnomAdDir/$Chrom.bed";
    open my $G, '<', $GnomAdBed or die "Can't open $GnomAdBed: $!\n";
    while (<$G>) {
        chomp;
        my @F = split(/\t/, $_);
        $GnomAdVars{$F[3]} = $F[4];
    }
    close $G;

    # Combined "All" file
    my $OutDirAll  = "$GwasCleanDir/All";
    my $AllOut     = "$OutDirAll/$Chrom.bed";
    my $SortedAll  = "$OutDirAll/$Chrom.sorted.bed";
    make_path($OutDirAll) unless -d $OutDirAll;
    open my $A, '>', $AllOut or die "Can't open $AllOut: $!\n";

    foreach my $ClassType (sort keys %{ $GwasVariants{$Chrom} }) {
        my $OutDirCT = "$GwasCleanDir/$ClassType";
        my $SubOut   = "$OutDirCT/$Chrom.bed";
        my $SortedSub = "$OutDirCT/$Chrom.sorted.bed";
        make_path($OutDirCT) unless -d $OutDirCT;
        open my $S, '>', $SubOut or die "Can't open $SubOut: $!\n";

        foreach my $Variant (sort keys %{ $GwasVariants{$Chrom}{$ClassType} }) {
            my $GnomAdFreq = $GnomAdVars{$Variant} // 0;
            print $A $GwasVariants{$Chrom}{$ClassType}{$Variant} . "\t$GnomAdFreq\n";
            print $S $GwasVariants{$Chrom}{$ClassType}{$Variant} . "\t$GnomAdFreq\n";
        }
        close $S;
        system("bedtools sort -i $SubOut > $SortedSub");
        system("bgzip -c $SortedSub > $SortedSub.gz");
        system("tabix -p bed $SortedSub.gz");
    }
    close $A;

    system("bedtools sort -i $AllOut > $SortedAll");
    system("bgzip -c $SortedAll > $SortedAll.gz");
    system("tabix -p bed $SortedAll.gz");
}

print "Done. Clean per-chromosome BEDs under $GwasCleanDir\n";
