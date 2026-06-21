#!/usr/bin/env perl
#===============================================================================
# 5_PrepareInputUnfazed.pl (per-trio worker)
#
# For one trio, produce the inputs Unfazed needs:
#   - <TrioId>.ped                       PED with sex / trio relationships
#   - <Trio>.input.chr.vcf               denovo VCF, NR/NV -> DP/AD, 'chr' prefix
#   - <RD>.input.chr.vcf.gz              full RD VCF (offspring) filtered & re-keyed
#   - <TrioId>.unfazed.vcf               Unfazed output (parent-of-origin VCF)
#
# Requires bcftools and unfazed in PATH.
#
# Inputs (from config <paths>):
#   denovo_cohort_file        Genomics England denovo LabKey TSV (V9)
#   labkey_rd_interpreted     LabKey rare-disease interpreted TSV
#                             (PlateKey in col 2, BamPath col 18, VcfPath col 19)
#   reference_fasta_full      single-file GRCh38 FASTA (chr-prefixed) for Unfazed
#   denovo_work_dir           top-level DNM work directory
#
# Usage:
#   perl 5_PrepareInputUnfazed.pl --TrioIdSelect=<TrioId>
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use IO::Zlib;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $TrioIdSelect = "";
GetOptions("TrioIdSelect=s" => \$TrioIdSelect);
die "--TrioIdSelect=<id> required\n" unless $TrioIdSelect;

my $TrioLabKeyFileP = $cfg{paths}{denovo_cohort_file};
my $RdLabKeyFileP   = $cfg{paths}{labkey_rd_interpreted};
my $WorkDir         = $cfg{paths}{denovo_work_dir};
my $RefFastaFull    = $cfg{paths}{reference_fasta_full};
my $FilteredDNMsBed = "$WorkDir/ListsOfFilteredDnms/_FilteredDnms.sorted.all.bed";
my $PhasingDir      = "$WorkDir/phasing";
my $OutputDir       = "$PhasingDir/$TrioIdSelect/";
my $PedFilePath     = "${OutputDir}${TrioIdSelect}.ped";

make_path($OutputDir) unless -d $OutputDir;

my %TrioMembers = ();
my %PlateKeys   = ();
my %RdVcfPaths  = ();
my %RdBamPaths  = ();
my %FilteredDNMs = ();
my %AllDNMs     = ();
my $TrioVcfPath = "";
my $ChildSex   = 1;

# Step 1: locate the trio in the denovo LabKey
print "Read In VCF paths\n";
my $LineNr = 0;
open my $C, '<', $TrioLabKeyFileP or die "Can't open $TrioLabKeyFileP: $!\n";
while (<$C>) {
    if ($LineNr) {
        chomp;
        my @F = split(/\t/, $_);
        if ($F[11] eq "GRCh38") {
            my $TrioId  = $F[0];
            my $PlateKey = $F[3];
            my $Member  = $F[5];
            my $Sex     = $F[10];
            my $VcfPath = $F[12];
            my $PedSex  = ($Sex eq "MALE") ? 1 : ($Sex eq "FEMALE" ? 2 : 1);
            if ($TrioId eq $TrioIdSelect) {
                $TrioMembers{$Member}   = $PlateKey;
                $TrioVcfPath           = $VcfPath;
                $ChildSex              = $PedSex;
                $PlateKeys{$PlateKey}  = undef;
                print "$PlateKey\t$Member\n";
            }
        }
    }
    $LineNr++;
}
close $C;

# Step 2: full VCF + BAM paths from RD-interpreted LabKey
print "Read In Complete VCF paths\n";
$LineNr = 0;
open my $C2, '<', $RdLabKeyFileP or die "Can't open $RdLabKeyFileP: $!\n";
while (<$C2>) {
    if ($LineNr) {
        chomp;
        my @F = split(/\t/, $_);
        my ($PlateKey, $BamPath, $VcfPath) = ($F[1], $F[17], $F[18]);
        if (exists $PlateKeys{$PlateKey}) {
            $RdVcfPaths{$PlateKey} = $VcfPath;
            $RdBamPaths{$PlateKey} = $BamPath;
        }
    }
    $LineNr++;
}
close $C2;

# Step 3: filtered DNMs targeted at this proband
print "Read In Filtered DNMs $FilteredDNMsBed\n";
open my $F, '<', $FilteredDNMsBed or die "Can't open $FilteredDNMsBed: $!\n";
while (<$F>) {
    chomp;
    my ($Chrom, undef, $VarPos, $Info) = split(/\t/, $_);
    my ($Change, $Samples) = split(/;/, $Info);
    $Change  =~ s/Change=//;
    $Samples =~ s/Samples=//;
    foreach my $SampleId (split(/~/, $Samples)) {
        if ($SampleId eq $TrioMembers{"Offspring"}) {
            $FilteredDNMs{"${Chrom}_${VarPos}_${Change}"} = undef;
        }
    }
}
close $F;

# Step 4: PED
print "Create PEDs\n";
open my $P, '>', $PedFilePath or die "Can't open $PedFilePath: $!\n";
print $P "#Family-ID\tIndividual-ID\tPaternal-ID\tMaternal-ID\tGender\n";
print $P "$TrioIdSelect\t$TrioMembers{Offspring}\t$TrioMembers{Father}\t$TrioMembers{Mother}\t$ChildSex\n";
print $P "$TrioIdSelect\t$TrioMembers{Father}\t0\t0\t1\n";
print $P "$TrioIdSelect\t$TrioMembers{Mother}\t0\t0\t2\n";
close $P;

chdir($OutputDir);

# Step 5: denovo trio VCF (subset to trio + rewrite NR/NV -> DP/AD)
print "Create Denovo VCF\n";
my $TrioChrVcf = basename($TrioVcfPath);
$TrioChrVcf =~ s/\.gz//;

system("bcftools view -s $TrioMembers{Offspring},$TrioMembers{Father},$TrioMembers{Mother} $TrioVcfPath > $TrioChrVcf");

my $FilteredTrioChrVcf = $TrioChrVcf;
$FilteredTrioChrVcf =~ s/\.vcf/\.input\.chr\.vcf/;

open my $O,  '>', $FilteredTrioChrVcf or die "Can't open $FilteredTrioChrVcf: $!\n";
open my $IH, '<', $TrioChrVcf         or die "Can't open $TrioChrVcf: $!\n";
while (<$IH>) {
    chomp;
    if (/^#/) {
        s/##FORMAT=<ID=NR,/##FORMAT=<ID=DP,/;
        s/##FORMAT=<ID=NV,/##FORMAT=<ID=AD,/;
        print $O "$_\n";
        next;
    }
    my @LV = split(/\t/, $_);
    my ($Chrom, $Pos, undef, $Ref, $Alt) = @LV[0..4];
    my @Base   = @LV[0..7];
    my $Base   = join("\t", @Base);
    my $Format = $LV[8];
    my ($S1, $S2, $S3) = @LV[9..11];
    my @S1 = split(/:/, $S1);
    my @S2 = split(/:/, $S2);
    my @S3 = split(/:/, $S3);

    $AllDNMs{"${Chrom}_${Pos}_${Ref}/${Alt}"} = undef;

    my @FF = split(/:/, $Format);
    my %FF; for (my $i = 0; $i < @FF; $i++) { $FF{$FF[$i]} = $i; }

    next unless exists $FilteredDNMs{"${Chrom}_${Pos}_${Ref}/${Alt}"};

    my $S1_DP = $S1[$FF{NR}];
    my $S1_AD = ($S1[$FF{NR}] - $S1[$FF{NV}]) . "," . $S1[$FF{NV}];
    my $S2_DP = $S2[$FF{NR}];
    my $S2_AD = ($S2[$FF{NR}] - $S2[$FF{NV}]) . "," . $S2[$FF{NV}];
    my $S3_DP = $S3[$FF{NR}];
    my $S3_AD = ($S3[$FF{NR}] - $S3[$FF{NV}]) . "," . $S3[$FF{NV}];

    $Format =~ s/^PS://;
    $Format =~ s/GT:GQ:GOF:NR:GL:NV:DE_NOVO_FLAG/GT:GQ:GOF:DP:GL:AD:DE_NOVO_FLAG/;

    die "comma in NV: $_" if $S1[$FF{NV}] =~ /,/;

    print $O join("\t",
        "chr$Base", $Format,
        "$S1[$FF{GT}]:$S1[$FF{GQ}]:$S1[$FF{GOF}]:$S1_DP:$S1[$FF{GL}]:$S1_AD:$S1[$FF{DE_NOVO_FLAG}]",
        "$S2[$FF{GT}]:$S2[$FF{GQ}]:$S2[$FF{GOF}]:$S2_DP:$S2[$FF{GL}]:$S2_AD:$S2[$FF{DE_NOVO_FLAG}]",
        "$S3[$FF{GT}]:$S3[$FF{GQ}]:$S3[$FF{GOF}]:$S3_DP:$S3[$FF{GL}]:$S3_AD:$S3[$FF{DE_NOVO_FLAG}]",
    ) . "\n";
}
close $IH;
close $O;

# Step 6: full sites VCF (offspring) — for nearby-phased-variant support
print "Create Complete VCF\n";
my $RdVcfName = basename($RdVcfPaths{$TrioMembers{"Offspring"}});
$RdVcfName =~ s/\.gz//;

system("zgrep -v 'JTFH01001724.1' $RdVcfPaths{$TrioMembers{Offspring}} | grep -v 'JTFH01' | grep -v 'KN707' | "
     . "bcftools view -s $TrioMembers{Offspring},$TrioMembers{Father},$TrioMembers{Mother} > $RdVcfName");

my $FilteredRdVcfName = $RdVcfName;
$FilteredRdVcfName    =~ s/\.vcf/\.input\.chr\.vcf/;

open my $OR, '>', $FilteredRdVcfName or die "Can't open $FilteredRdVcfName: $!\n";
open my $IR, '<', $RdVcfName         or die "Can't open $RdVcfName: $!\n";
while (<$IR>) {
    chomp;
    if (/^#/) {
        s/##FORMAT=<ID=NR,/##FORMAT=<ID=DP,/;
        s/##FORMAT=<ID=NV,/##FORMAT=<ID=AD,/;
        print $OR "$_\n";
        next;
    }
    my @LV = split(/\t/, $_);
    my ($Chrom, $Pos, undef, $Ref, $Alt) = @LV[0..4];
    my @Base   = @LV[0..7];
    my $Base   = join("\t", @Base);
    my $Format = $LV[8];
    my ($S1, $S2, $S3) = @LV[9..11];

    next if $Alt =~ /,/;
    next if exists $AllDNMs{"chr${Chrom}_${Pos}_${Ref}/${Alt}"};
    next unless ($Chrom =~ /^chr\d+$/ || $Chrom =~ /^chr[XY]$/);

    my @FF = split(/:/, $Format);
    my %FF; for (my $i = 0; $i < @FF; $i++) { $FF{$FF[$i]} = $i; }

    my @S1 = split(/:/, $S1);
    my @S2 = split(/:/, $S2);
    my @S3 = split(/:/, $S3);

    next if grep { /,/ } ($S1[$FF{NR}], $S2[$FF{NR}], $S3[$FF{NR}],
                          $S1[$FF{NV}], $S2[$FF{NV}], $S3[$FF{NV}]);

    my $mk_dp_ad = sub {
        my @S = @_;
        my $DP = $S[$FF{NR}];
        my $RefDP = ($S[$FF{NR}] =~ /\d+/ && $S[$FF{NV}] =~ /\d+/)
                    ? $S[$FF{NR}] - $S[$FF{NV}] : ".";
        my $AD = "$RefDP,$S[$FF{NV}]";
        return ($DP, $AD);
    };
    my ($S1_DP, $S1_AD) = $mk_dp_ad->(@S1);
    my ($S2_DP, $S2_AD) = $mk_dp_ad->(@S2);
    my ($S3_DP, $S3_AD) = $mk_dp_ad->(@S3);

    $Format =~ s/GT:GL:GOF:GQ:NR:NV/GT:GL:GOF:GQ:DP:AD/;

    print $OR join("\t",
        $Base, $Format,
        "$S1[$FF{GT}]:$S1[$FF{GL}]:$S1[$FF{GOF}]:$S1[$FF{GQ}]:$S1_DP:$S1_AD",
        "$S2[$FF{GT}]:$S2[$FF{GL}]:$S2[$FF{GOF}]:$S2[$FF{GQ}]:$S2_DP:$S2_AD",
        "$S3[$FF{GT}]:$S3[$FF{GL}]:$S3[$FF{GOF}]:$S3[$FF{GQ}]:$S3_DP:$S3_AD",
    ) . "\n";
}
close $IR;
close $OR;

system("bgzip -c $FilteredRdVcfName > $FilteredRdVcfName.gz");
system("tabix -p vcf $FilteredRdVcfName.gz");

# Step 7: run Unfazed
print "Unfazed\n";
my $UnfazedCmd = "unfazed -d $FilteredTrioChrVcf -s $FilteredRdVcfName.gz -p $PedFilePath "
               . "--bam-pairs "
               . "$TrioMembers{Offspring}:$RdBamPaths{$TrioMembers{Offspring}} "
               . "$TrioMembers{Father}:$RdBamPaths{$TrioMembers{Father}} "
               . "$TrioMembers{Mother}:$RdBamPaths{$TrioMembers{Mother}} "
               . "-t 1 -o vcf --outfile $TrioIdSelect.unfazed.vcf -g 38 -r $RefFastaFull";

print "$UnfazedCmd\n";
system($UnfazedCmd);
