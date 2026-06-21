#!/usr/bin/env perl
#===============================================================================
# 1_CreateListsOfDenovoVcfs.pl
#
# Build per-batch lists of denovo VCF file paths (one TrioId per line) and a
# matching LSF job file per sublist that runs 2_FilterDNMs.pl on the sublist.
#
# Families missing from the LabKey rare-disease cohort (FamiliesToExclude.txt)
# are skipped, as are families with anything other than exactly one Father and
# one Mother.
#
# Inputs  (paths read from config.GRCh38.ini):
#   denovo_cohort_file                  Genomics England denovo cohort TSV
#   denovo_work_dir/FamiliesToExclude.txt   produced by 0_FindTriosToExclude.pl
#
# Outputs:
#   denovo_work_dir/ListsOfDenovoVcfs/DenovoVcfs_<N>.txt   sublist files
#   denovo_work_dir/_Jobs/FilterDnms.<N>.job               LSF job per sublist
#   denovo_work_dir/ListsOfFilteredDnms/                   created (empty)
#
# Submit:
#   for f in <denovo_work_dir>/_Jobs/FilterDnms.*.job; do bsub < $f; done
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
my $FamiliesToExclude  = "$WorkDir/FamiliesToExclude.txt";
my $DenovoVcfDir       = "$WorkDir/ListsOfDenovoVcfs";
my $JobDirectory       = "$WorkDir/_Jobs";
my $FilteredDnmsDir    = "$WorkDir/ListsOfFilteredDnms";
my $CodeDir            = $FindBin::Bin;

my $LsfQueue           = $cfg{lsf}{queue};
my $LsfProject         = $cfg{lsf}{project};
my $LsfMem             = $cfg{lsf}{mem_filter};
my $PerlModule         = $cfg{modules}{perl_lsf};
my $TriosPerSublist    = $cfg{denovo}{trios_per_sublist} || 100;

my %FamiliesToExcludeH = ();
my %DE_NOVO_VCF_FILEPATHS = ();

for my $dir ($DenovoVcfDir, $JobDirectory, $FilteredDnmsDir) {
    make_path($dir) unless -d $dir;
}

# Clear any prior sublists / jobs so re-running gives a clean batch set
system("rm -f $DenovoVcfDir/DenovoVcfs_*.txt");
system("rm -f $JobDirectory/FilterDnms.*.job $JobDirectory/FilterDnms.*.e $JobDirectory/FilterDnms.*.o");

# Read families to exclude
open my $E, '<', $FamiliesToExclude or die "Can't open $FamiliesToExclude: $!\n";
while (<$E>) {
    next if /^#/;
    chomp;
    $FamiliesToExcludeH{$_} = undef;
}
close $E;

# Read denovo cohort: TrioId, FamilyId, PlateKey, Member, Assembly, VcfPath
open my $L, '<', $LabKeyFileP or die "Can't open $LabKeyFileP: $!\n";
while (<$L>) {
    chomp;
    my @F                  = split(/\t/, $_);
    my $TrioId             = $F[0];
    my $FamilyId           = $F[1];
    my $PlateKey           = $F[3];
    my $Member             = $F[5];
    my $Assembly           = $F[11];
    my $DenovoVcfFilePath  = $F[12];

    if ($Assembly eq "GRCh38" && !exists $FamiliesToExcludeH{$FamilyId}) {
        $DE_NOVO_VCF_FILEPATHS{$TrioId}{$DenovoVcfFilePath}{"$TrioId:$Member:$PlateKey"} = undef;
    }
}
close $L;

# Drop trios that don't have exactly one Father and one Mother
foreach my $TrioId (sort keys %DE_NOVO_VCF_FILEPATHS) {
    foreach my $DenovoVcfFilePath (sort keys %{$DE_NOVO_VCF_FILEPATHS{$TrioId}}) {
        my $NrOfFathers = 0;
        my $NrOfMothers = 0;
        foreach my $FamMember (sort keys %{$DE_NOVO_VCF_FILEPATHS{$TrioId}{$DenovoVcfFilePath}}) {
            $NrOfFathers++ if $FamMember =~ /Father/i;
            $NrOfMothers++ if $FamMember =~ /Mother/i;
        }
        if ($NrOfFathers != 1 || $NrOfMothers != 1) {
            delete $DE_NOVO_VCF_FILEPATHS{$TrioId};
        }
    }
}

# Write VCF-path sublists (TriosPerSublist trios per file)
my $TrioNr     = 0;
my $TrioListNr = 0;
foreach my $TrioId (sort keys %DE_NOVO_VCF_FILEPATHS) {
    foreach my $DenovoVcfFilePath (sort keys %{$DE_NOVO_VCF_FILEPATHS{$TrioId}}) {
        $TrioNr++;
        $TrioListNr      = int(($TrioNr - 1) / $TriosPerSublist);
        my $FamilyString = join("~", sort keys %{$DE_NOVO_VCF_FILEPATHS{$TrioId}{$DenovoVcfFilePath}});

        open my $OUT, '>>', "$DenovoVcfDir/DenovoVcfs_$TrioListNr.txt"
            or die "Can't open sublist $TrioListNr: $!\n";
        print $OUT $DenovoVcfFilePath . "\t" . $FamilyString . "\n";
        close $OUT;
    }
}

# Emit one LSF job per sublist
for (my $I = 0; $I <= $TrioListNr; $I++) {
    my $JobFilePath = "$JobDirectory/FilterDnms.$I.job";
    my $eFilePath   = "$JobDirectory/FilterDnms.$I.e";
    my $oFilePath   = "$JobDirectory/FilterDnms.$I.o";

    open my $JOB, '>', $JobFilePath or die "Can't open $JobFilePath: $!\n";
    print $JOB <<"JOB";
#!/bin/bash
#BSUB -q $LsfQueue
#BSUB -P $LsfProject
#BSUB -o $oFilePath
#BSUB -e $eFilePath
#BSUB -J FilterDnms.$I
#BSUB -R "rusage[mem=$LsfMem] span[hosts=1]"
#BSUB -n 1

module load $PerlModule

perl $CodeDir/2_FilterDNMs.pl \\
    --SubListNr=$I \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
    close $JOB;
    chmod 0755, $JobFilePath;
}

print "Wrote " . ($TrioListNr + 1) . " sublists to $DenovoVcfDir\n";
print "Wrote " . ($TrioListNr + 1) . " LSF jobs to $JobDirectory\n";
print "Submit: for f in $JobDirectory/FilterDnms.*.job; do bsub < \$f; done\n";
