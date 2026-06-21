#!/usr/bin/env perl
#===============================================================================
# 5_SubmitPrepareInputUnfazed.pl
#
# Emit one LSF job per trio that runs 5_PrepareInputUnfazed.pl (i.e. runs the
# Unfazed parent-of-origin phasing for that trio).  Reads the per-proband DNM
# count file produced by 3_ProcessDNMs.pl and includes only probands inside the
# 30-150 DNM window.
#
# Pass --Mode=fresh (default) for first submission, --Mode=retry to only
# re-emit jobs whose unfazed.vcf is missing or empty (uses .re.job suffix), or
# --Mode=retry2 for the second retry (.re.re.job suffix).
#
# Outputs:
#   denovo_work_dir/_Jobs/<TrioId>.phasing.job    (--Mode=fresh)
#   denovo_work_dir/_Jobs/<TrioId>.phasing.re.job (--Mode=retry)
#   denovo_work_dir/_Jobs/<TrioId>.phasing.re.re.job (--Mode=retry2)
#
# Submit: for f in <denovo_work_dir>/_Jobs/*.phasing*.job; do bsub < $f; done
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $Mode = "fresh";
GetOptions("Mode=s" => \$Mode);
die "--Mode must be fresh|retry|retry2\n" unless $Mode =~ /^(fresh|retry|retry2)$/;

my $WorkDir       = $cfg{paths}{denovo_work_dir};
my $NrOfDNMsFp    = "$WorkDir/ListsOfFilteredDnms/_DnmNumbersPerSample.txt";
my $JobDirectory  = "$WorkDir/_Jobs";
my $PhasingDir    = "$WorkDir/phasing";
my $CodeDir       = $FindBin::Bin;

my $LsfQueue      = $cfg{lsf}{queue};
my $LsfProject    = $cfg{lsf}{project};
my $LsfMem        = ($Mode eq "fresh") ? $cfg{lsf}{mem_phasing} : $cfg{lsf}{mem_phasing_retry};
my $PerlModule    = $cfg{modules}{perl_lsf};
my $BcfModule     = $cfg{modules}{bcftools_lsf};
my $MinNrOfDnms   = $cfg{denovo}{min_dnms} || 30;
my $MaxNrOfDnms   = $cfg{denovo}{max_dnms} || 150;

my %Suffix = (fresh => "phasing", retry => "phasing.re", retry2 => "phasing.re.re");
my $Suf    = $Suffix{$Mode};

make_path($JobDirectory) unless -d $JobDirectory;

my $LineNr = 0;
my $NrJobs = 0;
open my $N, '<', $NrOfDNMsFp or die "Can't open $NrOfDNMsFp: $!\n";
while (<$N>) {
    chomp;
    if ($LineNr) {
        my ($Proband, $TrioId, $NrOfDNMs) = split(/\t/, $_);
        next unless ($NrOfDNMs >= $MinNrOfDnms && $NrOfDNMs <= $MaxNrOfDnms);

        if ($Mode ne "fresh") {
            my $UnfazedVcf = "$PhasingDir/$TrioId/$TrioId.unfazed.vcf";
            next if (-f $UnfazedVcf && -s $UnfazedVcf);
        } else {
            make_path("$PhasingDir/$TrioId") unless -d "$PhasingDir/$TrioId";
        }

        my $JobFp = "$JobDirectory/$TrioId.$Suf.job";
        my $eFp   = "$JobDirectory/$TrioId.$Suf.e";
        my $oFp   = "$JobDirectory/$TrioId.$Suf.o";

        open my $JOB, '>', $JobFp or die "Can't open $JobFp: $!\n";
        print $JOB <<"JOB";
#!/bin/bash
#BSUB -q $LsfQueue
#BSUB -P $LsfProject
#BSUB -o $oFp
#BSUB -e $eFp
#BSUB -J $TrioId.$Suf
#BSUB -R "rusage[mem=$LsfMem] span[hosts=1]"
#BSUB -n 1

module load $PerlModule
module load $BcfModule

perl $CodeDir/5_PrepareInputUnfazed.pl \\
    --TrioIdSelect=$TrioId \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
        close $JOB;
        chmod 0755, $JobFp;
        $NrJobs++;
    }
    $LineNr++;
}
close $N;

print "Wrote $NrJobs LSF jobs to $JobDirectory (mode=$Mode, suffix=$Suf)\n";
print "Submit: for f in $JobDirectory/*.$Suf.job; do bsub < \$f; done\n";
