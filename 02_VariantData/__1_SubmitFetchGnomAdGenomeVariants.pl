#!/usr/bin/env perl
#===============================================================================
# Submitter for __1_FetchGnomAdGenomeVariants.pl
# One SLURM job per chromosome (1..22, X, Y).
#
# Reads paths from 00_Configuration/config.GRCh38.ini.
#
# Usage:
#   export PROJECT_ROOT=/path/to/project
#   perl __1_SubmitFetchGnomAdGenomeVariants.pl \
#       --ConfigFile=../00_Configuration/config.GRCh38.ini
#   for f in $PROJECT_ROOT/geneconv_complete/jobs/FetchGnomadGenome.*.job; do sbatch $f; done
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $CODE_DIR = $FindBin::Bin;
my $JOBDIR   = $cfg{paths}{jobs_dir};
my @CHROMS   = split /,/, $cfg{genome}{chromosomes};
my $PERL_MOD = $cfg{slurm}{perl_module};
my $HTS_MOD  = $cfg{modules}{htslib};

make_path($JOBDIR) unless -d $JOBDIR;

my $n = 0;
foreach my $chr (@CHROMS) {
    my $base = "$JOBDIR/FetchGnomadGenome.$chr";
    open my $jf, '>', "$base.job" or die "Can't open $base.job: $!";
    print $jf <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48400M
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=09:00:00

module load $PERL_MOD
module load $HTS_MOD

perl $CODE_DIR/__1_FetchGnomAdGenomeVariants.pl \\
    --CHROM=$chr \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
    close $jf;
    chmod 0755, "$base.job";
    $n++;
}

print "Wrote $n jobs to $JOBDIR\n";
print "Submit: for f in $JOBDIR/FetchGnomadGenome.*.job; do sbatch \$f; done\n";
