#!/usr/bin/env perl
#===============================================================================
# Submitter for __9_MakeConvFilesFromDiffs.pl
# One SLURM job per chromosome — builds the final per-k cumulative BED files
# from the gcdiffs set-difference outputs.
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();
my $CODE_DIR  = $FindBin::Bin;
my $JOBDIR    = $cfg{paths}{jobs_dir};
my @CHROMS    = split /,/, $cfg{genome}{chromosomes};
my $PERL_MOD  = $cfg{slurm}{perl_module};
my $HTS_MOD   = $cfg{modules}{htslib};
my $MEM       = $cfg{slurm}{mem_per_cpu};

make_path($JOBDIR) unless -d $JOBDIR;

my $n = 0;
foreach my $chr (@CHROMS) {
    my $base = "$JOBDIR/makeconv.chr${chr}";
    open my $jf, '>', "$base.job" or die $!;
    print $jf <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=$MEM
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=04:00:00

module load $PERL_MOD
module load $HTS_MOD

perl $CODE_DIR/__9_MakeConvFilesFromDiffs.pl \\
    --CHROM=$chr \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
    close $jf;
    chmod 0755, "$base.job";
    $n++;
}
print "Wrote $n jobs to $JOBDIR\n";
