#!/usr/bin/env perl
#===============================================================================
# Submitter for __3_GetPossibleConvPos.pl
# One SLURM job per (k x 5-nt seed) — produces *.repsum.geneconv.coords.txt.gz
# (per-seed file with paralog instances and complement instances).
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
my @KS        = split /,/, $cfg{repeat_lengths}{values};
my $PERL_MOD  = $cfg{slurm}{perl_module};
my $MEM       = $cfg{slurm}{mem_per_cpu};

make_path($JOBDIR) unless -d $JOBDIR;

my @ALPHABET = ('A','C','G','T');
my @SEEDS;
for my $a (@ALPHABET) { for my $b (@ALPHABET) { for my $c (@ALPHABET) {
    for my $d (@ALPHABET) { for my $e (@ALPHABET) {
        push @SEEDS, "$a$b$c$d$e";
    }}
}}}

my $n = 0;
foreach my $k (@KS) {
    foreach my $seed (@SEEDS) {
        my $base = "$JOBDIR/getconvpos.${k}.${seed}";
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

perl $CODE_DIR/__3_GetPossibleConvPos.pl \\
    --REP_LEN=$k --SEED=$seed \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
        close $jf;
        chmod 0755, "$base.job";
        $n++;
    }
}
print "Wrote $n jobs to $JOBDIR\n";
