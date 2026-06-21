#!/usr/bin/env perl
#===============================================================================
# Submitter for __7_FindGeneConvNotInShorter.pl
# One SLURM job per (chromosome x 10 Mb chunk x diff_pair).
#
# diff_pairs come from config (e.g. "17_19,19_21,21_31,...") and define the
# short_long pairs across which set-differences are computed.
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();
my $CODE_DIR    = $FindBin::Bin;
my $JOBDIR      = $cfg{paths}{jobs_dir};
my $REF_DICT    = $cfg{paths}{reference_dict};
my @DIFF_PAIRS  = split /,/, $cfg{repeat_lengths}{diff_pairs};
my @CHROMS      = split /,/, $cfg{genome}{chromosomes};
my $CHUNK       = $cfg{genome}{chunk_size};
my $PERL_MOD    = $cfg{slurm}{perl_module};
my $MEM         = $cfg{slurm}{mem_per_cpu};

make_path($JOBDIR) unless -d $JOBDIR;

# Read chromosome lengths
my %CL;
open D, '<', $REF_DICT or die "Can't open $REF_DICT: $!";
while (<D>) {
    chomp;
    my @L = split /\t/;
    next unless @L >= 3 && $L[1] =~ /^SN:(chr\w+)$/;
    my $c = $1; $c =~ s/^chr//;
    next unless grep { $_ eq $c } @CHROMS;
    my ($len) = $L[2] =~ /LN:(\d+)/;
    $CL{$c} = $len;
}
close D;

my $n = 0;
foreach my $pair (@DIFF_PAIRS) {
    foreach my $chr (@CHROMS) {
        next unless exists $CL{$chr};
        for (my $start = 1; $start <= $CL{$chr}; $start += $CHUNK) {
            my $end = $start + $CHUNK;
            my $base = "$JOBDIR/repdiff.chr${chr}.${start}_${end}.${pair}";
            open my $jf, '>', "$base.job" or die $!;
            print $jf <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=$MEM
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=06:00:00

module load $PERL_MOD

perl $CODE_DIR/__7_FindGeneConvNotInShorter.pl \\
    --CHROM=$chr --REP_LENGTHS=$pair --START=$start --END=$end \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
            close $jf;
            chmod 0755, "$base.job";
            $n++;
        }
    }
}
print "Wrote $n jobs to $JOBDIR\n";
