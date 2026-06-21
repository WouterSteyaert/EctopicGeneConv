#!/usr/bin/env perl
#===============================================================================
# Submitter for __1_FetchSequenceRepeatsFromFasta.pl
#
# Generates one SLURM job per (chromosome x 10 Mb chunk x k-mer length).
# Reads paths and parameter lists from 00_Configuration/config.GRCh38.ini.
#
# Usage:
#   export PROJECT_ROOT=/path/to/project
#   perl __1_SubmitFetchSequenceRepeatsFromFasta.pl \
#       --ConfigFile=../00_Configuration/config.GRCh38.ini
#
# Submit:
#   for f in $PROJECT_ROOT/geneconv_complete/jobs/chr*.rep.job; do sbatch $f; done
#===============================================================================
use strict;
use warnings;
use Cwd qw(getcwd);
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $CODE_DIR        = $FindBin::Bin;
my $JOBDIR          = $cfg{paths}{jobs_dir};
my $REF_DICT        = $cfg{paths}{reference_dict};
my @REPEAT_LENGTHS  = split /,/, $cfg{repeat_lengths}{values};
my @CHROMS          = split /,/, $cfg{genome}{chromosomes};
my $CHUNK           = $cfg{genome}{chunk_size};
my $PERL_MOD        = $cfg{slurm}{perl_module};
my $MEM             = $cfg{slurm}{mem_per_cpu};
my $TIME            = '02:00:00';

make_path($JOBDIR) unless -d $JOBDIR;

# --- read chromosome lengths from .dict ---
my %ChromLengths;
open D, '<', $REF_DICT or die "Can't open $REF_DICT: $!";
while (<D>) {
    chomp;
    my @L = split /\t/;
    next unless @L >= 3 && $L[1] =~ /^SN:(chr\w+)$/;
    my $chr = $1; $chr =~ s/^chr//;
    next unless grep { $_ eq $chr } @CHROMS;
    my ($len) = $L[2] =~ /LN:(\d+)/;
    $ChromLengths{$chr} = $len;
}
close D;

# --- emit one job per (chr x chunk x k) ---
my $n = 0;
foreach my $k (@REPEAT_LENGTHS) {
    foreach my $chr (@CHROMS) {
        next unless exists $ChromLengths{$chr};
        for (my $start = 0; $start <= $ChromLengths{$chr}; $start += $CHUNK) {
            my $end = $start + $CHUNK - 1;
            my $base = "$JOBDIR/chr${chr}_${k}.${start}_${end}.rep";
            open my $jf, '>', "$base.job" or die $!;
            print $jf <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=$MEM
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=$TIME

module load $PERL_MOD

perl $CODE_DIR/__1_FetchSequenceRepeatsFromFasta.pl \\
    --CHROM=$chr --REP_LEN=$k --QUERY_START=$start --QUERY_END=$end \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
            close $jf;
            chmod 0755, "$base.job";
            $n++;
        }
    }
}
print "Wrote $n jobs to $JOBDIR\n";
print "Submit: for f in $JOBDIR/chr*.rep.job; do sbatch \$f; done\n";
