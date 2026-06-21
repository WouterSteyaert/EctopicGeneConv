#!/usr/bin/env perl
#===============================================================================
# Submitter for the concat-and-tabix step on the gcdiffs output.
#
# For each (short_long k-pair, chromosome) it writes a SLURM job that:
#   1. concatenates all per-chunk .changes.gcdiffs.txt.gz files
#   2. concatenates all per-chunk .positions.gcdiffs.txt.gz files
#   3. sorts each, bgzips, and tabix-indexes
#
# Uses 00_Configuration/config.GRCh38.ini for paths and parameter lists.
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my @DIFF_PAIRS = split /,/, $cfg{repeat_lengths}{diff_pairs};
my @CHROMS     = split /,/, $cfg{genome}{chromosomes};
my $JOBDIR     = $cfg{paths}{jobs_dir};
my $GCDIFFDIR  = $cfg{paths}{gcdiffs_dir};
my $TEMPDIR    = "$cfg{paths}{data_root}/temp";
my $PERL_MOD   = $cfg{slurm}{perl_module};
my $HTS_MOD    = $cfg{modules}{htslib};
my $MEM        = $cfg{slurm}{mem_per_cpu};

make_path($JOBDIR)  unless -d $JOBDIR;
make_path($TEMPDIR) unless -d $TEMPDIR;

foreach my $pair (@DIFF_PAIRS) {
    my $SubDir = "$GCDIFFDIR/$pair/";
    next unless -d $SubDir;

    foreach my $chr (@CHROMS) {
        my $ConcatChanges      = "$SubDir/$chr.changes.gcdiffs.sorted.bed";
        my $ConcatChangesTmp   = "$SubDir/$chr.changes.gcdiffs.bed";
        my $ConcatPositions    = "$SubDir/$chr.positions.gcdiffs.sorted.bed";
        my $ConcatPositionsTmp = "$SubDir/$chr.positions.gcdiffs.bed";

        my $base = "$JOBDIR/$chr.$pair.gcdiffs";
        open my $jf, '>', "$base.job" or die "Can't open $base.job: $!";

        print $jf <<"HEADER";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=$MEM
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=09:00:00

module load $PERL_MOD
module load $HTS_MOD

# Clean any partial outputs from previous runs
[ -f $ConcatChangesTmp ]   && rm -f $ConcatChangesTmp
[ -f $ConcatPositionsTmp ] && rm -f $ConcatPositionsTmp
HEADER

        # Per-chunk file list assembled at submit-time (chromosome-specific)
        opendir(my $dh, $SubDir) or die "Can't open $SubDir: $!";
        my $chr_q = quotemeta($chr);
        while (my $file = readdir($dh)) {
            if ($file =~ /^${chr_q}\..*\.changes\.gcdiffs\.txt\.gz$/) {
                print $jf "zcat $SubDir/$file >> $ConcatChangesTmp\n";
            }
            if ($file =~ /^${chr_q}\..*\.positions\.gcdiffs\.txt\.gz$/) {
                print $jf "zcat $SubDir/$file >> $ConcatPositionsTmp\n";
            }
        }
        closedir($dh);

        print $jf <<"FOOTER";

sort -k1,1 -k2,2n -T $TEMPDIR $ConcatChangesTmp > $ConcatChanges
bgzip -f -c $ConcatChanges > $ConcatChanges.gz
tabix -p bed $ConcatChanges.gz

sort -k1,1 -k2,2n -T $TEMPDIR $ConcatPositionsTmp > $ConcatPositions
bgzip -f -c $ConcatPositions > $ConcatPositions.gz
tabix -p bed $ConcatPositions.gz

rm -f $ConcatChangesTmp $ConcatPositionsTmp $ConcatChanges $ConcatPositions
FOOTER
        close $jf;
        chmod 0755, "$base.job";
    }
}

print "Wrote SLURM jobs to $JOBDIR\n";
