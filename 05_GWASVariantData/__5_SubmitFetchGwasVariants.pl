#!/usr/bin/env perl
#===============================================================================
# __5_SubmitFetchGwasVariants.pl
#
# Split the full GWAS Catalog associations TSV into N sub-files and emit one
# SLURM job per sub-file that runs __5_FetchGwasVariants.pl on it.
#
# Inputs (config <paths>):
#   gwas_catalog_file   path to the full catalog TSV
# Resources:
#   <slurm>             default per-job SLURM resources
#   <gwas>n_sublists    number of sub-files (default 100)
#
# Outputs:
#   <gwas_raw_dir>/<basename>.NN.tsv   sub-files
#   <jobs_dir>/FetchGwasVars.NN.job    SLURM jobs
#
# Submit: for f in <jobs_dir>/FetchGwasVars.*.job; do sbatch $f; done
#
# After all jobs finish, run __5_WriteCleanGwasVariants.pl to merge,
# annotate with gnomAD AF, and write the final per-chromosome BEDs.
#===============================================================================
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $GwasCatalog = $cfg{paths}{gwas_catalog_file};
my $GwasRawDir  = $cfg{paths}{gwas_raw_dir};
my $JobDir      = $cfg{paths}{jobs_dir};
my $CodeDir     = $FindBin::Bin;
my $PerlMod     = $cfg{slurm}{perl_module};
my $NSub        = $cfg{gwas}{n_sublists} || 100;

make_path($GwasRawDir) unless -d $GwasRawDir;
make_path($JobDir)     unless -d $JobDir;

die "GWAS catalog file not found: $GwasCatalog\n" unless -f $GwasCatalog;

my $GwasBasis = $GwasCatalog;
$GwasBasis =~ s/\.tsv$//;

# Split the catalog: drop header, split into N-1 equal chunks, name *.NN.tsv
my $SplitCmd = "tail -n +2 $GwasCatalog | "
             . "split -l \$((( \$(wc -l < $GwasCatalog) - 1 ) / $NSub)) "
             . "-d -a 2 --additional-suffix=.tsv - $GwasBasis.";
system($SplitCmd) == 0 or die "split failed: $SplitCmd\n";

# Emit one SLURM job per sub-file
for (my $I = 0; $I <= $NSub; $I++) {
    my $NN     = sprintf("%02d", $I);
    my $Sub    = "$GwasBasis.$NN.tsv";
    next unless -f $Sub;

    my $JobFp  = "$JobDir/FetchGwasVars.$NN.job";
    my $eFp    = "$JobDir/FetchGwasVars.$NN.e";
    my $oFp    = "$JobDir/FetchGwasVars.$NN.o";

    open my $JOB, '>', $JobFp or die "Can't open $JobFp: $!\n";
    print $JOB <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8400M
#SBATCH --error=$eFp
#SBATCH --output=$oFp
#SBATCH --time=03:00:00

module load $PerlMod

perl $CodeDir/__5_FetchGwasVariants.pl \\
    --GwasFilePath=$Sub \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
    close $JOB;
    chmod 0755, $JobFp;
}

print "Split catalog into $NSub sub-files under " . (dirname_of($GwasBasis)) . "\n";
print "Wrote SLURM jobs to $JobDir\n";
print "Submit: for f in $JobDir/FetchGwasVars.*.job; do sbatch \$f; done\n";

sub dirname_of {
    my $p = shift;
    $p =~ s{/[^/]*$}{};
    return $p;
}
