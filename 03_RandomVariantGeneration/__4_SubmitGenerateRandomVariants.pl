#!/usr/bin/env perl
#===============================================================================
# __4_SubmitGenerateRandomVariants.pl
#
# Emit one SLURM job per random set that runs __4_GenerateRandomVariants on
# that set.  Random set names are read from <random>set_names in the config.
#
# Reproducibility note: by default each job uses Perl's default RNG state and
# the resulting random set is not bit-reproducible across reruns.  Pass
# --SetSeeds=<comma-list> matching the order of <random>set_names to make
# each job deterministic; for example
#     --SetSeeds=1,2,3
# gives random1 srand(1), random2 srand(2), random3 srand(3).
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $SetSeeds = "";
GetOptions("SetSeeds=s" => \$SetSeeds);

my $JobDir   = $cfg{paths}{jobs_dir};
my $CodeDir  = $FindBin::Bin;
my $PerlMod  = $cfg{slurm}{perl_module};
my $HtsMod   = $cfg{modules}{htslib};
my $BedMod   = $cfg{modules}{bedtools};
my @Sets     = split /,/, ($cfg{random}{set_names} // "random1,random2,random3");
my @Seeds    = $SetSeeds ? split(/,/, $SetSeeds) : ();

make_path($JobDir) unless -d $JobDir;

for (my $i = 0; $i < @Sets; $i++) {
    my $Set    = $Sets[$i];
    my $Seed   = $Seeds[$i];   # undef -> no --Seed= passed
    my $base   = "$JobDir/FetchRandom.$Set";
    my $SeedArg = defined $Seed ? " --Seed=$Seed" : "";

    open my $JOB, '>', "$base.job" or die "Can't open $base.job: $!\n";
    print $JOB <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16400M
#SBATCH --error=$base.e
#SBATCH --output=$base.o
#SBATCH --time=72:00:00

module load $PerlMod
module load $HtsMod
module load $BedMod

perl $CodeDir/__4_GenerateRandomVariants.pl \\
    --RandomSet=$Set$SeedArg \\
    --ConfigFile=$cfg{_meta}{config_file} --ProjectRoot=$cfg{_meta}{project_root}
JOB
    close $JOB;
    chmod 0755, "$base.job";
}

print "Wrote " . scalar(@Sets) . " SLURM jobs to $JobDir\n";
print "Submit: for f in $JobDir/FetchRandom.*.job; do sbatch \$f; done\n";
