#!/usr/bin/env perl
#===============================================================================
# CorrectForMht.pl
#
# Apply Bonferroni multiple-testing correction to the per-(k, AF-bin)
# p-values produced by __10_SummarizeContingencyStatistics (step 06).
#
# For each <StatsRoot>/stats_<context>/ directory, the script emits a single
# R driver that:
#   - reads every *.R.out.txt file (the per-VarSet summary tables);
#   - sets m = nrow(Table) (i.e. number of (k, AF-bin) rows actually present:
#       m = 90 for stratified gnomAD / DNM / GWAS (10 k × 9 AF bins),
#       m = 10 for unstratified random control (10 k, single bin));
#   - applies p.adjust(method="bonferroni") to both the positional chi-squared
#     and the concordance chi-squared p-values;
#   - shifts the log-scale p-value exponents by log10(m);
#   - writes the corrected table to *.corr.R.out.txt.
#
# Each variant set (gnomAD, GWAS, DNM ALL/FATHER/MOTHER, random) is corrected
# independently because each addresses an independent biological question.
#
# Usage:
#   # Genome / GWAS (use the default stats_dir from config):
#   perl CorrectForMht.pl --ConfigFile=../00_Configuration/config.GRCh38.ini
#
#   # DNM run (separate stats root):
#   perl CorrectForMht.pl --StatsRoot=<dnm_export_root> \
#       --ConfigFile=../00_Configuration/config.GRCh38.ini
#
# Then run the emitted R driver(s):
#   for f in <root>/stats_*/CorrectForMht.R; do Rscript $f; done
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $StatsRoot = "";
GetOptions("StatsRoot=s" => \$StatsRoot);
$StatsRoot ||= $cfg{paths}{stats_dir};

my @Contexts = split /,/, ($cfg{enrichment}{contexts} // "allmapp,segdupmapp,nosegdupmapp");

my $NEmitted = 0;
foreach my $Context (@Contexts) {
    my $StatDir   = "$StatsRoot/stats_$Context/";
    next unless -d $StatDir;

    my $RFilePath = "${StatDir}CorrectForMht.R";
    open my $R, '>', $RFilePath or die "Can't open $RFilePath: $!\n";

    opendir(my $SD, $StatDir) or die "Can't open $StatDir: $!\n";
    my $count = 0;
    while (my $File = readdir($SD)) {
        next unless $File =~ /^(.*)\.R\.out\.txt$/;
        next if     $File =~ /\.corr\.R\.out\.txt$/;
        my $FileBase         = $1;
        my $StatFilePath     = "$StatDir$File";
        my $StatCorrFilePath = "$StatDir$FileBase.corr.R.out.txt";

        print $R <<"R_BLOCK";
Table <- read.table("$StatFilePath", header = TRUE, sep = "\\t")
m <- nrow(Table)
Table\$log_p_value_exponent_bonferroni          <- Table\$log_p_value_exponent          + log10(m)
Table\$chisq_test_p_bonferroni                  <- p.adjust(Table\$chisq_test_p, method = "bonferroni")
Table\$log_varmatch_p_value_exponent_bonferroni <- Table\$log_p_varmatch_value_exponent + log10(m)
Table\$chisq_test_varmatch_p_bonferroni         <- p.adjust(Table\$chisq_test_varmatch_p, method = "bonferroni")
write.table(Table, file = "$StatCorrFilePath", quote = FALSE, row.names = FALSE, sep = "\\t")
R_BLOCK
        $count++;
    }
    closedir($SD);
    close $R;

    if ($count) {
        print "Wrote $RFilePath ($count tables to correct)\n";
        $NEmitted++;
    } else {
        unlink $RFilePath;   # nothing to correct in this context
    }
}

if ($NEmitted) {
    print "\nRun each R driver to apply the correction:\n";
    print "  for f in $StatsRoot/stats_*/CorrectForMht.R; do Rscript \$f; done\n";
} else {
    print "No *.R.out.txt summary tables found under $StatsRoot — run step 06 first.\n";
}
