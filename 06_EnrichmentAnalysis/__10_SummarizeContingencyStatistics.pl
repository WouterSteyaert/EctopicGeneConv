#!/usr/bin/env perl
#===============================================================================
# __10_SummarizeContingencyStatistics.pl
#
# Aggregate the per-window outputs of __1_ComputeContingencyStatistics across
# all chromosomes for each (k, AF-bin) combination, in each of the three
# mappability contexts.  Writes:
#   _All_<W>_<S>_<varlist>_StatsSummary.pos.txt
#       three blocks per context:
#         (1) PercDiff_Pos  = (var_fraction_conv - var_fraction_nonconv)
#                             / var_fraction_nonconv * 100   (k vs AF)
#         (2) ExcessPercPos = (observed - expected) / total_variants * 100
#                             where expected uses the non-conv variant density
#         (3) ExcessPercVar = (concordant - expected_from_discordant) / total_variants * 100
#   _All_<W>_<S>_<varlist>_StatsSummary.R
#       an R driver that performs the chi-squared, binomial, and concordance
#       chi-squared tests and writes _All_..._StatsSummary.R.out.txt with
#       per-(k, AF-bin) p-values + statistical power.
#   _All_<W>_<S>_<varlist>_StatsSummary.job
#       a SLURM job that runs the R script.
#
# Inputs (config):
#   <paths>stats_dir            base directory containing stats_<context>/<k>/...
#   <enrichment>window_size     window size of the worker outputs (default 1Mb)
#   <enrichment>step_size       step size of the worker outputs
#   <enrichment>contexts        comma-list (default allmapp,segdupmapp,nosegdupmapp)
#
# Usage:
#   perl __10_SummarizeContingencyStatistics.pl \
#       --VarListPrint=gnomad~genome \
#       --ConfigFile=../00_Configuration/config.GRCh38.ini
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

my $VarListPrint = "gnomad~genome";
my $AfMode       = "stratified";    # 'stratified' (9 AF bins) or 'unstratified' (single 0_all bin)
GetOptions(
    "VarListPrint=s" => \$VarListPrint,
    "AfMode=s"       => \$AfMode,
);
die "--AfMode must be stratified|unstratified\n" unless $AfMode =~ /^(stratified|unstratified)$/;

my $StatsDir   = $cfg{paths}{stats_dir};
my $WindowSize = $cfg{enrichment}{window_size} || 1_000_000;
my $StepSize   = $cfg{enrichment}{step_size}   || 1_000_000;
my @Contexts   = split /,/, ($cfg{enrichment}{contexts} // "allmapp,segdupmapp,nosegdupmapp");
my $PerlMod    = $cfg{slurm}{perl_module};
my $RMod       = $cfg{modules}{r}      // "R";
my @AfEdges    = split /,/, ($cfg{enrichment}{af_edges}
                             // "0,1e-05,1e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,2");

# Build the printed-form labels for the 9 AF bins so they match the worker
# output filenames (cf. __1_ComputeContingencyStatistics).
sub freq_label {
    my $v = shift;
    $v =~ s|/|~|; $v =~ s/\.//;
    return $v;
}
my @FreqIntervals;
if ($AfMode eq "unstratified") {
    # Single bin matching the submitter's unstratified mode (Min=0, Max=all)
    @FreqIntervals = ("0_all");
} else {
    for (my $i = 1; $i < @AfEdges; $i++) {
        push @FreqIntervals, freq_label($AfEdges[$i-1]) . "_" . freq_label($AfEdges[$i]);
    }
}

foreach my $Context (@Contexts) {
    my $RegionDir = "$StatsDir/stats_$Context";
    next unless -d $RegionDir;

    my $Tag           = "_All_${WindowSize}_${StepSize}_${VarListPrint}_StatsSummary";
    my $PosOut        = "$RegionDir/${Tag}.pos.txt";
    my $RScript       = "$RegionDir/${Tag}.R";
    my $ROutFile      = "$RegionDir/${Tag}.R.out.txt";
    my $JobFile       = "$RegionDir/${Tag}.job";
    my $eFile         = "$RegionDir/${Tag}.e";
    my $oFile         = "$RegionDir/${Tag}.o";

    my (%RepLengths, %GeneConvStats);

    opendir(my $SD, $RegionDir) or die "Can't open $RegionDir: $!\n";
    while (my $Subdir = readdir($SD)) {
        next unless $Subdir =~ /^\d+$/;
        next if $Subdir < 17;     # k=13,15 never used downstream
        my $KDir = "$RegionDir/$Subdir";
        opendir(my $KD, $KDir) or die "Can't open $KDir: $!\n";
        while (my $File = readdir($KD)) {
            next unless $File =~ /^(.*)_(\d+)_(\d+)_(.*)_(.*)_(.*)\.(\d+)\.(\d+)\.txt$/;
            my ($Chrom, $QW, $QS, $QVL, $QMin, $QMax, $IS, $IE)
                = ($1, $2, $3, $4, $5, $6, $7, $8);
            next unless ($QW == $WindowSize && $QS == $StepSize && $QVL eq $VarListPrint);

            my $fp = "$KDir/$File";
            open my $F, '<', $fp or die "Can't open $fp: $!\n";
            my $LineNr = 0;
            while (<$F>) {
                $LineNr++;
                next if $LineNr == 1;
                chomp;
                my @C = split /\t/, $_;
                my $Key = "${QMin}_${QMax}";
                $RepLengths{$Subdir} = undef;

                my @Cols = qw(
                    NrOfVarPosConvPos NrOfVarPosNoConvPos NrOfNoVarPosConvPos NrOfNoVarPosNoConvPos
                    NrOfVarMatchConv  NrOfVarMatchNoConv  NrOfNoVarMatchConv  NrOfNoVarMatchNoConv
                    NrOfConvPos       NrOfConvMatchVar    NrOfConvNoMatchVar  NrOfNoConvPos NrOfNoConvVar
                    ConcorVar         ConcorNoVar         DiscorVar           DiscorNoVar
                );
                # @C[0..2] are QueryChrom/WindowStart/WindowEnd; counts start at index 3
                for (my $ci = 0; $ci < @Cols; $ci++) {
                    $GeneConvStats{$Key}{$Subdir}{$Cols[$ci]} += $C[$ci + 3];
                }
            }
            close $F;
        }
        closedir($KD);
    }
    closedir($SD);

    # --- write three pos blocks ---------------------------------------------
    open my $OP, '>', $PosOut or die "Can't open $PosOut: $!\n";
    foreach my $FI (@FreqIntervals) { print $OP "\t$FI"; }
    print $OP "\n";

    # Block 1: PercDiff_Pos
    foreach my $K (sort { $a <=> $b } keys %RepLengths) {
        print $OP $K;
        foreach my $FI (@FreqIntervals) {
            my $NVP = $GeneConvStats{$FI}{$K}{NrOfVarPosNoConvPos}    // 0;
            my $NNN = $GeneConvStats{$FI}{$K}{NrOfNoVarPosNoConvPos}  // 0;
            foreach my $K2 (keys %RepLengths) {
                if ($K < $K2) {
                    $NVP -= $GeneConvStats{$FI}{$K2}{NrOfVarPosConvPos}   // 0;
                    $NNN -= $GeneConvStats{$FI}{$K2}{NrOfNoVarPosConvPos} // 0;
                }
            }
            my $VFConv    = $GeneConvStats{$FI}{$K}{NrOfVarPosConvPos}
                          / (($GeneConvStats{$FI}{$K}{NrOfVarPosConvPos}   // 0)
                           + ($GeneConvStats{$FI}{$K}{NrOfNoVarPosConvPos} // 0) || 1);
            my $VFNoConv  = $NVP / (($NVP + $NNN) || 1);
            my $PercDiff  = $VFNoConv ? sprintf("%.2f", ($VFConv - $VFNoConv) / $VFNoConv * 100) : "NA";
            print $OP "\t$PercDiff";
        }
        print $OP "\n";
    }
    print $OP "\n";

    # Block 2: ExcessPercPos
    foreach my $K (sort { $a <=> $b } keys %RepLengths) {
        print $OP $K;
        foreach my $FI (@FreqIntervals) {
            my $NVP = $GeneConvStats{$FI}{$K}{NrOfVarPosNoConvPos}    // 0;
            my $NNN = $GeneConvStats{$FI}{$K}{NrOfNoVarPosNoConvPos}  // 0;
            foreach my $K2 (keys %RepLengths) {
                if ($K < $K2) {
                    $NVP -= $GeneConvStats{$FI}{$K2}{NrOfVarPosConvPos}   // 0;
                    $NNN -= $GeneConvStats{$FI}{$K2}{NrOfNoVarPosConvPos} // 0;
                }
            }
            my $VFNoConv  = ($NVP + $NNN) ? $NVP / ($NVP + $NNN) : 0;
            my $ExpVar    = $VFNoConv * (($GeneConvStats{$FI}{$K}{NrOfVarPosConvPos}   // 0)
                                       + ($GeneConvStats{$FI}{$K}{NrOfNoVarPosConvPos} // 0));
            my $ExcessN   = ($GeneConvStats{$FI}{$K}{NrOfVarPosConvPos} // 0) - $ExpVar;
            my $Total     = (($GeneConvStats{$FI}{$K}{NrOfVarPosConvPos}   // 0)
                           + ($GeneConvStats{$FI}{$K}{NrOfVarPosNoConvPos} // 0)) || 1;
            print $OP "\t" . sprintf("%.3f", $ExcessN / $Total * 100);
        }
        print $OP "\n";
    }
    print $OP "\n";

    # Block 3: ExcessPercVar (concordance)
    foreach my $K (sort { $a <=> $b } keys %RepLengths) {
        print $OP $K;
        foreach my $FI (@FreqIntervals) {
            my $CV = $GeneConvStats{$FI}{$K}{ConcorVar}   // 0;
            my $CN = $GeneConvStats{$FI}{$K}{ConcorNoVar} // 0;
            my $DV = $GeneConvStats{$FI}{$K}{DiscorVar}   // 0;
            my $DN = $GeneConvStats{$FI}{$K}{DiscorNoVar} // 0;
            my $DiscorFreq = ($DV + $DN) ? $DV / ($DV + $DN) : 0;
            my $ConcorExp  = $DiscorFreq * ($CV + $CN);
            my $ConcorDiff = $CV - $ConcorExp;
            my $Total      = (($GeneConvStats{$FI}{$K}{NrOfVarPosConvPos}   // 0)
                            + ($GeneConvStats{$FI}{$K}{NrOfVarPosNoConvPos} // 0)) || 1;
            print $OP "\t" . sprintf("%.3f", $ConcorDiff / $Total * 100);
        }
        print $OP "\n";
    }
    print $OP "\n";
    close $OP;

    # --- write R driver -----------------------------------------------------
    open my $R, '>', $RScript or die "Can't open $RScript: $!\n";
    print $R <<'RHEAD';
library(pwr)
binom_test_approx <- function(x, n, p, alternative = "two.sided") {
  expected <- n * p
  se <- sqrt(n * p * (1 - p))
  z <- (x - expected) / se
  p_value <- if      (alternative == "two.sided") 2 * pnorm(-abs(z))
             else if (alternative == "greater")   pnorm(z, lower.tail = FALSE)
             else                                 pnorm(z)
  list(p.value = p_value)
}
RESULTS <- data.frame()
RHEAD

    foreach my $K (sort { $a <=> $b } keys %RepLengths) {
        foreach my $FI (@FreqIntervals) {
            my $G = $GeneConvStats{$FI}{$K};
            my ($A, $B, $C, $D) = ($G->{NrOfVarPosConvPos}    // 0,
                                   $G->{NrOfNoVarPosConvPos}  // 0,
                                   $G->{NrOfVarPosNoConvPos}  // 0,
                                   $G->{NrOfNoVarPosNoConvPos} // 0);
            # Set-difference correction
            foreach my $K2 (keys %RepLengths) {
                if ($K < $K2) {
                    $C -= $GeneConvStats{$FI}{$K2}{NrOfVarPosConvPos}   // 0;
                    $D -= $GeneConvStats{$FI}{$K2}{NrOfNoVarPosConvPos} // 0;
                }
            }
            my $Total = ($A + $B + $C + $D) || 1;
            my @Exp = (
                ($A + $B) * ($A + $C) / $Total,
                ($A + $B) * ($B + $D) / $Total,
                ($C + $D) * ($A + $C) / $Total,
                ($C + $D) * ($B + $D) / $Total,
            );
            my $Wsquared = 0;
            for my $cell ([$A, $Exp[0]], [$B, $Exp[1]], [$C, $Exp[2]], [$D, $Exp[3]]) {
                $Wsquared += ($cell->[0] - $cell->[1])**2 / ($cell->[1] || 1);
            }
            my $W = sqrt($Wsquared);

            my $NoConvVar  = $G->{NrOfNoConvVar}      // 0;
            my $NoConvPos  = $G->{NrOfNoConvPos}      || 1;
            my $ConvNoVar  = $G->{NrOfConvNoMatchVar} // 0;
            my $ConvMatch  = $G->{NrOfConvMatchVar}   // 0;
            my $ConvPos    = $G->{NrOfConvPos}        || 1;
            my $NonGcFreq  = (2 / 3) * $NoConvVar / $NoConvPos;
            my $GcFreq     = (1 / 3) * $NoConvVar / $NoConvPos;
            my $VarMatchNoConv = $G->{NrOfVarMatchNoConv} // 0;

            my $CV = $G->{ConcorVar} // 0; my $CN = $G->{ConcorNoVar} // 0;
            my $DV = $G->{DiscorVar} // 0; my $DN = $G->{DiscorNoVar} // 0;

            print $R <<"REVAL";
data_matrix <- matrix(c($A, $B, $C, $D), nrow = 2, byrow = TRUE)
chisq_test  <- chisq.test(data_matrix)
chisq_p     <- chisq_test\$p.value
chisq_p_log <- pchisq(chisq_test\$statistic, df = chisq_test\$parameter, lower.tail = FALSE, log.p = TRUE) / log(10)
power       <- pwr.chisq.test(w = $W, N = $Total, df = 1, sig.level = 0.005)\$power
binom_nongc <- binom_test_approx(x = $ConvNoVar, n = $ConvPos, p = $NonGcFreq, alternative = "two.sided")\$p.value
binom_gc    <- binom_test_approx(x = $ConvMatch, n = $ConvPos, p = $GcFreq,    alternative = "two.sided")\$p.value
vm_matrix   <- matrix(c($CV, $CN, $DV, $DN), nrow = 2, byrow = TRUE)
vm_test     <- chisq.test(vm_matrix)
vm_p        <- vm_test\$p.value
vm_p_log    <- pchisq(vm_test\$statistic, df = vm_test\$parameter, lower.tail = FALSE, log.p = TRUE) / log(10)
RESULTS     <- rbind(RESULTS, data.frame(
    RepLength="$K", FrequencyInterval="$FI",
    NrOfVarPosConvPos="$A", NrOfNoVarPosConvPos="$B",
    NrOfVarPosNoConvPos="$C", NrOfNoVarPosNoConvPos="$D",
    chisq_test_p=chisq_p, log_p_value_exponent=chisq_p_log,
    NrOfConvNoMatchVar="$ConvNoVar", NrOfVarMatchNoConv="$VarMatchNoConv",
    NrOfConvPos="$ConvPos", NonGcFrequency="$NonGcFreq",
    binom_test_nongc_p=binom_nongc,
    NrOfConvMatchVar="$ConvMatch", GcFrequency="$GcFreq", binom_test_gc_p=binom_gc,
    ConcorVar="$CV", ConcorNoVar="$CN", DiscorVar="$DV", DiscorNoVar="$DN",
    chisq_test_varmatch_p=vm_p, log_p_varmatch_value_exponent=vm_p_log,
    power=power))
REVAL
        }
    }

    print $R "write.table(RESULTS, file = \"$ROutFile\", sep = \"\\t\", "
           . "append = FALSE, col.names = TRUE, row.names = FALSE, quote = FALSE)\n";
    close $R;

    # --- SLURM job ----------------------------------------------------------
    open my $JOB, '>', $JobFile or die "Can't open $JobFile: $!\n";
    print $JOB <<"JOB";
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64400M
#SBATCH --error=$eFile
#SBATCH --output=$oFile
#SBATCH --time=03:00:00

module load $PerlMod
module load $RMod

Rscript $RScript
JOB
    close $JOB;
    chmod 0755, $JobFile;

    print "Wrote summary for $Context: $PosOut\n";
    print "Submit R driver: sbatch $JobFile\n";
}
