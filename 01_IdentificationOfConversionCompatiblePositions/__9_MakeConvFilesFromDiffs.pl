##############################################################################
##############################################################################
##############################################################################
########                                                            ##########
########               MAKE CONV FILES FROM DIFFS                   ##########
########                                                            ##########
##############################################################################
##############################################################################
##############################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use IO::Zlib;
use File::Path qw(make_path);
use File::Basename qw(dirname);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

# Load project config (consumes --ConfigFile= and --ProjectRoot= from @ARGV).
my %cfg = load_config();

my $TEMPDir       = "$cfg{paths}{data_root}/temp/";
my $Chromosome    = 22;
my @RepeatLengths = split /,/, $cfg{repeat_lengths}{values};
# Drop the longest k (no longer cumulative file to build for it)
pop @RepeatLengths;

make_path($TEMPDir) unless -d $TEMPDir;

##############################################################################
##############################################################################
################   				   Get Options   			  ################
##############################################################################
##############################################################################

GetOptions("Chromosome=s" => \$Chromosome,
           "CHROM=s"      => \$Chromosome);   # alias for consistency with submitters

##############################################################################
##############################################################################
################   		   Iterate Over RepeatLengths   	  ################
##############################################################################
##############################################################################

my $CONVDIR    = $cfg{paths}{conv_dir};
my $GCDIFFSDIR = $cfg{paths}{gcdiffs_dir};

for (my $I = $#RepeatLengths; $I >=0; $I--){

    my $RepLength 		= $RepeatLengths[$I];
    my $ConvRepLength 	= 0;

	if ($RepLength >=21)	{$ConvRepLength = $RepLength + 10;}
	else 					{$ConvRepLength = $RepLength + 2;}

    my $ConvFilePath  	= "$CONVDIR/$ConvRepLength/$Chromosome.repsum.geneconv.sorted.bed.gz";
	my $GcdDiffFilePath = "$GCDIFFSDIR/${RepLength}_${ConvRepLength}/$Chromosome.changes.gcdiffs.sorted.bed.gz";
	my $OutputFilePath  = "$GCDIFFSDIR/${RepLength}_${ConvRepLength}/$Chromosome.repsum.geneconv.sorted.bed.gz";
	my $CheckFilePath 	= "$CONVDIR/$RepLength/$Chromosome.repsum.geneconv.sorted.bed.gz";

    my $cmd = qq{LC_ALL=C sort -m -k1,1 -k2,2n -T "$TEMPDir" --parallel=1 <(zcat "$ConvFilePath") <(zcat "$GcdDiffFilePath") | bgzip -\@1 > "$OutputFilePath"};
	system("bash", "-c", $cmd) == 0 or die "Failed: $cmd\n";

    
    print "tabix -p bed $OutputFilePath\n";
    system("tabix -p bed $OutputFilePath");
	
	system("diff -s <(zcat $OutputFilePath | cut -f1-3)  <(zcat $CheckFilePath | cut -f1-3)");
}
