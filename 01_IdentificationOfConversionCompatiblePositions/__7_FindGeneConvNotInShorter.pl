#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
################																						 ################
################					 FIND GENE CONVERSIONS NOT IN SHORTER REPEAT LENGTH					 ################
################																						 ################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

#######################################################################################
#######################################################################################
######     		     			  	   HPC Modules  	  		   				 ######
#######################################################################################
#######################################################################################

# module load Perl/5.34.0-GCCcore-11.2.0

#######################################################################################
#######################################################################################
######     		     			  Load Other Libraries  	  		   			 ######
#######################################################################################
#######################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use IO::Zlib;
use List::Util qw(min max);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

#######################################################################################
#######################################################################################
######     		     			Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$REP_LENGTHS			= "21_19";
my 	$CHROM					= "22";
my 	$START					= 1;
my 	$END					= 1;
my 	$FileHandle				= new IO::Zlib;

#######################################################################################
#######################################################################################
######     		     			  	  Get Options   	  		   			 	 ######
#######################################################################################
#######################################################################################

GetOptions ("CHROM=s"				=>	\$CHROM,
			"START=i"				=> 	\$START,
			"END=i"					=> 	\$END,
			"REP_LENGTHS=s"			=>	\$REP_LENGTHS);

#######################################################################################
#######################################################################################
######     		     			Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my @REP_LENGTHS					= split(m/\_/, $REP_LENGTHS);
my $Short						= min(@REP_LENGTHS);
my $Long						= max(@REP_LENGTHS);
my $ShortGCs					= "$cfg{paths}{conv_dir}/$Short/$CHROM.repsum.geneconv.sorted.bed.gz";
my $LongGCs						= "$cfg{paths}{conv_dir}/$Long/$CHROM.repsum.geneconv.sorted.bed.gz";
my $DiffSubDir					= "$cfg{paths}{gcdiffs_dir}/${Short}_${Long}";
my $OutputPositions				= "$DiffSubDir/$CHROM.$START.$END.positions.gcdiffs.txt";
my $OutputChanges				= "$DiffSubDir/$CHROM.$START.$END.changes.gcdiffs.txt";
my %ShortGCs					= ();
my %LongGCs						= ();
my %Differences					= ();

unless (-d $DiffSubDir) { system("mkdir -p $DiffSubDir"); }

#######################################################################################
#######################################################################################
######     		     	  		Read In Short Gene Conversions  	  		 	 ######
#######################################################################################
#######################################################################################

print "Read in " . $ShortGCs . "\n";

if ($FileHandle->open($ShortGCs, "rb")){
	
	while(<$FileHandle>){
		
		$_ =~ s/\n//g;		
		
		(undef, undef, my $Position, my $Info) = split(m/\t/, $_);
		(my $Change, my $Recurrence) = split(m/;/, $Info);
		
		if ($Position >=$START && $Position < $END){
			
			$ShortGCs{$Position}{$Change} = $Recurrence;
		}
	}
}
$FileHandle->close;

#######################################################################################
#######################################################################################
######     		     	  		 Read In Long Gene Conversions  	  		 	 ######
#######################################################################################
#######################################################################################

print "Read in " . $LongGCs . "\n";

if ($FileHandle->open($LongGCs, "rb")){
	
	while(<$FileHandle>){
		
		$_ =~ s/\n//g;
		
		(undef, undef, my $Position, my $Info) = split(m/\t/, $_);
		(my $Change, my $Recurrence) = split(m/;/, $Info);
		
		if ($Position >=$START && $Position < $END){
		
			$LongGCs{$Position}{$Change} = $Recurrence;
		}
	}
}
$FileHandle->close;

#######################################################################################
#######################################################################################
######     		     	  		 		Differences  	  		 	 			 ######
#######################################################################################
#######################################################################################

open C, ">$OutputChanges" or die ("Can't open $OutputChanges\n");
open P, ">$OutputPositions" or die ("Can't open $OutputPositions\n");

foreach my $Position (sort {$a <=> $b} keys %ShortGCs){
	
	my $BedStart = $Position-1;
	
	if (!exists $LongGCs{$Position}){
		
		foreach my $Change (sort keys %{$ShortGCs{$Position}}){
		
			print P $CHROM . "\t" . $BedStart . "\t" . $Position . "\n";
			print C $CHROM . "\t" . $BedStart . "\t" . $Position . "\t" . $Change . ";" . $ShortGCs{$Position}{$Change} . "\n";
		}
	}
	else {
	
		foreach my $Change (sort keys %{$ShortGCs{$Position}}){
			
			if (!exists $LongGCs{$Position}{$Change}){
				
				print C $CHROM . "\t" . $BedStart . "\t" . $Position . "\t" . $Change . ";" . $ShortGCs{$Position}{$Change} . "\n";
			}
		}
	}
}

close P;
close C;

if (-f "$OutputChanges.gz"){system("rm $OutputChanges.gz");}
if (-f "$OutputPositions.gz"){system("rm $OutputPositions.gz");}

system("gzip $OutputChanges");
system("gzip $OutputPositions");

print "The script has finished with success\n"; 


