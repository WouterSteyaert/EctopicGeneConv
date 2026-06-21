#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
################																						 ################
################					 		   COMPARE TRACT LENGTH WITH GC DIFF						 ################
################																						 ################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

##################################################################################################
##################################################################################################
##################################################################################################
#######     		     	  	  		  Read Config File  	  		   			   	   #######
##################################################################################################
##################################################################################################
##################################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use IO::Zlib;
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

#######################################################################################
#######################################################################################
######     		     			 Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$GcDiffDir					= "$cfg{paths}{gcdiffs_dir}/";
my 	$ConvDir					= "$cfg{paths}{conv_dir}/";
my 	$TractLengthDir				= "$cfg{paths}{tractlengths_wgs_dir}/";
my 	$TractLengthSumDir			= "$cfg{paths}{tractlengths_sum_dir}/";
my 	$CHROM						= "10";
my 	$LineNr						= 0;
my 	$FileHandle					= new IO::Zlib;
my 	@REP_LENGTHS				= ("91_81", "81_71", "71_61", "61_51", "51_41", "41_31", "31_21", "21_19", "19_17");
my 	%GcDiffs					= ();
my 	%Variants					= ();
my 	%AllVariants				= ();

#######################################################################################
#######################################################################################
######     		     			  	  Get Options   	  		   			 	 ######
#######################################################################################
#######################################################################################

GetOptions ("CHROM=s"			=>	\$CHROM);

#######################################################################################
#######################################################################################
######     		     			 Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$TractLengthVsGcDiffFp		= "$cfg{paths}{tractlengths_sum_dir}/_TractLengthsVsGcChangesDiff.$CHROM.txt";

#######################################################################################
#######################################################################################
######     		     			    Read In GcDiffs   	  		   			 	 ######
#######################################################################################
#######################################################################################

foreach my $RepLength (@REP_LENGTHS){
	
	(my $MaxLength, my $MinLength) = split(m/_/, $RepLength);
	
	 my $GcDiffFp = $GcDiffDir . "$MinLength\_$MaxLength\/" . "$CHROM.changes.gcdiffs.sorted.bed.gz";
	 
	 print $GcDiffFp . "\n";
	
	if ($FileHandle->open($GcDiffFp, "rb")){
	
		while(<$FileHandle>){
			
			$_ =~ s/\n//g;	
			
			(undef, undef, my $Position, my $Info) = split(m/\t/, $_);
			(my $Change, my $Recurrence) = split(m/;/, $Info);
			
			$GcDiffs{$Position}{$Change}{$MinLength} = undef;
			
		}
	}
	$FileHandle->close;	
}

# Add 91

my $ConvFp = $ConvDir . "91/$CHROM.repsum.geneconv.sorted.bed.gz";

if ($FileHandle->open($ConvFp, "rb")){
	
	while(<$FileHandle>){
		
		$_ =~ s/\n//g;	
		
		(undef, undef, my $Position, my $Info) = split(m/\t/, $_);
		(my $Change, my $Recurrence) = split(m/;/, $Info);
		
		$GcDiffs{$Position}{$Change}{"91"} = undef;
		
	}
}
$FileHandle->close;	

#######################################################################################
#######################################################################################
######     		     			    Read In Variants   	  		   			 	 ######
#######################################################################################
#######################################################################################

print "Read In Variants\n";

opendir(SUM, $TractLengthSumDir) or die ("Can't open $TractLengthSumDir\n");
while(my $File = readdir(SUM)){
	
	if ($File =~ m/^_TractLengthsVsFreq\.\d+\.scatter$/){

		my $FilePath = $TractLengthSumDir . $File;
		
		open V, "$FilePath" or die ("Can't open $FilePath\n");
		while(<V>){
			
			if ($LineNr){
				
				$_ =~ s/\n//g;
				
				(my $Variant, my $AlleleFreq, my $TractLength, my $TractAnalysisLength) = split(m/\t/, $_);
				(my $Chrom, my $VarPosition, undef, my $NucChange) = split(m/_/, $Variant);
				(my $RefNuc, my $AltNuc) = split(m/\//, $NucChange);
				
				if ($Chrom eq $CHROM){
					
					$AllVariants{$_} = undef;
					
					if (exists $GcDiffs{$VarPosition}{"$RefNuc>$AltNuc"}){
						
						my 	@Lengths 	= sort {$a <=> $b} keys %{$GcDiffs{$VarPosition}{"$RefNuc>$AltNuc"}};
						my 	$MaxLength	= $Lengths[-1];
					
						$Variants{$_} = $MaxLength;
					}
				}
			}
			$LineNr++;
		}
		close V;
	}
}
closedir(SUM);


#######################################################################################
#######################################################################################
######     		     			    	Write Out   	  		   			 	 ######
#######################################################################################
#######################################################################################

open O, ">$TractLengthVsGcDiffFp" or die ("Can't open $TractLengthVsGcDiffFp\n");
print O "Variant\tFrequency\tTractLength\tTractAnalysisLength\tMaxGcLength\n";

foreach my $VariantLine (sort keys %AllVariants){
	
	if ($Variants{$VariantLine}){
		
		print O $VariantLine . "\t" . $Variants{$VariantLine} . "\n";
	}
	else {
		
		print O $VariantLine . "\t" . "/" . "\n";
	}
}

close O;




