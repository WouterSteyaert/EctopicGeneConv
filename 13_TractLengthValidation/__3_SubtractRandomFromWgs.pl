#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
################																						 ################
################								   SUBTRACT RANDOM FROM WGS								 ################
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
use List::Util qw(min max);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

#######################################################################################
#######################################################################################
######     		     			 Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$TractLengthSumDirWgs			= "$cfg{paths}{tractlengths_sum_dir}/";
my 	$TractLengthSumDirRandom		= "$cfg{paths}{data_root}/tractlengthsSumRandom3/";
my 	$OutputFilePath 				= $TractLengthSumDirWgs . "_TractLengthAnalysis.txt";
my 	@Frequencies					= split /,/, $cfg{enrichment}{af_edges};
my 	$BinWidth 						= 5;
my 	$LineNr 						= 0;
my 	%RandomDistribution				= ();
my 	%WgsDistribution				= ();
my 	%ExistingBins 					= ();

#######################################################################################
#######################################################################################
######     		     			   Read In Wgs Scatters 	  		   			 ######
#######################################################################################
#######################################################################################

opendir (WGS, $TractLengthSumDirWgs) or die ("Can't open $TractLengthSumDirWgs");
while(my $File = readdir(WGS)){
	
	$LineNr = 0;
	
	if ($File =~ m/.*\.scatter$/){
		
		my $FilePath = $TractLengthSumDirWgs . $File;
		
		open F, "$FilePath" or die ("Can't open $FilePath\n");
		while(<F>){
			
			if ($LineNr){
			
				$_ =~ s/\n//g;
			
				my 	@LineValues 	= split(m/\t/, $_);
				my 	$Frequency		= $LineValues[1];
				my 	$TractLength	= $LineValues[2]; # raw BLAST tract length (asymmetric); minimum tract length 21 bp filters out short random matches
				next if ($TractLength < 21);
				my 	$Bin 			= int($TractLength/$BinWidth)+1;
				my 	$FreqInterval	= "";
				
				$ExistingBins{$Bin} = undef;
				
				for (my $I = 0; $I < $#Frequencies; $I++) {
					
					if ($Frequency >= $Frequencies[$I] && $Frequency < $Frequencies[$I + 1]) {
				
						$FreqInterval = "$Frequencies[$I]\_$Frequencies[$I + 1]";
					}
				}
				
				if ($WgsDistribution{$FreqInterval} 		&&
					$WgsDistribution{$FreqInterval}{$Bin}){
				
					$WgsDistribution{$FreqInterval}{$Bin}++;
				}
				else {
					$WgsDistribution{$FreqInterval}{$Bin}=1;
				}
			}
			$LineNr++;
		}
		close F;
	}
}
closedir(WGS);

#print Dumper %ExistingBins;
#die($!);

#######################################################################################
#######################################################################################
######     		     			   Read In Random Scatters 	  		   			 ######
#######################################################################################
#######################################################################################

opendir (RAND, $TractLengthSumDirRandom) or die ("Can't open $TractLengthSumDirRandom");
while(my $File = readdir(RAND)){
	
	$LineNr = 0;
	
	if ($File =~ m/.*\.scatter$/){
		
		my $FilePath = $TractLengthSumDirRandom . $File;
		
		open F, "$FilePath" or die ("Can't open $FilePath\n");
		while(<F>){
			
			if ($LineNr){
			
				$_ =~ s/\n//g;
			
				my 	@LineValues 	= split(m/\t/, $_);
				my 	$TractLength	= $LineValues[2]; # raw BLAST tract length (asymmetric); minimum tract length 21 bp filters out short random matches
				next if ($TractLength < 21);
				my 	$Bin 			= int($TractLength/$BinWidth)+1;
				
				$ExistingBins{$Bin} = undef;
				
				if ($RandomDistribution{$Bin})	{$RandomDistribution{$Bin}++;}
				else 							{$RandomDistribution{$Bin}=1;}
			}
			$LineNr++;
		}
		close F;
	}
}
closedir(RAND);

#######################################################################################
#######################################################################################
######     		     			   		  Subtract 	  		 	  				 ######
#######################################################################################
#######################################################################################

open O, ">$OutputFilePath" or die ("Can't open $OutputFilePath\n");
print O "FrequencyInterval\tAbsDifference\tRelDifference\tSum_Random\tSum_GnomAD\n";

for (my $I = 0; $I < $#Frequencies; $I++) {
	
	my 	$FrequencyInterval 	= "$Frequencies[$I]\_$Frequencies[$I + 1]";
	my 	$Difference 		= 0;
	my 	$RandomSum			= 0;
	my 	$WgsSum				= 0;
	
	foreach my $Bin (sort {$a <=> $b} keys %ExistingBins){
		
		my $RandomNr 	= 0;
		my $WgsNr		= 0;
		
		if ($WgsDistribution{$FrequencyInterval}{$Bin})	{$WgsNr = $WgsDistribution{$FrequencyInterval}{$Bin};}
		if ($RandomDistribution{$Bin})					{$RandomNr = $RandomDistribution{$Bin};}
		
		$RandomSum+=$RandomNr;
		$WgsSum+=$WgsNr;
		
		$Difference += ($WgsNr-$RandomNr);
		
		#print "\t" . $WgsNr . "\t" . $RandomNr . "\t" . $Difference . "\n";
	}
	
	my $RelDifference = $Difference/50;
	
	printf "%-20s %-15.1f %-15.1f %-15.1f %-15.1f\n",
    $FrequencyInterval, $Difference, $RelDifference, $RandomSum, $WgsSum;
	
	print O "$FrequencyInterval\t$Difference\t$RelDifference\t$RandomSum\t$WgsSum\n";
}

close O;