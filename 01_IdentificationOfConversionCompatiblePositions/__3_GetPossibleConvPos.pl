#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
################																						 ################
################					 		 	GET POSSIBLE CONVERSION POSITIONS						 ################
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
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

#######################################################################################
#######################################################################################
######     		     			Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$REP_LEN				= 51;
my 	$FileHandle				= new IO::Zlib;
my 	$Seed					= "AAAAA";
my 	@Bases					= ("A", "C", "G", "T");
my 	%Repeats				= ();

#######################################################################################
#######################################################################################
######     		     			  	  Get Options   	  		   			 	 ######
#######################################################################################
#######################################################################################

GetOptions ("SEED=s"				=>	\$Seed,
			"REP_LEN=i"				=>	\$REP_LEN);

#######################################################################################
#######################################################################################
######     		     			  Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$RepSumDir					= "$cfg{paths}{repsum_dir}/$REP_LEN/";
my 	$ConvDir					= "$cfg{paths}{conv_dir}/$REP_LEN/";
my 	$RepSumFilePath				= $RepSumDir . $Seed . ".repsum.txt.gz";
my 	$RepGeneConvCoordFilePath	= $ConvDir . $Seed . ".repsum.geneconv.coords.txt";
my 	$RepGeneConvCountFilePath	= $ConvDir . $Seed . ".repsum.geneconv.counts.txt";
my 	$FlankLength				= ($REP_LEN-1)/2;

unless (-d $ConvDir){system("mkdir -p $ConvDir");}

#######################################################################################
#######################################################################################
######     		     			  	 Read In Repeats  	  		   				 ######
#######################################################################################
#######################################################################################

print "Read In Repeats\n";

if ($FileHandle->open($RepSumFilePath, "rb")){
	
	while(<$FileHandle>){
		
		$_ =~ s/\n//g;
	
		(my $Repeat, my $Coords_0, my $Coords_1) = split(m/\t/, $_);
		
		if ($Coords_0 ne "/"){
			
			my @Coords_0 = split(m/\,/, $Coords_0);
			
			foreach my $Coord_0 (@Coords_0){
			
				(my $Chrom, my $Position) = split(m/\:/, $Coord_0);
				
				$Repeats{$Repeat}{0}{$Chrom}{$Position} = undef;
			}
		}
		if ($Coords_1 ne "/"){
			
			my @Coords_1 = split(m/\,/, $Coords_1);
			
			foreach my $Coord_1 (@Coords_1){
			
				(my $Chrom, my $Position) = split(m/\:/, $Coord_1);
				
				$Repeats{$Repeat}{1}{$Chrom}{$Position} = undef;
			}
		}
	}
}

$FileHandle->close();

#######################################################################################
#######################################################################################
######     		     			  	Iterate Over Repeats  	  		   			 ######
#######################################################################################
#######################################################################################

open C, ">$RepGeneConvCountFilePath" or die ("Can't open $RepGeneConvCountFilePath\n");
open O, ">$RepGeneConvCoordFilePath" or die ("Can't open $RepGeneConvCoordFilePath\n");

foreach my $Repeat (sort keys %Repeats){

	my 	$CentralBase 		= substr($Repeat,$FlankLength,1);
	my 	$ComplCentralBase	= $CentralBase;
		$ComplCentralBase	=~ tr/ATGCatgc/TACGtacg/;
	my 	$RepCoords_0		= "";
	my 	$RepCoords_1		= "";
	my 	$CompRepCoords_0 	= "";
	my 	$CompRepCoords_1 	= "";
	my 	$RepeatInstances	= 0;
	my 	%ConvInstances		= ();
	
	if (exists $Repeats{$Repeat}{0}){
	
		foreach my $Chrom (keys %{$Repeats{$Repeat}{0}}){
			foreach my $Position (keys %{$Repeats{$Repeat}{0}{$Chrom}}){
				
				$Position+=$FlankLength;
				
				$RepCoords_0 .= "$CentralBase=$Chrom:$Position,";
				
				$RepeatInstances++;
			}
		}
		
		$RepCoords_0 =~ s/,$//g;
	}
	
	if (exists $Repeats{$Repeat}{1}){
	
		foreach my $Chrom (keys %{$Repeats{$Repeat}{1}}){
			foreach my $Position (keys %{$Repeats{$Repeat}{1}{$Chrom}}){
				
				$Position+=$FlankLength;
				
				$RepCoords_1 .= "$ComplCentralBase=$Chrom:$Position,";
				
				$RepeatInstances++;
			}
		}
		
		$RepCoords_1 =~ s/,$//g;
	}
	
	foreach my $Base (@Bases){
		
		if ($Base ne $CentralBase){

			my $ComplementRepeat = substr($Repeat,0,$FlankLength) . $Base . substr($Repeat,($FlankLength+1),$FlankLength);

			if (exists $Repeats{$ComplementRepeat} && 
				exists $Repeats{$ComplementRepeat}{0}){
				
				foreach my $Chrom (keys %{$Repeats{$ComplementRepeat}{0}}){
					foreach my $Position (keys %{$Repeats{$ComplementRepeat}{0}{$Chrom}}){
						
						$Position+=$FlankLength;
						
						$CompRepCoords_0 .= "$Base=$Chrom:$Position,";
						
						if ($ConvInstances{$Base})	{$ConvInstances{$Base}++;}
						else 						{$ConvInstances{$Base} = 1;}
					}
				}
			}
			
			if (exists $Repeats{$ComplementRepeat} &&
				exists $Repeats{$ComplementRepeat}{1}){
				
				foreach my $Chrom (keys %{$Repeats{$ComplementRepeat}{1}}){
					foreach my $Position (keys %{$Repeats{$ComplementRepeat}{1}{$Chrom}}){
						
						$Position+=$FlankLength;

						my 	$ComplBase 	= $Base;
							$ComplBase	=~ tr/ATGCatgc/TACGtacg/; # this is the base in the refseq
							
						$CompRepCoords_1 .= "$ComplBase=$Chrom:$Position,";
						
						if ($ConvInstances{$Base})	{$ConvInstances{$Base}++;}
						else 						{$ConvInstances{$Base} = 1;}
					}
				}
			}	
		}
	}
	
	if ($CompRepCoords_0  ne "" || $CompRepCoords_1 ne ""){
		
		$CompRepCoords_0 =~ s/,$//g;
		$CompRepCoords_1 =~ s/,$//g;
		
		if ($RepCoords_0 eq "")			{$RepCoords_0 = "/";}
		if ($RepCoords_1 eq "")			{$RepCoords_1 = "/";}
		if ($CompRepCoords_0 eq "")		{$CompRepCoords_0 = "/";}
		if ($CompRepCoords_1 eq "")		{$CompRepCoords_1 = "/";}
				
		print O $Repeat . "\t" . $RepCoords_0 . "\t" . $RepCoords_1 . "\t" . $CompRepCoords_0 . "\t" . $CompRepCoords_1 . "\n";
		
		foreach my $ConvBase (sort keys %ConvInstances){
			
			print C $Repeat . "\t" . $CentralBase  . "\t" . $RepeatInstances . "\t" . $ConvBase . "\t" . $ConvInstances{$ConvBase} . "\n";
		}
	}	
}

close O;
close C;

if (-f "$RepGeneConvCoordFilePath.gz"){system("rm $RepGeneConvCoordFilePath.gz");}
if (-f "$RepGeneConvCountFilePath.gz"){system("rm $RepGeneConvCountFilePath.gz");}

system("gzip $RepGeneConvCoordFilePath");
system("gzip $RepGeneConvCountFilePath");


print "The script has finished with success\n";








