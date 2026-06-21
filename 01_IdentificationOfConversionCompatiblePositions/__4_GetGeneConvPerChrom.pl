#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
################																						 ################
################					 		 	  GET GENE CONVERSIONS PER CHROM						 ################
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

my 	$CHROM					= "1";
my 	$REP_LEN				= 11;
my 	$FileHandle				= new IO::Zlib;
my 	$SEED					= "A";
my 	$SUBSEED				= "";

#######################################################################################
#######################################################################################
######     		     			  	    Get Options   	  		   			 	 ######
#######################################################################################
#######################################################################################

GetOptions ("CHROM=s"		=>	\$CHROM,
			"REP_LEN=i"		=>	\$REP_LEN,
			"SEED=s"		=>	\$SEED,
			"SUBSEED=s"		=>	\$SUBSEED);

#######################################################################################
#######################################################################################
######     		     			Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$RepSumDir				= "$cfg{paths}{conv_dir}/$REP_LEN/";
my 	$SeedQM					= quotemeta($SEED);
my 	$SeedSubseedQM			= quotemeta("$SEED$SUBSEED");
my 	%ConvCoords				= ();
my 	$RepGeneConvCoordFP		= "";

if ($SUBSEED ne ""){
	
	$RepGeneConvCoordFP		= $RepSumDir . "$CHROM.$SEED.$SUBSEED.repsum.geneconv.txt";
}
else {
	
	$RepGeneConvCoordFP		= $RepSumDir . "$CHROM.$SEED.repsum.geneconv.txt";
}


#######################################################################################
#######################################################################################
######     		     			  	 Read In Repeats  	  		   				 ######
#######################################################################################
#######################################################################################

print "Read In Coords\n";

opendir (DIR, $RepSumDir) or die ("Can't open $RepSumDir");
while (my $File = readdir(DIR)){
	
	if ($File =~ m/^$SeedQM.*\.repsum\.geneconv\.coords\.txt\.gz$/){
		
		my $FilePath = $RepSumDir . $File;
		
		if ($FileHandle->open($FilePath, "rb")){
	
			while(<$FileHandle>){
				
				$_ =~ s/\n//g;
							
				(my $Repeat, 
				 my $RepCoords_0, 
				 my $RepCoords_1, 
				 my $CompRepCoords_0, 
				 my $CompRepCoords_1) 
				= split(m/\t/, $_);
				
				if ($SUBSEED eq "/" || $Repeat =~ m/^$SeedSubseedQM/){
					
					my @RepCoords_0 		= ();
					my @RepCoords_1 		= ();
					my @CompRepCoords_0		= ();
					my @CompRepCoords_1		= ();
					
					if ($RepCoords_0 ne "/")			{@RepCoords_0 			= split(m/\,/, $RepCoords_0);}
					if ($RepCoords_1 ne "/")			{@RepCoords_1 			= split(m/\,/, $RepCoords_1);}
					if ($CompRepCoords_0 ne "/")		{@CompRepCoords_0 		= split(m/\,/, $CompRepCoords_0);}
					if ($CompRepCoords_1 ne "/")		{@CompRepCoords_1 		= split(m/\,/, $CompRepCoords_1);}
					
					foreach my $RepCoord_0 (@RepCoords_0){
						
						$RepCoord_0 =~ m/^([A|C|G|T])\=(.*):(\d+)$/;
						
						my $Base		= $1; #base in ref
						my $Chrom 		= $2;
						my $Position 	= $3;
						
						if ($Chrom eq $CHROM){

							foreach my $CompRepCoord_0 (@CompRepCoords_0){ # no revcom needed
								
								$CompRepCoord_0 =~ m/^([A|C|G|T])\=(.*):(\d+)$/;
						
								my 	$CompBase = $1; # base in ref
								
								if ($ConvCoords{$Position}{"$Base>$CompBase"})	{$ConvCoords{$Position}{"$Base>$CompBase"}++;}
								else 											{$ConvCoords{$Position}{"$Base>$CompBase"}=1;}
							}

							foreach my $CompRepCoord_1 (@CompRepCoords_1){ # take revcom
				
								$CompRepCoord_1 =~ m/^([A|C|G|T])\=(.*):(\d+)$/;
						
								my 	$CompBase = $1;
									$CompBase =~ tr/ATGCatgc/TACGtacg/;
									
								if ($ConvCoords{$Position}{"$Base>$CompBase"})	{$ConvCoords{$Position}{"$Base>$CompBase"}++;}
								else 											{$ConvCoords{$Position}{"$Base>$CompBase"}=1;}
							}	
						}
					}
					
					# foreach my $RepCoords_1 (@RepCoords_1){
						
						# $RepCoords_1 =~ m/^([A|C|G|T])\=(.*):(\d+)$/;
						
						# my $Base		= $1; #base in ref
						# my $Chrom 		= $2;
						# my $Position 	= $3;
						
						# if ($Chrom eq $CHROM){
							
							# if ($Position == 14657){
								# print "$FilePath\n";
								# print "$_\n";
							# }
							
							# foreach my $CompRepCoord_0 (@CompRepCoords_0){ # revcom needed
								
								# $CompRepCoord_0 =~ m/^([A|C|G|T])\=(.*):(\d+)$/;
						
								# my 	$CompBase = $1;
									# $CompBase =~ tr/ATGCatgc/TACGtacg/;
									
								# if ($ConvCoords{$Position}{"$Base>$CompBase"})	{$ConvCoords{$Position}{"$Base>$CompBase"}++;}
								# else 											{$ConvCoords{$Position}{"$Base>$CompBase"}=1;}
							# }

							# foreach my $CompRepCoord_1 (@CompRepCoords_1){ # no revcom
				
								# $CompRepCoord_1 =~ m/^([A|C|G|T])\=(.*):(\d+)$/;
						
								# my 	$CompBase = $1;
								
								# if ($ConvCoords{$Position}{"$Base>$CompBase"})	{$ConvCoords{$Position}{"$Base>$CompBase"}++;}
								# else 											{$ConvCoords{$Position}{"$Base>$CompBase"}=1;}	
							# }
						# }
					# }					
				}
			}
		}

		$FileHandle->close();
	}
}
closedir(DIR);

open O, ">$RepGeneConvCoordFP" or die ("Can't open $RepGeneConvCoordFP\n");

foreach my $Position (sort {$a <=> $b} keys %ConvCoords){
	foreach my $Change (sort keys %{$ConvCoords{$Position}}){
		
		my $Recurrence = $ConvCoords{$Position}{$Change};
		
		print O $Position . "\t" . $Change . "\t" . $Recurrence . "\n";
	}
}

close O;

if (-f "$RepGeneConvCoordFP.gz"){system("rm $RepGeneConvCoordFP.gz");}

system("gzip $RepGeneConvCoordFP");

print "The script has finished with success\n";

