#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
################																						 ################
################							FETCH SEQUENCE REPEATS FROM FASTA						 	 ################
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
######     		     			  	Load Other Libraries  	  		   			 ######
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

# Load project config (consumes --ConfigFile= and --ProjectRoot= from @ARGV).
my %cfg = load_config();

#######################################################################################
#######################################################################################
######     		     			Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$CHROM								= "22";
my 	$QUERY_START						= 10000000;
my 	$QUERY_END							= 19999999;
my 	$REP_LEN							= 41;

#######################################################################################
#######################################################################################
######     		     			  	  Get Options   	  		   			 	 ######
#######################################################################################
#######################################################################################

GetOptions ("CHROM=s"				=>	\$CHROM,
			"QUERY_START=i"			=>	\$QUERY_START,
			"QUERY_END=i"			=>	\$QUERY_END,
			"REP_LEN=i"				=>	\$REP_LEN);

#######################################################################################
#######################################################################################
######     		     			Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$FastaFilePath						= sprintf($cfg{paths}{reference_fasta_pattern}, $CHROM);
my 	$LineNr								= 1;
my 	$First								= 1;
my	$Position							= 1;
my 	$Rest								= "";
my	$Slice 								= "";
my 	$Segment							= "";
my 	$RevComSlice						= "";
my 	$FileHandle							= new IO::Zlib;
my 	$RepeatFolder						= "$cfg{paths}{repeats_dir}/$REP_LEN/";
my 	$PositionFilePath					= $RepeatFolder . $CHROM . ".rep.$QUERY_START\_$QUERY_END.pos";
my 	%RepeatPositions					= ();

unless (-d $RepeatFolder){system("mkdir -p $RepeatFolder");};
			
#######################################################################################
#######################################################################################
######     		     			 	 Loop Over Chromosome 	  		   		 	 ######
#######################################################################################
#######################################################################################	

if ($FileHandle->open($FastaFilePath, "rb")){
	
	while(<$FileHandle>){
		
		if ($_ !~ m/>/){

			$_ =~ s/\n//g;
			
			if ($First == 0) {
				$Rest .= $_;
			}
			if 		(length($Segment) <= $REP_LEN){									# enlarge segment
				$Segment .= $_;
			}
			elsif 	($Rest eq "" && $First == 1){											# split segment for first time
			
				$Slice 		 	= substr($Segment, 0, $REP_LEN);
				$Rest	 		= substr($Segment, $REP_LEN, length($Segment)-$REP_LEN);
				$Rest			.= $_;
				$First	 		= 0;
				
				$Slice 	 		= uc ($Slice);
				$RevComSlice	= revcom($Slice);
				
				if ($Slice !~ m/N/) {
					
					if ($Position >= $QUERY_START && $Position <= $QUERY_END){
						
						# Positions
						
						if 		($RepeatPositions{$Slice} && 
								 $RepeatPositions{$Slice}{0})			{$RepeatPositions{$Slice}{0}.= "$Position\,";}			
						else 											{$RepeatPositions{$Slice}{0} = "$Position,";}
						
						if 		($RepeatPositions{$RevComSlice} &&
								 $RepeatPositions{$RevComSlice}{1})		{$RepeatPositions{$RevComSlice}{1}.= "$Position\,";}		
						else 											{$RepeatPositions{$RevComSlice}{1} = "$Position,";}
					}
				}
			}		
			if 		($Rest ne ""){
			
				do {
																	# glide over segment												
					$Position++;
				
					$Slice 	 		= 	substr($Slice, 1, $REP_LEN-1);
					$Slice 			.= 	substr($Rest, 0, 1);
					$Slice 			= 	uc ($Slice);
					$RevComSlice	= revcom($Slice);

					if ($Slice !~ m/N/) {
						
						if ($Position >= $QUERY_START && $Position <= $QUERY_END){
							
							# Positions
						
							if 		($RepeatPositions{$Slice} && 
									 $RepeatPositions{$Slice}{0})			{$RepeatPositions{$Slice}{0}.= "$Position\,";}			
							else 											{$RepeatPositions{$Slice}{0} = "$Position,";}
							
							if 		($RepeatPositions{$RevComSlice} &&
									 $RepeatPositions{$RevComSlice}{1})		{$RepeatPositions{$RevComSlice}{1}.= "$Position\,";}	
							else 											{$RepeatPositions{$RevComSlice}{1} = "$Position,";}
						}
					}
					
					if (length($Rest)==1){
						$Rest = "";
					}
					else {
						$Rest = substr($Rest, 1, length($Rest)-1);
					}
				}
				until ($Rest eq "");
			}
		}
		
		$LineNr++;
		
		if ($LineNr%25000==0){
			
			print $LineNr . "\n";
		}
	}
}

$FileHandle->close();

#######################################################################################
#######################################################################################
######     		     			 	 Write Out Repeats  	  		   		 	 ######
#######################################################################################
#######################################################################################

open P, ">$PositionFilePath" or die ("Can't open $PositionFilePath\n");

foreach my $Repeat (keys %RepeatPositions){

	print P "$Repeat";
	
	if (exists $RepeatPositions{$Repeat}{0}){
		
		my 	$Positions 			= $RepeatPositions{$Repeat}{0};
			$Positions			=~ s/,$//g;
			
		print P "\t$Positions";
	}
	else {
		
		print P "\t/";
	}
	
	if (exists $RepeatPositions{$Repeat}{1}){
		
		my 	$Positions 			= $RepeatPositions{$Repeat}{1};
			$Positions			=~ s/,$//g;
			
		print P "\t$Positions";
	}
	else {
		
		print P "\t/";
	}
	
	print P "\n";
}

close P;

system("gzip $PositionFilePath");
print "The script has finished with success\n";

#######################################################################################
#######################################################################################
######     		     			 	 	SUBROUTINES 	  		   		 	 	 ######
#######################################################################################
#######################################################################################

sub revcom {

    my($dna) = @_;
    my $revcom = reverse $dna;
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}
