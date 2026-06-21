#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
################																						 ################
################					 		FETCH GNOMAD VARIANTS GENOME						 		 ################
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
use List::Util qw(min max);
use IO::Zlib;
use File::Path qw(make_path);
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

#######################################################################################
#######################################################################################
######     		     			  Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$CHROM							= "";
my 	$FileHandle						= new IO::Zlib;

#######################################################################################
#######################################################################################
######     		     			  	    Get Options   	  		   			 	 ######
#######################################################################################
#######################################################################################

GetOptions ("CHROM=s"				=>	\$CHROM);

#######################################################################################
#######################################################################################
######     		     			  Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$VarFilePath					= "$cfg{paths}{gnomad_genome_source}/gnomad.genomes.v4.1.sites.chr$CHROM.vcf.bgz";
my 	$OutFilePath					= "$cfg{paths}{gnomad_genome_dir}/$CHROM.bed";
my 	$SortedOutFilePath				= "$cfg{paths}{gnomad_genome_dir}/$CHROM.sorted.bed";
my 	$TempDir						= "$cfg{paths}{data_root}/temp";

make_path($cfg{paths}{gnomad_genome_dir}) unless -d $cfg{paths}{gnomad_genome_dir};
make_path($TempDir)                       unless -d $TempDir;

#######################################################################################
#######################################################################################
######     		     			  		Iterate  	  		   			 		 ######
#######################################################################################
#######################################################################################

open O, ">$OutFilePath" or die ("Can't open $OutFilePath\n");

if ($FileHandle->open($VarFilePath, "rb")){
	
	while(<$FileHandle>){

		$_ =~ s/\n//g;
		
		next if ($_ =~ m/^#/);
		
		(my $Chromosome, 
		 my $Position, 
		 undef, 
		 my $ReferenceAllele, 
		 my $VariantAlleles, 
		 my $Qual, 
		 my $Filter, 
		 my $Info) 
		= split(m/\t/, $_);
		
		$Chromosome =~ s/^chr//g;
		
		my $Lcr = 0;
		
		if ($Filter eq "PASS"){
		
			my $AN			= "/";
			my $AC			= "/";
			my $AF			= "/";
			my @AF			= ();
			my @AC			= ();
			
			# Parse Allele Frequency #
			
			if ($Info =~ m/^.*;AF=(.+?);.*$/){
				$AF = $1;
				@AF = split (",", $AF);
			}
			else {die ("AF\t$_\n");}
			
			if ($Info =~ m/AC=(.+?);.*$/){
				$AC = $1;
				@AC = split (",", $AC);
			}
			else {die ("AC\t$_\n");}			
			
			if ($Info =~ m/^.*;AN=(.+?);.*$/){
				$AN = $1;
			}
			else {die ("AN\t$_\n");}
			
			# Lcr
			
			if ($Info =~ m/;lcr;/){$Lcr = 1;}
			
			my @VariantAlleles = split (",", $VariantAlleles);
			
			for (my $J = 0; $J < scalar @VariantAlleles; $J++){
			
				my $StartOfVariant		= 0;
				my $EndOfVariant		= 0;
				my $Variant				= "";
				
				###  FOR NOW WE ONLY TAKE SUBSTITUTIONS ###
				
				if (length ($ReferenceAllele) == 1 && length ($VariantAlleles[$J]) == 1){ 
					
					$EndOfVariant 	= $Position + length ($ReferenceAllele) - 1;						
					$Variant 		= "$Chromosome\_$Position\_$EndOfVariant\_$ReferenceAllele\/$VariantAlleles[$J]";
					
					my $Min 	= min($Position, $EndOfVariant); 
					my $Max 	= max($Position, $EndOfVariant);
					my $BedMin	= $Min-1;
					
					
					print O $Chromosome . "\t" . $BedMin . "\t" . $Max  . "\t" . $Variant . "\t" . "$AF\t$AC\t$AN\n";
				}
			}
		}
	}
}

$FileHandle->close();

close O;

system("TMPDIR=$TempDir && bedtools sort -i $OutFilePath > $SortedOutFilePath");
system("bgzip -c $SortedOutFilePath > $SortedOutFilePath.gz");
system("tabix -p bed $SortedOutFilePath.gz");

















