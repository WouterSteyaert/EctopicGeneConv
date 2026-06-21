#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
################																						 ################
################			  MERGE CHUNKED RECOMBHOTSPOTS OUTPUT FILES								 	 ################
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
######     		     			   Load Other Libraries  	  		   			 ######
#######################################################################################
#######################################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Glob ':glob';
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

#######################################################################################
#######################################################################################
######     		     			 Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my %cfg = load_config();

my $BASE         = $cfg{paths}{stats_dir};
my $Hg38DictFp   = $cfg{paths}{reference_dict};
my %ChromLengths = ();

#######################################################################################
#######################################################################################
######     		     	 Initialize And Declare - Query Parameters 	  			 ######
#######################################################################################
#######################################################################################

my @RepLengths      = split /,/, $cfg{repeat_lengths}{values};
my @QueryVarLists   = ("gnomad/genome");
my @QueryVarFreqs   = split /,/, $cfg{enrichment}{af_edges};
my $QueryWindowSize = $cfg{enrichment}{window_size} || 1_000_000;
my $QueryStepSize   = $cfg{enrichment}{step_size}   || 1_000_000;
my $RegionType      = "RecombHotspots";

my @StatsDirs       = map { "stats_$_" } split /,/, $cfg{enrichment}{contexts};

#######################################################################################
#######################################################################################
######     		     			    Read In Chrom Sizes  	  		   	 		 ######
#######################################################################################
#######################################################################################

print "Read In Chrom Sizes\n";

open D, "$Hg38DictFp" or die ("Can't open $Hg38DictFp\n");
while(<D>){

	$_ =~ s/\n//g;

	my 	@LineValues 	= split(m/\t/, $_);
	my 	$Chrom			= $LineValues[1];
		$Chrom			=~ s/SN\://g;
		$Chrom			=~ s/chr//g;

	if ($Chrom =~ m/^\d+$/ || $Chrom eq "X" || $Chrom eq "Y"){

		my 	$Length = $LineValues[2];
			$Length =~ s/LN\://g;

		$ChromLengths{$Chrom} = $Length;
	}
}
close D;

#######################################################################################
#######################################################################################
######     		     			  		  Iterate  	  		   			 		 ######
#######################################################################################
#######################################################################################

my $MergedCount = 0;
my $SkippedCount = 0;

foreach my $StatsDir (@StatsDirs){

	print "\n=== Processing $StatsDir ===\n";

	foreach my $RepLength (@RepLengths){

		my $InputDir = "$BASE/$StatsDir/$RepLength/";

		unless (-d $InputDir){
			print "Directory does not exist: $InputDir\n";
			next;
		}

		for (my $I = 1; $I <= 24; $I++){

			my 	$Chrom = $I;

			if 		($I == 23){$Chrom = "X";}
			elsif 	($I == 24){$Chrom = "Y";}

			foreach my $QueryVarList (@QueryVarLists){

				my 	$QueryVarListPrint 	= $QueryVarList;
					$QueryVarListPrint 	=~ s/\//~/g;

				for (my $FreqI = 1; $FreqI < scalar (@QueryVarFreqs); $FreqI++){

					my 	$QueryVarMinFreq 			= $QueryVarFreqs[$FreqI-1];
					my 	$QueryVarMaxFreq 			= $QueryVarFreqs[$FreqI];

					my 	$QueryVarMinFreqPrint		= $QueryVarMinFreq;
						$QueryVarMinFreqPrint		=~ s/\//~/;
						$QueryVarMinFreqPrint		=~ s/\.//;
					my 	$QueryVarMaxFreqPrint		= $QueryVarMaxFreq;
						$QueryVarMaxFreqPrint		=~ s/\//~/;
						$QueryVarMaxFreqPrint		=~ s/\.//;

					# Find all chunked files for this combination
					my $FilePattern = "$InputDir/$Chrom\_$QueryWindowSize\_$QueryStepSize\_$QueryVarListPrint\_$QueryVarMinFreqPrint\_$QueryVarMaxFreqPrint.*.$RegionType.txt";
					my @ChunkedFiles = bsd_glob($FilePattern);

					# Filter to only include chunked files (those with start.end pattern like 1.30000000)
					my @RealChunkedFiles = grep { /\.\d+\.\d+\.$RegionType\.txt$/ } @ChunkedFiles;

					if (scalar(@RealChunkedFiles) == 0){
						next;  # No chunked files found
					}

					# Sort files by chunk start position
					@RealChunkedFiles = sort {
						my ($a_start) = $a =~ /\.(\d+)\.\d+\.$RegionType\.txt$/;
						my ($b_start) = $b =~ /\.(\d+)\.\d+\.$RegionType\.txt$/;
						$a_start <=> $b_start;
					} @RealChunkedFiles;

					# Output file (merged, with full chromosome range)
					my $ChromLength = $ChromLengths{$Chrom};
					my $OutputFile = "$InputDir/$Chrom\_$QueryWindowSize\_$QueryStepSize\_$QueryVarListPrint\_$QueryVarMinFreqPrint\_$QueryVarMaxFreqPrint.1.$ChromLength.$RegionType.txt";

					# Check if merged file already exists and is newer than all chunks
					if (-e $OutputFile){
						my $OutputMtime = (stat($OutputFile))[9];
						my $AllOlder = 1;
						foreach my $ChunkFile (@RealChunkedFiles){
							my $ChunkMtime = (stat($ChunkFile))[9];
							if ($ChunkMtime > $OutputMtime){
								$AllOlder = 0;
								last;
							}
						}
						if ($AllOlder){
							$SkippedCount++;
							next;  # Skip, merged file is up to date
						}
					}

					print "Merging $Chrom RepLen=$RepLength Freq=$QueryVarMinFreq-$QueryVarMaxFreq (" . scalar(@RealChunkedFiles) . " chunks)\n";

					# Merge files
					open OUT, ">$OutputFile" or die ("Can't open $OutputFile for writing\n");

					my $HeaderWritten = 0;
					foreach my $ChunkFile (@RealChunkedFiles){

						open IN, "<$ChunkFile" or die ("Can't open $ChunkFile\n");
						while(<IN>){
							# Skip header lines except for first file
							if (/^QueryChrom\t/){
								if (!$HeaderWritten){
									print OUT $_;
									$HeaderWritten = 1;
								}
								next;
							}
							print OUT $_;
						}
						close IN;
					}

					close OUT;
					$MergedCount++;

					# Optionally delete chunk files after successful merge
					# Uncomment the following lines to enable automatic cleanup:
					# foreach my $ChunkFile (@RealChunkedFiles){
					# 	unlink $ChunkFile;
					# }
				}
			}
		}
	}
}

print "\n=== Summary ===\n";
print "Merged: $MergedCount combinations\n";
print "Skipped (up to date): $SkippedCount combinations\n";
