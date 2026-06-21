#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
################																						 ################
################					 		 COMPUTE TRACT LENGTH GNOMAD GENOME						 	 ################
################																						 ################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

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
use FindBin;
use lib "$FindBin::Bin/../00_Configuration";
use ProjectConfig qw(load_config);

my %cfg = load_config();

#######################################################################################
#######################################################################################
######     		     			Initialize And Declare  	  		   			 ######
#######################################################################################
#######################################################################################

my 	$VARIANT				= "11_70560_70560_C/T";
my 	$FREQ					= "/";
my 	@TractLengthRanges		= split /,/, $cfg{tract}{window_sizes};
my 	@RepeatLengths			= split /,/, $cfg{repeat_lengths}{values};
my 	$Hg38_2bitFp			= $cfg{paths}{reference_2bit};
my 	$RefSeqHg38				= $cfg{paths}{reference_fasta};
my 	$TractLengthDir			= "$cfg{paths}{tractlengths_wgs_dir}/";
my 	$LineNr 				= 0;
my 	%AllTractLengths		= ();
my 	%AllRepeatLengths		= ();

#######################################################################################
#######################################################################################
######     		     			  	  Get Options   	  		   			 	 ######
#######################################################################################
#######################################################################################

GetOptions ("VARIANT=s"		=>	\$VARIANT,
			"FREQ=s"		=>	\$FREQ);

#######################################################################################
#######################################################################################
######     		     		    Iterate Over TractLengths  		   			 	 ######
#######################################################################################
#######################################################################################

(my $Chrom, my $VarPos, undef, my $Change) 	= split(m/_/, $VARIANT);
(my $Ref, my $Alt) 							= split(m/\//, $Change);

my 	$ResultFilePath			= $TractLengthDir . $Chrom . "_" . $VarPos . ".result.txt";

foreach my $TractLengthRange (@TractLengthRanges){
	
	my $AcceptorStart			= $VarPos - $TractLengthRange;
	my $AcceptorEnd				= $VarPos + $TractLengthRange;
	my $AcceptorTwoBitStart		= $AcceptorStart-1;
	my $AcceptorTwoBitEnd		= $AcceptorEnd+1;
	
	### MAKE REFSEQ SLICE ###
	
	my @RefSeqSlice 			= `twoBitToFa $Hg38_2bitFp:chr$Chrom:$AcceptorTwoBitStart-$AcceptorTwoBitEnd stdout`;
	my $RefSeqSlice				= "";
	
	foreach my $Line (@RefSeqSlice){
		
		if ($Line !~ m/^>/){
		
			$Line =~ s/\n//g;
			
			$RefSeqSlice .= $Line;
		}
	}
	
	### MAKE MUTATED SLICE ###
	
	my $MutatedSlice 			= IntroduceMutationInSlice($RefSeqSlice, $VARIANT, $TractLengthRange);
	
	### BLAST MUTATED SLICE ###
	
	unless (-d "$TractLengthDir\/$Chrom\_$VarPos\/"){
		system("mkdir -p $TractLengthDir\/$Chrom\_$VarPos\/");
	}
						
	my $QueryIn 				= "$TractLengthDir\/$Chrom\_$VarPos\/$TractLengthRange.fa";
	my $QueryOut 				= "$TractLengthDir\/$Chrom\_$VarPos\/$TractLengthRange.out";
	my $OutputFileP				= "$TractLengthDir\/$Chrom\_$VarPos\/$TractLengthRange.txt";

	open O, ">$QueryIn" or die ("Can't open $QueryIn");
	print O ">$Chrom\_$VarPos\n";
	print O "$MutatedSlice\n";
	close O;
	
	print "blastn -db $RefSeqHg38 -query \"$QueryIn\" -out \"$QueryOut\" -perc_identity $cfg{tract}{blast_perc_identity} -dust no -word_size $cfg{tract}{blast_word_size}\n";
	system("blastn -db $RefSeqHg38 -query \"$QueryIn\" -out \"$QueryOut\" -perc_identity $cfg{tract}{blast_perc_identity} -dust no -word_size $cfg{tract}{blast_word_size}");
	
	### COMPUTE TRACT LENGTH ###
	
	my 	$NewQueryBlock			= 0;
	my 	$NewSubjectBlock		= 0;
	my 	$Strand					= "";
	my 	$Subject				= "";
	my 	$Query					= "";
	my 	$ChromName				= "";
	my 	$SubjectCoord			= 0;
	my 	$QueryCoord				= 0;
	my 	$FirstSubjectCoord		= 0;
	my 	$LastSubjectCoord		= 0;
	my 	$FirstQueryCoord		= 0;
	my 	$LastQueryCoord			= 0;
	
	open O, ">$OutputFileP" or die ("Can't open $OutputFileP\n");
	
	print O "QUERY_VARIANT\tFREQ\tChrom\tAcceptorStart\tAcceptorEnd\tChromName\tFirstSubjectCoord\tLastSubjectCoord\tNrOfMatches\tSUNs\tHomSUNs\tTractLengthRange\tStrand\tLargestSmallerHom\tSmallestLargerHom\tLargestSmallerSunHom\tSmallestLargerSunHom\tLargestSmallerSun\tSmallestLargerSun\tLargestSmaller\tSmallestLarger\tTractLength\tTractAnalysisLength\tLeftTractFlank\tRightTractFlank\tFlankRatio\n";
	
	open B, "$QueryOut" or die ("Can't open $QueryOut\n"); 
	while(<B>){
		
		$_ =~ s/\n//g;
		
		if ($_ =~ m/^>/){
			
			# write out previous block if present #
			# Subject = potential donor
			
			if ($Subject ne "" && $ChromName !~ m/_/){
				
				my $FirstMatch 						= 0;
				my @SubjectBases 					= split(m//, $Subject);
				my @QueryBases 						= split(m//, $Query);
				my $SubstitutionRef 				= "";
				my $SubstitutionAlt					= "";
				my $DeletionRef						= "";
				my $InsertionAlt					= "";
				my $NoMatch							= 0;
				my $PosInQuery						= $AcceptorStart-1 + $QueryCoord;
				my $NrOfSubjectBases				= 0;
				my %HomoloVarMappings				= ();
				my %RevHomoloVarMappings			= ();
				my %SUNs							= ();
				my $AcgtLen 						= ($Query =~ tr/ACGTacgt//);

				if ($AcgtLen>$TractLengthRange){
				
					# Check if central base matches with subject base #
					
					my $I 	= $QueryCoord;
					my $J 	= 0;
					
					foreach my $QueryBase (@QueryBases){

						if ($I >=$TractLengthRange && $I < $TractLengthRange+length($Alt)){
							
							if ($SubjectBases[$J] ne $QueryBases[$J]){
								$NoMatch = 1;
								last;
							}
						}
						
						if ($QueryBase =~ m/A|C|G|T/i){$I++;}
						$J++;
					}
					
					# Compute Number Of Bases In Subject
					
					foreach my $SubjectBase (@SubjectBases){
						
						if ($SubjectBase =~ m/A|C|G|T/i){$NrOfSubjectBases++;}
					}
					
					if ($Strand !~ m/Minus/){
						
						$LastSubjectCoord = $FirstSubjectCoord+$NrOfSubjectBases;
					}
					else {
						
						$LastSubjectCoord = $FirstSubjectCoord-$NrOfSubjectBases;
					}
					
					if ($NoMatch == 0){ # Check other bases # Central base matches
					
						my $NrOfMatches = 0;

						for (my $I = 0; $I < scalar @QueryBases; $I++){ # acceptor (mutated reference sequence)
							
							if ($QueryBases[$I] eq $SubjectBases[$I])			{$NrOfMatches++;}
							if ($QueryBases[$I] =~ m/A|C|G|T/i)					{$PosInQuery++;}
							
							# we dont trust sequence differences at the edges of the query
							if ($QueryBases[$I] =~ m/A|C|G|T/i 	&& 
								$SubjectBases[$I] eq $QueryBases[$I])			{$FirstMatch = 1;} 
							
							if ($FirstSubjectCoord < $LastSubjectCoord 	&& 
								$SubjectBases[$I] =~ m/A|C|G|T/i){ # Plus strand
								
								$SubjectCoord++;
							}
							elsif ($SubjectBases[$I] =~ m/A|C|G|T/i){
								
								$SubjectCoord--;
							}
						
							if ($FirstMatch == 1){ # we only execute the code below from the first match onwards
								
								# SUB
							
								if (	$SubjectBases[$I] 	=~ m/A|C|G|T/i 	&& # pot. donor
										$QueryBases[$I] 	=~ m/A|C|G|T/i 	&& # pot. acceptor
										$SubjectBases[$I] 	ne $QueryBases[$I]){
									
									$SubstitutionRef	.= $QueryBases[$I];
									$SubstitutionAlt	.= $SubjectBases[$I];
								}
								elsif ($SubstitutionRef ne ""){ # write out when difference ends
									
									my $VariantEnd 		= $PosInQuery-1;
									my $VariantStart	= $VariantEnd-length($SubstitutionRef)+1;
									my $EnsemblVariant	= $Chrom . "_" . $VariantStart . "_" . $VariantEnd . "_" . $SubstitutionRef . "/" . $SubstitutionAlt;
									
									$SUNs{$EnsemblVariant} = undef;
									
									# find homolo variant 
									
									if ($Strand !~ m/Minus/){ # plus
										
										my $HomVariantEnd 		= $SubjectCoord-1;
										my $HomVariantStart		= $HomVariantEnd-length($SubstitutionRef)+1;
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $SubstitutionAlt . "/" . $SubstitutionRef;
										
										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									else {
										
										my $HomVariantStart 	= $SubjectCoord+1;
										my $HomVariantEnd		= $HomVariantStart+length($SubstitutionRef)-1;
										my $HomSubAlt			= scalar reverse($SubstitutionRef);
										my $HomSubRef			= scalar reverse($SubstitutionAlt);
										
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $HomSubRef . "/" . $HomSubAlt;
										
										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									
									$SubstitutionRef	= "";
									$SubstitutionAlt	= "";
								}
							
								# DEL
								
								if ($SubjectBases[$I] eq "-"){
									
									$DeletionRef		.= $QueryBases[$I];
								}
								elsif ($DeletionRef ne ""){
									
									my $VariantEnd 			= $PosInQuery-1;
									my $VariantStart		= $VariantEnd-length($DeletionRef)+1;
									my $EnsemblVariant		= $Chrom . "_" . $VariantStart . "_" . $VariantEnd . "_" . $DeletionRef . "/" . "-";
									
									$SUNs{$EnsemblVariant} = undef;
									
									# it becomes an insertion
									
									if ($Strand !~ m/Minus/){ # plus

										my $HomVariantStart 	= $SubjectCoord;
										my $HomVariantEnd		= $SubjectCoord-1;
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . "-" . "/" . $DeletionRef;
										
										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									else {
										
										my $HomVariantStart 	= $SubjectCoord+1;
										my $HomVariantEnd		= $SubjectCoord;
										my $HomIns				= scalar reverse($DeletionRef);
										
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . "-" . "/" . $HomIns;

										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}

									$DeletionRef 		= "";
								}
								
								# INS
								
								if 		($QueryBases[$I] eq "-"){
									
									$InsertionAlt		.= $SubjectBases[$I];
								}
								elsif 	($InsertionAlt ne ""){
									
									my $VariantStart 		= $PosInQuery;
									my $VariantEnd			= $VariantStart-1;
									my $EnsemblVariant		= $Chrom . "_" . $VariantStart . "_" . $VariantEnd . "_" . "-" . "/" . $InsertionAlt;
									
									$SUNs{$EnsemblVariant} 	= undef;
									
									# it becomes a deletion
									
									if ($Strand !~ m/Minus/){ # plus

										my $HomVariantEnd 		= $SubjectCoord-1;
										my $HomVariantStart		= $HomVariantEnd-length($SubstitutionRef)+1;
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $InsertionAlt . "/" . "-";
										
										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									else {
										
										my $HomVariantStart 	= $SubjectCoord+1;
										my $HomVariantEnd		= $HomVariantStart+length($SubstitutionRef)-1;
										my $HomDel				= scalar reverse($InsertionAlt);
										
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $HomDel . "/" . "-";

										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									
									$InsertionAlt 		= "";
								}
							}	
						}
						
						my $SUNs 						= join(', ', sort keys %SUNs);
						my @SUNs						= keys %SUNs;
						my $HomSUNs						= "";
						
						foreach my $Sun (sort keys %SUNs){
							
							$HomSUNs 				   .= $HomoloVarMappings{$Sun} . ", ";
						}

						$HomSUNs						=~ s/, $//g;
						
						# calculate donor interval #

						# smallest non introduced sun > DNM
						# largest non introduced sun < DNM
						# this is always the case since this is based on acceptor coordinates
						
						my $LargestSmaller			= $AcceptorStart-1 + $QueryCoord;
						my $SmallestLarger			= $PosInQuery + 1;
						my $SmallestLargerSun		= "";
						my $LargestSmallerSun		= "";
						my $LargestSmallerHom		= $FirstSubjectCoord; # beginning of interval
						my $SmallestLargerHom		= $LastSubjectCoord;  # end of interval
						
						foreach my $Sun (@SUNs){
							
							(undef, my $SunPosStart, my $SunPosEnd, undef) = split(m/_/, $Sun);
																					
							if 	($SunPosEnd < $VarPos){
								
								if ($SunPosEnd > $LargestSmaller){
									
									$LargestSmaller 	= $SunPosEnd;
									$LargestSmallerSun	= $Sun;
								}
							}
							elsif ($SunPosStart > $VarPos){
								
								if ($SunPosStart < $SmallestLarger){
									
									$SmallestLarger 	= $SunPosStart;
									$SmallestLargerSun	= $Sun;
								}
							}
						}
						
						my $SmallestLargerSunHom = "/";
						my $LargestSmallerSunHom = "/";
						my $TractLength			 = $SmallestLarger-$LargestSmaller-1;
						
						if ($SmallestLargerSun ne ""){
							
							$SmallestLargerSunHom = $HomoloVarMappings{$SmallestLargerSun};
							
							(undef, $SmallestLargerHom, undef) = split(m/_/, $SmallestLargerSunHom);
							
						}
						
						if ($LargestSmallerSun ne ""){
							
							$LargestSmallerSunHom = $HomoloVarMappings{$LargestSmallerSun};
							
							(undef, $LargestSmallerHom, undef) = split(m/_/, $LargestSmallerSunHom);
						}
						
						my $MinFinalDonor		= min($SmallestLargerHom,$LargestSmallerHom);
						my $MaxFinalDonor		= max($SmallestLargerHom,$LargestSmallerHom);
						
						# Compute TractAnalysisLength #
						
						my $TractAnalysisLength		= 0;
						my $LeftDiff 				= $VarPos - $LargestSmaller;
						my $RightDiff 				= $SmallestLarger - $VarPos;
						
						my $LeftTractFlank			= $LeftDiff-1;
						my $RightTractFlank			= $RightDiff-1;
						
						
						my $Min 					= min($LeftDiff, $RightDiff);
						my $FlankRatio				= $Min/($LeftTractFlank+1+$RightTractFlank);
						
						
						foreach my $RepeatLength (@RepeatLengths){ # find largest which is still smaller than repeatlength
						
							my $Flank = ($RepeatLength-1)/2;
							
							if($Flank < $Min){
								$TractAnalysisLength = $RepeatLength;
							}
						}
						
						# print #
							
						print O $VARIANT . "\t" . $FREQ . "\t" . $Chrom . "\t" . $AcceptorStart . "\t" . $AcceptorEnd  . "\t" .  $ChromName  . "\t" . $FirstSubjectCoord  . "\t" . $LastSubjectCoord . "\t" . $NrOfMatches . "\t" . $SUNs . "\t" . $HomSUNs . "\t" . $TractLengthRange. "\t" . $Strand . "\t" . $LargestSmallerHom . "\t" . $SmallestLargerHom . "\t" . $LargestSmallerSunHom . "\t" . $SmallestLargerSunHom . "\t" . $LargestSmallerSun . "\t" . $SmallestLargerSun . "\t" . $LargestSmaller . "\t" . $SmallestLarger . "\t" . $TractLength . "\t" . $TractAnalysisLength . "\t" . $LeftTractFlank . "\t" . $RightTractFlank . "\t" . $FlankRatio . "\n";
						
					}
				}
			}
		
			# set for new block
			
			$NewQueryBlock			= 1;
			$NewSubjectBlock		= 1;
			$Subject				= "";
			$Query					= "";
			
			# new chromosome
			
			$ChromName 	= $_;
			$ChromName	=~ s/^> //g;
			$ChromName	=~ s/chr//g;
		}
		
		if ($_ =~ m/Score/){
			
			# check previous block
			
			if ($Subject ne "" && $ChromName !~ m/_/){
				
				my $FirstMatch 						= 0;
				my @SubjectBases 					= split(m//, $Subject);
				my @QueryBases 						= split(m//, $Query);
				my $SubstitutionRef 				= "";
				my $SubstitutionAlt					= "";
				my $DeletionRef						= "";
				my $InsertionAlt					= "";
				my $NoMatch							= 0;
				my $PosInQuery						= $AcceptorStart-1 + $QueryCoord;
				my $NrOfSubjectBases				= 0;
				my %HomoloVarMappings				= ();
				my %RevHomoloVarMappings			= ();
				my %SUNs							= ();
				my $AcgtLen 						= ($Query =~ tr/ACGTacgt//);

				if ($AcgtLen>$TractLengthRange){
				
					# Check if central base matches with subject base #
					
					my $I 	= $QueryCoord;
					my $J 	= 0;
					
					foreach my $QueryBase (@QueryBases){

						if ($I >=$TractLengthRange && $I < $TractLengthRange+length($Alt)){
							
							if ($SubjectBases[$J] ne $QueryBases[$J]){
								$NoMatch = 1;
								last;
							}
						}
						
						if ($QueryBase =~ m/A|C|G|T/i){$I++;}
						$J++;
					}
					
					# Compute Number Of Bases In Subject
					
					foreach my $SubjectBase (@SubjectBases){
						
						if ($SubjectBase =~ m/A|C|G|T/i){$NrOfSubjectBases++;}
					}
					
					if ($Strand !~ m/Minus/){
						
						$LastSubjectCoord = $FirstSubjectCoord+$NrOfSubjectBases;
					}
					else {
						
						$LastSubjectCoord = $FirstSubjectCoord-$NrOfSubjectBases;
					}
					
					if ($NoMatch == 0){ # Check other bases # Central base matches
					
						my $NrOfMatches = 0;

						for (my $I = 0; $I < scalar @QueryBases; $I++){ # acceptor (mutated reference sequence)
							
							if ($QueryBases[$I] eq $SubjectBases[$I])			{$NrOfMatches++;}
							if ($QueryBases[$I] =~ m/A|C|G|T/i)					{$PosInQuery++;}
							
							# we dont trust sequence differences at the edges of the query
							if ($QueryBases[$I] =~ m/A|C|G|T/i 	&& 
								$SubjectBases[$I] eq $QueryBases[$I])			{$FirstMatch = 1;} 
							
							if ($FirstSubjectCoord < $LastSubjectCoord 	&& 
								$SubjectBases[$I] =~ m/A|C|G|T/i){ # Plus strand
								
								$SubjectCoord++;
							}
							elsif ($SubjectBases[$I] =~ m/A|C|G|T/i){
								
								$SubjectCoord--;
							}
						
							if ($FirstMatch == 1){ # we only execute the code below from the first match onwards
								
								# SUB
							
								if (	$SubjectBases[$I] 	=~ m/A|C|G|T/i 	&& # pot. donor
										$QueryBases[$I] 	=~ m/A|C|G|T/i 	&& # pot. acceptor
										$SubjectBases[$I] 	ne $QueryBases[$I]){
									
									$SubstitutionRef	.= $QueryBases[$I];
									$SubstitutionAlt	.= $SubjectBases[$I];
								}
								elsif ($SubstitutionRef ne ""){ # write out when difference ends
									
									my $VariantEnd 		= $PosInQuery-1;
									my $VariantStart	= $VariantEnd-length($SubstitutionRef)+1;
									my $EnsemblVariant	= $Chrom . "_" . $VariantStart . "_" . $VariantEnd . "_" . $SubstitutionRef . "/" . $SubstitutionAlt;
									
									$SUNs{$EnsemblVariant} = undef;
									
									# find homolo variant 
									
									if ($Strand !~ m/Minus/){ # plus
										
										my $HomVariantEnd 		= $SubjectCoord-1;
										my $HomVariantStart		= $HomVariantEnd-length($SubstitutionRef)+1;
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $SubstitutionAlt . "/" . $SubstitutionRef;
										
										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									else {
										
										my $HomVariantStart 	= $SubjectCoord+1;
										my $HomVariantEnd		= $HomVariantStart+length($SubstitutionRef)-1;
										my $HomSubAlt			= scalar reverse($SubstitutionRef);
										my $HomSubRef			= scalar reverse($SubstitutionAlt);
										
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $HomSubRef . "/" . $HomSubAlt;
										
										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									
									$SubstitutionRef	= "";
									$SubstitutionAlt	= "";
								}
							
								# DEL
								
								if ($SubjectBases[$I] eq "-"){
									
									$DeletionRef		.= $QueryBases[$I];
								}
								elsif ($DeletionRef ne ""){
									
									my $VariantEnd 			= $PosInQuery-1;
									my $VariantStart		= $VariantEnd-length($DeletionRef)+1;
									my $EnsemblVariant		= $Chrom . "_" . $VariantStart . "_" . $VariantEnd . "_" . $DeletionRef . "/" . "-";
									
									$SUNs{$EnsemblVariant} = undef;
									
									# it becomes an insertion
									
									if ($Strand !~ m/Minus/){ # plus

										my $HomVariantStart 	= $SubjectCoord;
										my $HomVariantEnd		= $SubjectCoord-1;
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . "-" . "/" . $DeletionRef;
										
										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									else {
										
										my $HomVariantStart 	= $SubjectCoord+1;
										my $HomVariantEnd		= $SubjectCoord;
										my $HomIns				= scalar reverse($DeletionRef);
										
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . "-" . "/" . $HomIns;

										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}

									$DeletionRef 		= "";
								}
								
								# INS
								
								if 		($QueryBases[$I] eq "-"){
									
									$InsertionAlt		.= $SubjectBases[$I];
								}
								elsif 	($InsertionAlt ne ""){
									
									my $VariantStart 		= $PosInQuery;
									my $VariantEnd			= $VariantStart-1;
									my $EnsemblVariant		= $Chrom . "_" . $VariantStart . "_" . $VariantEnd . "_" . "-" . "/" . $InsertionAlt;
									
									$SUNs{$EnsemblVariant} 	= undef;
									
									# it becomes a deletion
									
									if ($Strand !~ m/Minus/){ # plus

										my $HomVariantEnd 		= $SubjectCoord-1;
										my $HomVariantStart		= $HomVariantEnd-length($SubstitutionRef)+1;
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $InsertionAlt . "/" . "-";
										
										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									else {
										
										my $HomVariantStart 	= $SubjectCoord+1;
										my $HomVariantEnd		= $HomVariantStart+length($SubstitutionRef)-1;
										my $HomDel				= scalar reverse($InsertionAlt);
										
										my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $HomDel . "/" . "-";

										$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
										$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
									}
									
									$InsertionAlt 		= "";
								}
							}	
						}
						
						my $SUNs 						= join(', ', sort keys %SUNs);
						my @SUNs						= keys %SUNs;
						my $HomSUNs						= "";
						
						foreach my $Sun (sort keys %SUNs){
							
							$HomSUNs 				   .= $HomoloVarMappings{$Sun} . ", ";
						}

						$HomSUNs						=~ s/, $//g;
						
						# calculate donor interval #

						# smallest non introduced sun > DNM
						# largest non introduced sun < DNM
						# this is always the case since this is based on acceptor coordinates
						
						my $LargestSmaller			= $AcceptorStart-1 + $QueryCoord;
						my $SmallestLarger			= $PosInQuery + 1;
						my $SmallestLargerSun		= "";
						my $LargestSmallerSun		= "";
						my $LargestSmallerHom		= $FirstSubjectCoord; # beginning of interval
						my $SmallestLargerHom		= $LastSubjectCoord;  # end of interval
						
						foreach my $Sun (@SUNs){
							
							(undef, my $SunPosStart, my $SunPosEnd, undef) = split(m/_/, $Sun);
																					
							if 	($SunPosEnd < $VarPos){
								
								if ($SunPosEnd > $LargestSmaller){
									
									$LargestSmaller 	= $SunPosEnd;
									$LargestSmallerSun	= $Sun;
								}
							}
							elsif ($SunPosStart > $VarPos){
								
								if ($SunPosStart < $SmallestLarger){
									
									$SmallestLarger 	= $SunPosStart;
									$SmallestLargerSun	= $Sun;
								}
							}
						}
						
						my $SmallestLargerSunHom = "/";
						my $LargestSmallerSunHom = "/";
						my $TractLength			 = $SmallestLarger-$LargestSmaller-1;
						
						if ($SmallestLargerSun ne ""){
							
							$SmallestLargerSunHom = $HomoloVarMappings{$SmallestLargerSun};
							
							(undef, $SmallestLargerHom, undef) = split(m/_/, $SmallestLargerSunHom);
							
						}
						
						if ($LargestSmallerSun ne ""){
							
							$LargestSmallerSunHom = $HomoloVarMappings{$LargestSmallerSun};
							
							(undef, $LargestSmallerHom, undef) = split(m/_/, $LargestSmallerSunHom);
						}
						
						my $MinFinalDonor		= min($SmallestLargerHom,$LargestSmallerHom);
						my $MaxFinalDonor		= max($SmallestLargerHom,$LargestSmallerHom);
						
						# Compute TractAnalysisLength #
						
						my $TractAnalysisLength		= 0;
						my $LeftDiff 				= $VarPos - $LargestSmaller;
						my $RightDiff 				= $SmallestLarger - $VarPos;
						
						my $LeftTractFlank			= $LeftDiff-1;
						my $RightTractFlank			= $RightDiff-1;
						
						
						my $Min 					= min($LeftDiff, $RightDiff);
						my $FlankRatio				= $Min/($LeftTractFlank+1+$RightTractFlank);
						
						foreach my $RepeatLength (@RepeatLengths){ # find largest which is still smaller than repeatlength
						
							my $Flank = ($RepeatLength-1)/2;
							
							if($Flank < $Min){
								$TractAnalysisLength = $RepeatLength;
							}
						}
						
						# print #
							
						print O $VARIANT . "\t" . $FREQ . "\t" . $Chrom . "\t" . $AcceptorStart . "\t" . $AcceptorEnd  . "\t" .  $ChromName  . "\t" . $FirstSubjectCoord  . "\t" . $LastSubjectCoord . "\t" . $NrOfMatches . "\t" . $SUNs . "\t" . $HomSUNs . "\t" . $TractLengthRange. "\t" . $Strand . "\t" . $LargestSmallerHom . "\t" . $SmallestLargerHom . "\t" . $LargestSmallerSunHom . "\t" . $SmallestLargerSunHom . "\t" . $LargestSmallerSun . "\t" . $SmallestLargerSun . "\t" . $LargestSmaller . "\t" . $SmallestLarger . "\t" . $TractLength . "\t" . $TractAnalysisLength . "\t" . $LeftTractFlank . "\t" . $RightTractFlank . "\t" . $FlankRatio . "\n";
						
					}
				}
			}
			
			# set for new block
			
			$NewQueryBlock			= 1;
			$NewSubjectBlock		= 1;
			$Subject				= "";
			$Query					= "";
		}
		
		if ($_ =~ m/Strand/){
			
			$Strand = $_;
			$Strand =~ s/ Strand\=//g;
		}
		
		if ($_ =~ m/^Sbjct/){
			
			$_ =~ s/\s+/ /g;
								
			my @LineValues = split(m/ /, $_);
			$Subject .= uc($LineValues[2]);
			
			if ($NewSubjectBlock == 1){
				
				$FirstSubjectCoord = $LineValues[1];
				$NewSubjectBlock = 0;
			}
			
			if ($Strand =~ m/Minus/){
				$SubjectCoord = $FirstSubjectCoord+1;
			}
			else {
				$SubjectCoord = $FirstSubjectCoord-1;
			}
		}

		if ($_ =~ m/^Query/ && $_ !~ m/^Query\=/){
			
			$_ =~ s/\s+/ /g;
								
			my @LineValues = split(m/ /, $_);
			$Query .= uc($LineValues[2]);
			
			if ($NewQueryBlock == 1){
				
				$FirstQueryCoord = $LineValues[1];
				$NewQueryBlock = 0;
			}
			
			$QueryCoord = $FirstQueryCoord-1;
		}
	}
	close B;
	
	if ($Subject ne "" && $ChromName !~ m/_/){
				
		my $FirstMatch 						= 0;
		my @SubjectBases 					= split(m//, $Subject);
		my @QueryBases 						= split(m//, $Query);
		my $SubstitutionRef 				= "";
		my $SubstitutionAlt					= "";
		my $DeletionRef						= "";
		my $InsertionAlt					= "";
		my $NoMatch							= 0;
		my $PosInQuery						= $AcceptorStart-1 + $QueryCoord;
		my $NrOfSubjectBases				= 0;
		my %HomoloVarMappings				= ();
		my %RevHomoloVarMappings			= ();
		my %SUNs							= ();
		my $AcgtLen 						= ($Query =~ tr/ACGTacgt//);

		if ($AcgtLen>$TractLengthRange){
		
			# Check if central base matches with subject base #
			
			my $I 	= $QueryCoord;
			my $J 	= 0;
			
			foreach my $QueryBase (@QueryBases){

				if ($I >=$TractLengthRange && $I < $TractLengthRange+length($Alt)){
					
					if ($SubjectBases[$J] ne $QueryBases[$J]){
						$NoMatch = 1;
						last;
					}
				}
				
				if ($QueryBase =~ m/A|C|G|T/i){$I++;}
				$J++;
			}
			
			# Compute Number Of Bases In Subject
			
			foreach my $SubjectBase (@SubjectBases){
				
				if ($SubjectBase =~ m/A|C|G|T/i){$NrOfSubjectBases++;}
			}
			
			if ($Strand !~ m/Minus/){
				
				$LastSubjectCoord = $FirstSubjectCoord+$NrOfSubjectBases;
			}
			else {
				
				$LastSubjectCoord = $FirstSubjectCoord-$NrOfSubjectBases;
			}
			
			if ($NoMatch == 0){ # Check other bases # Central base matches
			
				my $NrOfMatches = 0;

				for (my $I = 0; $I < scalar @QueryBases; $I++){ # acceptor (mutated reference sequence)
					
					if ($QueryBases[$I] eq $SubjectBases[$I])			{$NrOfMatches++;}
					if ($QueryBases[$I] =~ m/A|C|G|T/i)					{$PosInQuery++;}
					
					# we dont trust sequence differences at the edges of the query
					if ($QueryBases[$I] =~ m/A|C|G|T/i 	&& 
						$SubjectBases[$I] eq $QueryBases[$I])			{$FirstMatch = 1;} 
					
					if ($FirstSubjectCoord < $LastSubjectCoord 	&& 
						$SubjectBases[$I] =~ m/A|C|G|T/i){ # Plus strand
						
						$SubjectCoord++;
					}
					elsif ($SubjectBases[$I] =~ m/A|C|G|T/i){
						
						$SubjectCoord--;
					}
				
					if ($FirstMatch == 1){ # we only execute the code below from the first match onwards
						
						# SUB
					
						if (	$SubjectBases[$I] 	=~ m/A|C|G|T/i 	&& # pot. donor
								$QueryBases[$I] 	=~ m/A|C|G|T/i 	&& # pot. acceptor
								$SubjectBases[$I] 	ne $QueryBases[$I]){
							
							$SubstitutionRef	.= $QueryBases[$I];
							$SubstitutionAlt	.= $SubjectBases[$I];
						}
						elsif ($SubstitutionRef ne ""){ # write out when difference ends
							
							my $VariantEnd 		= $PosInQuery-1;
							my $VariantStart	= $VariantEnd-length($SubstitutionRef)+1;
							my $EnsemblVariant	= $Chrom . "_" . $VariantStart . "_" . $VariantEnd . "_" . $SubstitutionRef . "/" . $SubstitutionAlt;
							
							$SUNs{$EnsemblVariant} = undef;
							
							# find homolo variant 
							
							if ($Strand !~ m/Minus/){ # plus
								
								my $HomVariantEnd 		= $SubjectCoord-1;
								my $HomVariantStart		= $HomVariantEnd-length($SubstitutionRef)+1;
								my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $SubstitutionAlt . "/" . $SubstitutionRef;
								
								$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
								$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
							}
							else {
								
								my $HomVariantStart 	= $SubjectCoord+1;
								my $HomVariantEnd		= $HomVariantStart+length($SubstitutionRef)-1;
								my $HomSubAlt			= scalar reverse($SubstitutionRef);
								my $HomSubRef			= scalar reverse($SubstitutionAlt);
								
								my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $HomSubRef . "/" . $HomSubAlt;
								
								$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
								$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
							}
							
							$SubstitutionRef	= "";
							$SubstitutionAlt	= "";
						}
					
						# DEL
						
						if ($SubjectBases[$I] eq "-"){
							
							$DeletionRef		.= $QueryBases[$I];
						}
						elsif ($DeletionRef ne ""){
							
							my $VariantEnd 			= $PosInQuery-1;
							my $VariantStart		= $VariantEnd-length($DeletionRef)+1;
							my $EnsemblVariant		= $Chrom . "_" . $VariantStart . "_" . $VariantEnd . "_" . $DeletionRef . "/" . "-";
							
							$SUNs{$EnsemblVariant} = undef;
							
							# it becomes an insertion
							
							if ($Strand !~ m/Minus/){ # plus

								my $HomVariantStart 	= $SubjectCoord;
								my $HomVariantEnd		= $SubjectCoord-1;
								my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . "-" . "/" . $DeletionRef;
								
								$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
								$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
							}
							else {
								
								my $HomVariantStart 	= $SubjectCoord+1;
								my $HomVariantEnd		= $SubjectCoord;
								my $HomIns				= scalar reverse($DeletionRef);
								
								my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . "-" . "/" . $HomIns;

								$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
								$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
							}

							$DeletionRef 		= "";
						}
						
						# INS
						
						if 		($QueryBases[$I] eq "-"){
							
							$InsertionAlt		.= $SubjectBases[$I];
						}
						elsif 	($InsertionAlt ne ""){
							
							my $VariantStart 		= $PosInQuery;
							my $VariantEnd			= $VariantStart-1;
							my $EnsemblVariant		= $Chrom . "_" . $VariantStart . "_" . $VariantEnd . "_" . "-" . "/" . $InsertionAlt;
							
							$SUNs{$EnsemblVariant} 	= undef;
							
							# it becomes a deletion
							
							if ($Strand !~ m/Minus/){ # plus

								my $HomVariantEnd 		= $SubjectCoord-1;
								my $HomVariantStart		= $HomVariantEnd-length($SubstitutionRef)+1;
								my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $InsertionAlt . "/" . "-";
								
								$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
								$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
							}
							else {
								
								my $HomVariantStart 	= $SubjectCoord+1;
								my $HomVariantEnd		= $HomVariantStart+length($SubstitutionRef)-1;
								my $HomDel				= scalar reverse($InsertionAlt);
								
								my $HomoloVariant		= $ChromName . "_" . $HomVariantStart . "_" . $HomVariantEnd . "_" . $HomDel . "/" . "-";

								$HomoloVarMappings{$EnsemblVariant} = $HomoloVariant;
								$RevHomoloVarMappings{$HomoloVariant} = $EnsemblVariant;
							}
							
							$InsertionAlt 		= "";
						}
					}	
				}
				
				my $SUNs 						= join(', ', sort keys %SUNs);
				my @SUNs						= keys %SUNs;
				my $HomSUNs						= "";
				
				foreach my $Sun (sort keys %SUNs){
					
					$HomSUNs 				   .= $HomoloVarMappings{$Sun} . ", ";
				}

				$HomSUNs						=~ s/, $//g;
				
				# calculate donor interval #

				# smallest non introduced sun > DNM
				# largest non introduced sun < DNM
				# this is always the case since this is based on acceptor coordinates
				
				my $LargestSmaller			= $AcceptorStart-1 + $QueryCoord;
				my $SmallestLarger			= $PosInQuery + 1;
				my $SmallestLargerSun		= "";
				my $LargestSmallerSun		= "";
				my $LargestSmallerHom		= $FirstSubjectCoord; # beginning of interval
				my $SmallestLargerHom		= $LastSubjectCoord;  # end of interval
				
				foreach my $Sun (@SUNs){
							
					(undef, my $SunPosStart, my $SunPosEnd, undef) = split(m/_/, $Sun);
																			
					if 	($SunPosEnd < $VarPos){
						
						if ($SunPosEnd > $LargestSmaller){
							
							$LargestSmaller 	= $SunPosEnd;
							$LargestSmallerSun	= $Sun;
						}
					}
					elsif ($SunPosStart > $VarPos){
						
						if ($SunPosStart < $SmallestLarger){
							
							$SmallestLarger 	= $SunPosStart;
							$SmallestLargerSun	= $Sun;
						}
					}
				}
				
				my $SmallestLargerSunHom = "/";
				my $LargestSmallerSunHom = "/";
				my $TractLength			 = $SmallestLarger-$LargestSmaller-1;
				
				if ($SmallestLargerSun ne ""){
					
					$SmallestLargerSunHom = $HomoloVarMappings{$SmallestLargerSun};
					
					(undef, $SmallestLargerHom, undef) = split(m/_/, $SmallestLargerSunHom);
					
				}
				
				if ($LargestSmallerSun ne ""){
					
					$LargestSmallerSunHom = $HomoloVarMappings{$LargestSmallerSun};
					
					(undef, $LargestSmallerHom, undef) = split(m/_/, $LargestSmallerSunHom);
				}
				
				my $MinFinalDonor		= min($SmallestLargerHom,$LargestSmallerHom);
				my $MaxFinalDonor		= max($SmallestLargerHom,$LargestSmallerHom);
				
				# Compute TractAnalysisLength #
				
				my $TractAnalysisLength		= 0;
				my $LeftDiff 				= $VarPos - $LargestSmaller;
				my $RightDiff 				= $SmallestLarger - $VarPos;
				
				my $LeftTractFlank			= $LeftDiff-1;
				my $RightTractFlank			= $RightDiff-1;
				
				my $Min 					= min($LeftDiff, $RightDiff);
				my $FlankRatio				= $Min/($LeftTractFlank+1+$RightTractFlank);
				
				foreach my $RepeatLength (@RepeatLengths){ # find largest which is still smaller than repeatlength
				
					my $Flank = ($RepeatLength-1)/2;
					
					if($Flank < $Min){
						$TractAnalysisLength = $RepeatLength;
					}
				}
				
				# print #
					
				print O $VARIANT . "\t" . $FREQ . "\t" . $Chrom . "\t" . $AcceptorStart . "\t" . $AcceptorEnd  . "\t" .  $ChromName  . "\t" . $FirstSubjectCoord  . "\t" . $LastSubjectCoord . "\t" . $NrOfMatches . "\t" . $SUNs . "\t" . $HomSUNs . "\t" . $TractLengthRange. "\t" . $Strand . "\t" . $LargestSmallerHom . "\t" . $SmallestLargerHom . "\t" . $LargestSmallerSunHom . "\t" . $SmallestLargerSunHom . "\t" . $LargestSmallerSun . "\t" . $SmallestLargerSun . "\t" . $LargestSmaller . "\t" . $SmallestLarger . "\t" . $TractLength . "\t" . $TractAnalysisLength . "\t" . $LeftTractFlank . "\t" . $RightTractFlank . "\t" . $FlankRatio . "\n";
			}
		}
	}
	
	close O;
}

#######################################################################################
#######################################################################################
######     		     	SELECT LARGEST REPEAT LENGTH AND WRITE OUT 		   	 	 ######
#######################################################################################
#######################################################################################

chdir ("$TractLengthDir\/$Chrom\_$VarPos\/");
system("find $TractLengthDir/$Chrom\_$VarPos/ -name \"*txt\" | xargs cat | grep -v '^QUERY_VARIANT' > TEMP_TRACTLENGTH_TABLE");
#print "find $TractLengthDir/$Chrom\_$VarPos/ -name *txt | xargs cat | grep -v '^QUERY_VARIANT' > TEMP_TRACTLENGTH_TABLE\n";

### Iterate Over Table

open T, "TEMP_TRACTLENGTH_TABLE" or die ("Can't open TEMP_TRACTLENGTH_TABLE\n");
while(<T>){
		
	$_ =~ s/\n//g;
	
	my 	@LineValues 			= split(m/\t/, $_);
	my 	$Variant				= $LineValues[0];
	my 	$Freq					= $LineValues[1];
	my 	$TractLength			= $LineValues[21];
	my 	$TractAnalysisLength	= $LineValues[22];

	$AllTractLengths{"$Variant\~$Freq"}{$TractLength} = undef;
	$AllRepeatLengths{"$Variant\~$Freq"}{$TractAnalysisLength} = undef;
	
}
close T;

open O, ">$ResultFilePath" or die ("Can't open $ResultFilePath\n");
open T, "TEMP_TRACTLENGTH_TABLE" or die ("Can't open TEMP_TRACTLENGTH_TABLE\n");
while(<T>){
	
	$_ =~ s/\n//g;
		
	my 	@LineValues 			= split(m/\t/, $_);
	my 	$Variant				= $LineValues[0];
	my 	$Freq					= $LineValues[1];
	my 	$TractLength			= $LineValues[21];
	my 	$TractAnalysisLength	= $LineValues[22];
	
	my 	@TractLengths 		= keys %{$AllTractLengths{"$Variant\~$Freq"}};
	my 	$MaxTractLength		= max(@TractLengths);
	
	my 	@RepeatLengths 		= keys %{$AllRepeatLengths{"$Variant\~$Freq"}};
	my 	$MaxRepeatLength	= max(@RepeatLengths);
	
	if ($TractLength == $MaxTractLength || $TractAnalysisLength == $MaxRepeatLength){
	
		print O $_ . "\t" . $MaxTractLength . "\t" . $MaxRepeatLength . "\n";
	}	
}
close T;
close O;

chdir($TractLengthDir);
# system("rm -rf $Chrom\_$VarPos\/");
# if (-f "$ResultFilePath.gz"){system("rm $ResultFilePath.gz");}
# system("gzip $ResultFilePath");

print "$VARIANT has finished with success\n";

sub IntroduceMutationInSlice {
		
	(my $RefSeqSlice,
	 my $Variant,
	 my $Range)
	= @_;
	
	my $MutatedSlice = "";
	
	(my $Chrom, my $StartPos, my $EndPos, my $Alt) = split(m/_/, $Variant);
	(my $RefAllele, my $AltAllele) = split(m/\//, $Alt);
	
	my $LeftSeq 	= substr($RefSeqSlice, 0, $Range);
	my $RightSeq 	= substr($RefSeqSlice, $Range+1, $Range);
	
	$MutatedSlice	= $LeftSeq . $AltAllele . $RightSeq;
	
	return $MutatedSlice;
}
