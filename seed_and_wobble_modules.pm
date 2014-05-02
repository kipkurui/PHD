#!/usr/bin/perl
#
# (c)2006-2008 The Brigham and Women's Hospital, Inc and President and Fellows of Harvard College
#

use warnings;
use strict;

package seed_and_wobble_modules;

use lib './';

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(
			elementbuilder
			get_kmer_truncatedarea
			truncated_enrichmentU_onearray
			convert_to_letters
			median
			get_top_X_elements
			extend_element
			wobble_seed_rerank
			get_variants
			add_element
			collapse_rankings
			get_reranked_truncated_enrichmentU
			find_minimum_info_pos
			extend_seed_rerank
			extend_seed_allpatterns_rerank
			gapped_convert_to_letters
			rc
			trim_seed
			rerank_spots_unexpected_signals
			score_sequence_against_pwm
			rerank_spots_unexpected_signals_gomer
			score_sequence_against_pwm_gomer
			collapse_dinucleotide_rankings
		     );

#################################################################
# Recursive function to build all gapped patterns of a certain order for a given input word.
# Here, $currentorder is the number of non-dot positions at the current function call, $desiredorder
# is the ultimate number of non-dot positions that are desired
#################################################################
# written by A. Philippakis, 2006
#################################################################

sub elementbuilder{
	my $word=shift; #word used in current function call
	my $currentorder=shift; #number of non-. positions in current word
	my $startorder=shift; #minimum number of non-. positions in order for it to be stored
	my $stoporder=shift; #maximum number of non-. positions in order for it to be stored
	my $lastdotpos=shift; #location of where last dot is
	my $wordlength=shift; #length of starting word
	my $spacedwordsAref=shift; #array that holds all of the spaced words

	my $i; my $j;
	my $currentword;


	if($currentorder == $stoporder){
		push @{$spacedwordsAref},$word;
	}
	else{
		if($currentorder <= $startorder){push @{$spacedwordsAref},$word;}

		my @splitword=split//,$word;

		for($i=$lastdotpos+1; $i<($wordlength-1); $i++){
			$currentword="";
			for($j=0; $j<$wordlength; $j++){
				if($i != $j){ $currentword = $currentword.$splitword[$j];}
				else{$currentword = $currentword.".";}
			}
			elementbuilder($currentword,$currentorder-1,$startorder,$stoporder,$i,$wordlength,$spacedwordsAref);
		}
	}
}

##############################################################################
# Program that takes all k-mers and their ranks and, for each k-mer,
# calculates the "area" statistic after dropping a certain fraction from
# the foreground and background.
# Recall that  area = (1/(B+F))(averagebackgroundrank-averageforegroundrank).
##############################################################################
# written by A. Philippakis, 2006
##############################################################################

sub get_kmer_truncatedarea{
	my $kmerranksHref=shift; #ranks of every kmer
	my $kmertruncareasHref=shift; #will store truncated areas for every kmer
	my $numberspots=shift; #number of spots on the array
	my $keepfraction=shift; #fraction of spots to keep.  For example 0.75 drops bottom quarter, 0.5 drops bottom half
	my $key;

	foreach $key (keys %{$kmerranksHref}){
		$kmertruncareasHref->{$key}=truncated_enrichmentU_onearray(\@{$kmerranksHref->{$key}},$numberspots,$keepfraction);
		#print $key
	}
}

##################################################################################
# Function to compute the area statistic in the case of one array, but with the
# proviso that a fraction of the bottom of the foreground and background will be dropped.
# This is based on the idea that there will be outliers in the bottom that
# one wants to avoid.  Inputs are an array (not necessarily sorted) of the foreground ranks,
# the number of data points, and the fraction of the foreground/background to consider.
##################################################################################
# written by A. Philippakis, 2006
##################################################################################

sub truncated_enrichmentU_onearray{
	my $foregroundAref=shift; #array holding foreground elements
	my $Npoints=shift; #number of data points
	my $keepfraction=shift; #fraction of data points to drop. For example, if we want to drop bottom half do 0.5, bottom quarter do 0.75.

	@{$foregroundAref} = sort {$a<=>$b} @{$foregroundAref};

	my $i; my $j; my $k;
	my $originalfgsize=$#{$foregroundAref}+1;
	my $fgsize = int($originalfgsize*$keepfraction); #foregroundsize after truncation;
	my $bgsize = int(($Npoints-$originalfgsize)*$keepfraction); #bgsize after truncation
	my $ranksum=0;
	my $result;
	
	my $lastrank=$bgsize;

	for($i=0; $i<$fgsize; $i++){
		if(($foregroundAref->[$i]-$i)>$bgsize){
			$ranksum+= ($lastrank+1);
			$lastrank++;
		}
		else{
			$ranksum += $foregroundAref->[$i];
			$lastrank++;
		}
	}
	if ($fgsize > 0) {
	    $result= (($fgsize*($fgsize+1)/2) + ($fgsize*$bgsize/2) - $ranksum)/($fgsize*$bgsize);
	}
	else {$result=0;}
	return($result);
}

#####################################################################################
# Function takes a number and a "k" and returns the ACGT-equivalent of that number
# for length k.
#####################################################################################
# written by A. Philippakis, 2006
#####################################################################################

sub convert_to_letters {
	my $number = shift;
	my $k = shift;
	my @letters = ("A", "C", "G", "T");
	my $string = "";
	my $i;
	
	for ($i = ($k-1)*2; $i >= 0; $i-=2) {
		$string .= $letters[($number >> $i) & 0x3];
	}
	return $string;
}

########################################################
# Function returns the median of an array of numbers.
########################################################

sub median
{
	my $Aref = shift; #ref to input array;
	my $median;
	my $i; my $j; my $center;
	my @tmpA; #will store the sorted $Aref;

	my $oddeven = (($#{$Aref}+1)%2);
	@tmpA = sort {$a <=> $b} @{$Aref};


	if($oddeven){
		$center = $#tmpA/2;
		$median = $tmpA[$center];
	}
	else{
		$center = int($#tmpA/2);
		$median = ($tmpA[$center]+$tmpA[$center+1])/2;
	}
	return($median);
}

###############################################################################
# Function to take the hash that stores the area for every k-mer and find the
# top X values (where "X" is the passed variable $topX).  These are stored in an
# array of hashes that is sorted by these areas.
###############################################################################
# written by A. Philippakis, 2006
###############################################################################

sub get_top_X_elements{
	my $elementHref=shift; #hash storing area of every element
	my $topX=shift; #top X elements to return
	my $topelementsAref=shift; #array 

	@{$topelementsAref}=();

	my $i;
	my $key;
	my @elementvalues; #basically, restores hash as an array of hashes for sorting

	foreach $key (keys %{$elementHref}){
		push @elementvalues, {
			element=>$key,
			value=>$elementHref->{$key}
		};
	}

 	@elementvalues = sort {$b->{value} <=> $a->{value}} @elementvalues;

	for($i=0; $i<$topX; $i++){
		push @{$topelementsAref},{
			element=>$elementvalues[$i]{element},
			value=>$elementvalues[$i]{value}
		};
	}
}

###################################################################################
# Simple function to take an element of a spaced seed and extend it by adding .'s.
# The number of positions added is such that it is possible to wobble all the
# elements up to a given max-width
###################################################################################
# written by A. Philippakis, 2006
###################################################################################

sub extend_element{
	my $element=shift; #word to extend
	my $maxseedorder=shift; #max seed order collected in this run

	my $len=length($element);
	my $iterations=$maxseedorder-$len;
	my $i;

	for($i=0; $i<$iterations; $i++){
		$element = ".".$element.".";
	}
	return($element);
}

################################################################
# Function to build a pwm based on area medians
################################################################
# written by A. Philippakis, 2006
################################################################

sub wobble_seed_rerank{
	my $kmerranksHref=shift;
	my $center=shift; #element to wobble
	my $pwmHref=shift; #ref to HofA storing PWM entries
	my $arrayfile=shift;
	my $spotlength=shift;
	my $startpos=shift;
	my $keepfraction=shift;

	my $i; my $j; my $k;
	my @splitcenter=split//,$center;
	my $Avariant; my $Cvariant; my $Gvariant; my $Tvariant;
	my $Aval; my $Cval; my $Gval; my $Tval;
	my $medianrank;
	my %collapsedranks;
	my $numbercollapsedobs;

	undef %{$pwmHref};

	for($i=0; $i<=$#splitcenter; $i++){
		if($splitcenter[$i] ne "."){

			($Avariant,$Cvariant,$Gvariant,$Tvariant) = get_variants($center,$i);

			if(!exists($kmerranksHref->{$Avariant})){
				add_element($Avariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
			}
			if(!exists($kmerranksHref->{$Cvariant})){
				add_element($Cvariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
			}
			if(!exists($kmerranksHref->{$Gvariant})){
				add_element($Gvariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
			}
			if(!exists($kmerranksHref->{$Tvariant})){
				add_element($Tvariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
			}

			$numbercollapsedobs=collapse_rankings($Avariant,$Cvariant,$Gvariant,$Tvariant,$kmerranksHref,\%collapsedranks);

			$Aval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Avariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
			$Cval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Cvariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
			$Gval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Gvariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
			$Tval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Tvariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
			$pwmHref->{A}[$i]=$Aval;
			$pwmHref->{C}[$i]=$Cval;
			$pwmHref->{G}[$i]=$Gval;
			$pwmHref->{T}[$i]=$Tval;
		}
		else{
			$pwmHref->{A}[$i]=0;
			$pwmHref->{C}[$i]=0;
			$pwmHref->{G}[$i]=0;
			$pwmHref->{T}[$i]=0;
		}	
	}	
}

###############################################################################
# Quick function to get all of the variants at a given position for an
# input word;  NOTE THAT THIS RETURNS VARIANTS THAT ARE LOWER IN LEXICOGRAPHIC
# ORDER WITH RESPECT TO REVERSE COMPLEMENTATION.  NOTE THAT THIS ALSO TRIMS
# OFF .'S THAT BEGIN AND END THE STRING
###############################################################################
# written by A. Philippakis, 2006
###############################################################################

sub get_variants{

	my $center=shift;
	my $wobbleposition=shift; #position to wobble

	my @splitcenter=split//,$center;
	my $Avariant=""; my $rcAvariant;
	my $Cvariant=""; my $rcCvariant;
	my $Gvariant=""; my $rcGvariant;
	my $Tvariant=""; my $rcTvariant; 

	my $j;
	my $startpos=$wobbleposition; my $stoppos=$wobbleposition;
	
	for($j=0; $j<=$#splitcenter;$j++){
		if($splitcenter[$j] ne "."){
			if($j<$startpos){$startpos=$j;}
			if($j>$stoppos){$stoppos=$j;}
		}
	}

	for($j=$startpos; $j<=$stoppos; $j++){
		if($wobbleposition!=$j){
			$Avariant = $Avariant.$splitcenter[$j];
			$Cvariant = $Cvariant.$splitcenter[$j];
			$Gvariant = $Gvariant.$splitcenter[$j];
			$Tvariant = $Tvariant.$splitcenter[$j];
		}
		else{
			$Avariant = $Avariant."A";
			$Cvariant = $Cvariant."C";
			$Gvariant = $Gvariant."G";
			$Tvariant = $Tvariant."T";
		}
	}

	$rcAvariant = rc($Avariant);
	$rcCvariant = rc($Cvariant);
	$rcGvariant = rc($Gvariant);
	$rcTvariant = rc($Tvariant);

	if($Avariant gt $rcAvariant){$Avariant = $rcAvariant;}
	if($Cvariant gt $rcCvariant){$Cvariant = $rcCvariant;}
	if($Gvariant gt $rcGvariant){$Gvariant = $rcGvariant;}
	if($Tvariant gt $rcTvariant){$Tvariant = $rcTvariant;}
	
	return($Avariant,$Cvariant,$Gvariant,$Tvariant);
}

#########################################################################################
# Function takes a k-mer that may have been ignored in the first pass, scans through the
# array, and adds it as needed.
#########################################################################################
# written by A. Philippakis, 2006
# modified by M. Berger (4/24/07) to discard positions closest to end of spot
#########################################################################################

sub add_element{
	my $element=shift; #element which it is desired to get the rank orderings of
	my $kmerranksHref=shift; #place where rank orderings are stored
	my $arrayfile=shift; #name of arrayfile
	my $spotlength=shift; #length of spot to consider
	my $startpos=shift; #position from end of strand to consider, starting at 1

	my $rank=1;
	my @tmp;
	my $i; my $j; my $k;
	my $lastval;	
	my $rcelement = rc($element);
	my $bool;
	my $sequence;

	open(ARRAYFILE,"<$arrayfile") || die "couldn't open $arrayfile\n";
	while(<ARRAYFILE>){
		chomp($_);
		@tmp=split/\t/,$_;
		$lastval=$#tmp;
		$sequence=substr($tmp[$lastval],$startpos-1,$spotlength-$startpos+1);

		if($sequence=~/$element/i){
			$bool=1;		
		}
		elsif($sequence=~/$rcelement/i){
			$bool=1;
		}
		else{
			$bool=0;
		}
		if($bool){
			if($element lt $rcelement){push @{$kmerranksHref->{$element}},$rank;}
			else{push @{$kmerranksHref->{$rcelement}},$rank;}
		}
		$rank++;
	}
	close(ARRAYFILE);
}

###########################################################################
# Function to take ranks of A variant, C variant G variant and T variant
# and collapse them into a hash that stores their rankings with respect to each other
###########################################################################
# written by A. Philippakis, 2006
###########################################################################

sub collapse_rankings{
	my $Avariant=shift;
	my $Cvariant=shift;
	my $Gvariant=shift;
	my $Tvariant=shift;
	my $kmerranksHref=shift; #ref to ranks
	my $collapsedranksHref=shift;

	undef %{$collapsedranksHref};
	my @observedranks;
	my $i; my $j; my $key;
	my $tester;

	for($i=0; $i<=$#{$kmerranksHref->{$Avariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$Avariant}[$i] == $observedranks[$j]){
				$tester=0;
				last;				
			}
		}
		if($tester){
			push @observedranks,$kmerranksHref->{$Avariant}[$i];
		}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$Cvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$Cvariant}[$i] == $observedranks[$j]){
				$tester=0;
				last;				
			}
		}
		if($tester){
			push @observedranks,$kmerranksHref->{$Cvariant}[$i];
		}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$Gvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$Gvariant}[$i] == $observedranks[$j]){
				$tester=0;
				last;				
			}
		}
		if($tester){
			push @observedranks,$kmerranksHref->{$Gvariant}[$i];
		}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$Tvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$Tvariant}[$i] == $observedranks[$j]){
				$tester=0;
				last;				
			}
		}
		if($tester){
			push @observedranks,$kmerranksHref->{$Tvariant}[$i];
		}	
	}	

	@observedranks = sort {$a<=>$b} @observedranks;

	for($i=0; $i<=$#observedranks; $i++){
		$collapsedranksHref->{$observedranks[$i]} = $i+1;
	}

	return($#observedranks+1);
}

############################################################################################
# Function to get the truncated area for a given variant AFTER DOING THE RERANKING PROCEDURE
############################################################################################
# written by A. Philippakis, 2006
############################################################################################

sub get_reranked_truncated_enrichmentU{
	my $variantranksAref=shift; #array that stores the ranks in original array of a given word variant
	my $collapsedranksHref=shift; #stores the ranks for the wobbles at that position
	my $numberobs=shift; #total number of data points in the foreground and background
	my $keepfraction=shift;

	my @collapsedmedian=();
	my $i; my $j=0; my $key; 
	my $returnval;

	for($i=0; $i<=$#{$variantranksAref}; $i++){
		push @collapsedmedian, $collapsedranksHref->{$variantranksAref->[$i]};
	}
	$returnval=truncated_enrichmentU_onearray(\@collapsedmedian,$numberobs,$keepfraction);
	return($returnval);
}

#################################################################################
# Function to take a nascent pwm (i.e., the entries are in areas), and find the
# position of minimum information content. Here it converts areas to a boltzman
# using a scaling factor that is passed to the function.  It then computes entropy at 
# each position using a uniform background, and returns the position that is the least info.
# Note that this also takes as input a seed that tells it which positions to ignore.
#################################################################################
# written by A. Philippakis, 2006
#################################################################################

sub find_minimum_info_pos{
	my $areapwmHref=shift; #nascent pwm that has areas as entries
	my $seed=shift;
	my $scalingfactor=shift;

	my $i; my $j; my $k;
	my $key; my $key2;
	my %tmppwm;
	my @splitseed=split//,$seed;
	my $seedlength=$#splitseed;
	my $denominator;	
	my $entropy;
	my $minentropy=5; #value greater than 2 to start it at;
	my $minentropypos;

	if($#splitseed != $#{$areapwmHref->{A}}){die "pwm not same size as input seed in find_minimum_info_pos\n";}
	if($#{$areapwmHref->{A}} != $#{$areapwmHref->{C}}){die "pwm rows do not have same number of entries in find_minimum_info_pos\n";}
	if($#{$areapwmHref->{A}} != $#{$areapwmHref->{G}}){die "pwm rows do not have same number of entries in find_minimum_info_pos\n";}
	if($#{$areapwmHref->{A}} != $#{$areapwmHref->{T}}){die "pwm rows do not have same number of entries in find_minimum_info_pos\n";}

	for($i=0; $i<=$#splitseed; $i++){
		if($splitseed[$i] =~ m/[ACGT]/i){
			$denominator=0;
			$entropy=0;
			foreach $key (keys %{$areapwmHref}){
				$tmppwm{$key}[$i] = exp($areapwmHref->{$key}[$i]*$scalingfactor);
				$denominator += $tmppwm{$key}[$i];
			}
			foreach $key (keys %{$areapwmHref}){
				$tmppwm{$key}[$i] = $tmppwm{$key}[$i]/$denominator;
				$entropy+=$tmppwm{$key}[$i]*log($tmppwm{$key}[$i])/log(2);
			}
			$entropy = 2+$entropy;
			if($entropy<$minentropy){
				$minentropy=$entropy;
				$minentropypos=$i;
			}
		}
	}	
	return($minentropypos);
}

############################################################################
# Variation on wobble_seed_rerank, but this takes a pwm that has already
# been built and extends it.
############################################################################
# written by A. Philippakis, 2006
############################################################################

sub extend_seed_rerank{
	my $kmerranksHref=shift;
	my $center=shift; #element to wobble
	my $minimuminfopos=shift; #position to now wobble in center when looking at the others
	my $pwmHref=shift; #ref to HofA storing PWM entries
	my $arrayfile=shift; #ref to arrayfile
	my $spotlength=shift;
	my $startpos=shift;
	my $keepfraction=shift; #fraction of foreground and background to keep in getting area

	my $i; my $j; my $k;
	my @splitcenter=split//,$center;
	my $newcenter="";
	my $Avariant; my $Cvariant; my $Gvariant; my $Tvariant;
	my $Aval; my $Cval; my $Gval; my $Tval;
	my $medianrank;
	my %collapsedranks;
	my $numbercollapsedobs;

	for($i=0; $i<=$#splitcenter; $i++){
		if($i!=$minimuminfopos){
			$newcenter = $newcenter.$splitcenter[$i];
		}
		else{
			$newcenter = $newcenter.".";
		}
	}

	for($i=0; $i<=$#splitcenter; $i++){
		if(($splitcenter[$i] eq ".")&&($i!=$minimuminfopos)){

			($Avariant,$Cvariant,$Gvariant,$Tvariant) = get_variants($newcenter,$i);

			if(!exists($kmerranksHref->{$Avariant})){
				add_element($Avariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
			}
			if(!exists($kmerranksHref->{$Cvariant})){
				add_element($Cvariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
			}
			if(!exists($kmerranksHref->{$Gvariant})){
				add_element($Gvariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
			}
			if(!exists($kmerranksHref->{$Tvariant})){
				add_element($Tvariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
			}

			$numbercollapsedobs=collapse_rankings($Avariant,$Cvariant,$Gvariant,$Tvariant,$kmerranksHref,\%collapsedranks);

			$Aval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Avariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
			$Cval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Cvariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
			$Gval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Gvariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
			$Tval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Tvariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
			$pwmHref->{A}[$i]=$Aval;
			$pwmHref->{C}[$i]=$Cval;
			$pwmHref->{G}[$i]=$Gval;
			$pwmHref->{T}[$i]=$Tval;
		}
	}	
}

############################################################################
# Variation on extend_seed_rerank, but this first discards the minimum
#  information position and examines all other positions corresponding to
#  spaced seeds covered in the array design.
# Requires a separate file of gapped patterns specific to that array.
############################################################################
# adapted from extend_seed_rerank by M. Berger (4/24/07)
############################################################################

sub extend_seed_allpatterns_rerank{
	my $kmerranksHref=shift;
	my $center=shift; #element to wobble
	my $minimuminfopos=shift; #position to now wobble in center when looking at the others
	my $pwmHref=shift; #ref to HofA storing PWM entries
	my $arrayfile=shift; #ref to arrayfile
	my $spotlength=shift;
	my $startpos=shift;
	my $keepfraction=shift; #fraction of foreground and background to keep in getting area
	my $patternsAref=shift; #reference to array of patterns covered [e.g. 11111111, 1.1111111, 11.111111, etc.]

	my $i; my $j; my $k;
	my @splitcenter=split//,$center;
	my $newcenter="";
	my $centerseed="";
	my $Avariant; my $Cvariant; my $Gvariant; my $Tvariant;
	my $Aval; my $Cval; my $Gval; my $Tval;
	my $medianrank;
	my %collapsedranks;
	my $numbercollapsedobs;

	for($i=0; $i<=$#splitcenter; $i++){
		if($i!=$minimuminfopos){
			$newcenter = $newcenter.$splitcenter[$i];
		}
		else{
			$newcenter = $newcenter.".";
		}
	}

	$centerseed = $newcenter;
	$centerseed =~ tr/ACGTacgt/11111111/;

	for($i=0; $i<=$#splitcenter; $i++){

		if(($splitcenter[$i] eq ".")&&($i!=$minimuminfopos)){

			my @centerseed_array = undef;
			my $query_seed = "";
			@centerseed_array = split//,$centerseed;
			$centerseed_array[$i] = "1";
			for (my $x=0; $x<=$#centerseed_array; $x++) {
				$query_seed = $query_seed.$centerseed_array[$x];
			}
			$query_seed = trim_seed($query_seed);
#			print "$query_seed\n";

			my $seed_present = "false";
			my $counter = 0;
			while ($$patternsAref[$counter]) {
				my $reverse_pattern = reverse $$patternsAref[$counter];
				if ($query_seed eq $$patternsAref[$counter]) {$seed_present = "true";}
				if ($query_seed eq $reverse_pattern) {$seed_present = "true";}
				$counter++;
			}

#			print "$seed_present\n";
			if ($seed_present eq "true") {

				($Avariant,$Cvariant,$Gvariant,$Tvariant) = get_variants($newcenter,$i);

				if(!exists($kmerranksHref->{$Avariant})){
					add_element($Avariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
				}
				if(!exists($kmerranksHref->{$Cvariant})){
					add_element($Cvariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
				}
				if(!exists($kmerranksHref->{$Gvariant})){
					add_element($Gvariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
				}
				if(!exists($kmerranksHref->{$Tvariant})){
					add_element($Tvariant,$kmerranksHref,$arrayfile,$spotlength,$startpos);
				}

				$numbercollapsedobs=collapse_rankings($Avariant,$Cvariant,$Gvariant,$Tvariant,$kmerranksHref,\%collapsedranks);

				$Aval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Avariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
				$Cval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Cvariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
				$Gval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Gvariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
				$Tval=get_reranked_truncated_enrichmentU(\@{$kmerranksHref->{$Tvariant}},\%collapsedranks,$numbercollapsedobs,$keepfraction);
				$pwmHref->{A}[$i]=$Aval;
				$pwmHref->{C}[$i]=$Cval;
				$pwmHref->{G}[$i]=$Gval;
				$pwmHref->{T}[$i]=$Tval;
			}
		}
	}	
}


######################################################################
# Very similar to convert to letters, but also takes a gapped pattern
# to generate the corresponding word
######################################################################

sub gapped_convert_to_letters
{
	my $number=shift;
	my $k=shift; 
	my $seed=shift; #gapped pattern represented as a bit string
	my $seedwidth=shift; #width of gapped pattern

	my @letters = ("A", "C", "G", "T");
	my $string = "";
	my $i;
	my $seedpositionholder=0;


	for ($i = ($k-1)*2; $i >= 0; $i-=2) {
		while(!((1<<$seedpositionholder)&$seed)){
			$string .= ".";
			$seedpositionholder += 1;
		}

		$string .= $letters[($number >> $i) & 0x3];
		$seedpositionholder += 1;
	}
  
  	return $string;
}

##########################################################
# Takes a string and returns its reverse complement
##########################################################

sub rc {
  my $string = shift;  
  $string =~ tr/acgtACGT/tgcaTGCA/;
  $string = reverse($string);
  return($string);
}

##########################################################
# Takes a spaced seed element (e.g. ..111.1..111..) and 
# trims uninformative positions from end.
##########################################################
# written by M. Berger (04/24/07)
##########################################################

sub trim_seed {
	my $fullseed = shift;
	my $return;
	my @splitseed = split//,$fullseed;
	while ($splitseed[0] eq ".") {
		shift @splitseed;
	}
	while ($splitseed[-1] eq ".") {
		pop @splitseed;
	}
	for (my $i=0; $i<=$#splitseed; $i++) {
		$return = $return.$splitseed[$i];
	}
	return $return;
}

###################################################################
# Takes a given PWM as well as a ranked list of sequences and finds
# the max-scoring match to the motif on each spot of the array.
#
# Let i be the rank (1=highest) of a given spot by signal
# intensity, and let j be the rank by sequence socre.
# We re-sort the spots in ascending order by (i/j).
###################################################################
# written by A. Philippakis, 2006; adapted by M. Berger (5/14/07)
###################################################################

sub rerank_spots_unexpected_signals{
	my $rankedseqAref = shift; #array of sequences and signals, ranked in descending order
	my $rerankedseqAref = shift; #array where output is sent after this re-ranking procedure
	my $pwmHref=shift; #pwm being used on this iteration
	my $spotlength=shift; #amount of sequence excluding primer

	my @arraydata; #array of hashes that holds sequence, observed rank, predicted rank by sequence, and observed/predicted ratio
	my $numberspots=0;
	my @tmp;
	my $sequence;
	my $sequencescore;
	my $i; my $j; my $k;
	my $key; my $key2;

	for ($i=0; $i<=$#{$rankedseqAref}; $i++) {
		$numberspots++;
		$sequence=$$rankedseqAref[$i][1];
		$sequencescore=score_sequence_against_pwm($pwmHref,$sequence,$spotlength);
		push @arraydata,{
			sequence=>$sequence,
			observedrank=>$numberspots,
			sequencescore=>$sequencescore,
			rerankvalue=>0};
	
		if(($numberspots%1000)==0){print "reviewing spot $numberspots\n";}
	}

	@arraydata = sort {$b->{sequencescore} <=> $a->{sequencescore}} @arraydata;
	for($i=0; $i<=$#arraydata; $i++){
		$arraydata[$i]{rerankvalue}=$arraydata[$i]{observedrank}/($i+1);

	}

	@arraydata = sort {$a->{rerankvalue} <=> $b->{rerankvalue}} @arraydata;

	for($i=0; $i<=$#arraydata; $i++) {
		$$rerankedseqAref[$i][0]=$arraydata[$i]{rerankvalue};
		$$rerankedseqAref[$i][1]=$arraydata[$i]{sequence};
	}
}

#####################################################
# Takes input PWM and sequence, and finds the score
# of the maximum scoring sub-sequence.
#####################################################
# Written by A. Philippakis, 2006
#####################################################

sub score_sequence_against_pwm{
	my $pwmHref=shift;
	my $sequence=shift;
	my $spotlength=shift; # amount of spot to consider after removing primer

	my $i; my $j;

	my $max; my $rcmax;
	my $sum; my $rcsum;
	my @seqarray=split//,$sequence;

	my $pwmlen;#stores length of pwm;
	if(!exists($pwmHref->{A})){die "pwm was not appropriately formatted in score_sequence_against_pwm\n";}
	else{$pwmlen=$#{$pwmHref->{A}}+1;}

	my %rcpwm; #build RC of pwm
	for($i=0; $i<$pwmlen; $i++){
		$rcpwm{A}[$i]=$pwmHref->{T}[$pwmlen-$i-1];
		$rcpwm{C}[$i]=$pwmHref->{G}[$pwmlen-$i-1];
		$rcpwm{G}[$i]=$pwmHref->{C}[$pwmlen-$i-1];
		$rcpwm{T}[$i]=$pwmHref->{A}[$pwmlen-$i-1];
	}

	for($i=($pwmlen-1);$i>=0; $i--){
		$sum=0;
		$rcsum=0;
		for($j=$i; $j<$pwmlen; $j++){
		    if(($j-$i)<$spotlength){
			$sum += $pwmHref->{$seqarray[$j-$i]}[$j];
			$rcsum += $rcpwm{$seqarray[$j-$i]}[$j];
		    }
		}
		if($i==($pwmlen-1)){
			$max=$sum;
			$rcmax=$rcsum;
		}
		else{
			if($sum>$max){$max=$sum;}
			if($rcsum>$rcmax){$rcmax=$rcsum;}
		}
	}

	for($i=0; $i<$spotlength; $i++){
		$sum=0;
		$rcsum=0;
		for($j=0; $j<$pwmlen; $j++){
			if(($i+$j)<$spotlength){
				$sum += $pwmHref->{$seqarray[$i+$j]}[$j];
				$rcsum += $rcpwm{$seqarray[$i+$j]}[$j];
			}
		}
		if($sum>$max){$max=$sum;}
		if($rcsum>$rcmax){$rcmax=$rcsum;}
	}
	if($max>$rcmax){return($max);}
	else{return($rcmax);}
}


###################################################################
# Takes a given *Probability* PWM as well as a ranked list of
# sequences and finds the max-scoring match to the motif on each spot
# of the array *Using the GOMER framework* adapted by Chen and Morris
# in Bioinformatics 2007 (RankMotif).
#
# Let i be the rank (1=highest) of a given spot by signal
# intensity, and let j be the rank by sequence socre.
# We re-sort the spots in ascending order by (i/j).
###################################################################
# written by A. Philippakis, 2006; adapted by M. Berger (9/5/07)
###################################################################

sub rerank_spots_unexpected_signals_gomer{
	my $rankedseqAref = shift; #array of sequences and signals, ranked in descending order
	my $rerankedseqAref = shift; #array where output is sent after this re-ranking procedure
	my $pwmHref=shift; #pwm being used on this iteration
	my $spotlength=shift; #amount of sequence to consider

	my @arraydata; #array of hashes that holds sequence, observed rank, predicted rank by sequence, and observed/predicted ratio
	my $numberspots=0;
	my @tmp;
	my $sequence;
	my $sequencescore;
	my $i; my $j; my $k;
	my $key; my $key2;

	for ($i=0; $i<=$#{$rankedseqAref}; $i++) {
		$numberspots++;
		$sequence=$$rankedseqAref[$i][1];
		$sequencescore=score_sequence_against_pwm_gomer($pwmHref,$sequence,$spotlength);
		push @arraydata,{
			sequence=>$sequence,
			observedrank=>$numberspots,
			sequencescore=>$sequencescore,
			rerankvalue=>0};
	
		if(($numberspots%1000)==0){print "reviewing spot $numberspots\n";}
	}

	@arraydata = sort {$b->{sequencescore} <=> $a->{sequencescore}} @arraydata;
	for($i=0; $i<=$#arraydata; $i++){
		$arraydata[$i]{rerankvalue}=$arraydata[$i]{observedrank}/($i+1);

	}

	@arraydata = sort {$a->{rerankvalue} <=> $b->{rerankvalue}} @arraydata;

	for($i=0; $i<=$#arraydata; $i++) {
		$$rerankedseqAref[$i][0]=$arraydata[$i]{rerankvalue};
		$$rerankedseqAref[$i][1]=$arraydata[$i]{sequence};
	}
}

##################################################################
# Takes input *Probability* PWM and sequence, and
# finds the score of the maximum scoring sub-sequence
# using the GOMER framework adapted by Chen and Morris
# Bioinformatics 2007 (RankMotif).
##################################################################
# Adapted by M. Berger (9/5/07) from S. Jaeger and A. Philippakis
##################################################################

sub score_sequence_against_pwm_gomer{
	my $pwmHref=shift;
	my $sequence=shift;
	my $spotlength=shift; # amount of spot to consider

	my $i; my $j;
	my $prod_gomer; my $prod_gomer_rc;
	my $value_gomer=1;

	my $max; my $rcmax;
	my $sum; my $rcsum;
	my @seqarray=split//,$sequence;

	my $pwmlen; #stores length of pwm
	if(!exists($pwmHref->{A})){die "pwm was not appropriately formatted in score_sequence_against_pwm\n";}
	else{$pwmlen=$#{$pwmHref->{A}}+1;}

	my %rcpwm; #build RC of pwm
	for($i=0; $i<$pwmlen; $i++){
		$rcpwm{A}[$i]=$pwmHref->{T}[$pwmlen-$i-1];
		$rcpwm{C}[$i]=$pwmHref->{G}[$pwmlen-$i-1];
		$rcpwm{G}[$i]=$pwmHref->{C}[$pwmlen-$i-1];
		$rcpwm{T}[$i]=$pwmHref->{A}[$pwmlen-$i-1];
	}

	for($i=($pwmlen-1); $i>0; $i--){
		$prod_gomer=1;
		$prod_gomer_rc=1;
		for($j=0; $j<$pwmlen; $j++){
			if ($j<$i) {
				$prod_gomer = $prod_gomer*0.25;
				$prod_gomer_rc = $prod_gomer_rc*0.25;
			}
			elsif ($j-$i>$#seqarray) {
				$prod_gomer = $prod_gomer*0.25;
				$prod_gomer_rc = $prod_gomer_rc*0.25;
			}
			else {
				$prod_gomer = $prod_gomer*$pwmHref->{$seqarray[$j-$i]}[$j];
				$prod_gomer_rc = $prod_gomer_rc*$rcpwm{$seqarray[$j-$i]}[$j];
			}
		}
		$value_gomer = $value_gomer*(1-$prod_gomer)*(1-$prod_gomer_rc);
	}

	for($i=0; $i<$spotlength; $i++){
		$prod_gomer=1;
		$prod_gomer_rc=1;
		for($j=0; $j<$pwmlen; $j++){
			if ($j+$i>=$spotlength) {
				$prod_gomer = $prod_gomer*0.25;
				$prod_gomer_rc = $prod_gomer_rc*0.25;
			}
			else {
				$prod_gomer = $prod_gomer*$pwmHref->{$seqarray[$j+$i]}[$j];
				$prod_gomer_rc = $prod_gomer_rc*$rcpwm{$seqarray[$j+$i]}[$j];
			}
		}
		$value_gomer = $value_gomer*(1-$prod_gomer)*(1-$prod_gomer_rc);
	}

	$value_gomer = 1 - $value_gomer;
	return ($value_gomer);
}


###########################################################################
# Function to take ranks of 16 dinucleotide variants and collapse
# them into a hash that stores their rankings with respect to each other
###########################################################################
# adapted by M. Berger (06/11/07) from A. Philippakis, 2006
###########################################################################

sub collapse_dinucleotide_rankings{
	my $AAvariant=shift;
	my $ACvariant=shift;
	my $AGvariant=shift;
	my $ATvariant=shift;
	my $CAvariant=shift;
	my $CCvariant=shift;
	my $CGvariant=shift;
	my $CTvariant=shift;
	my $GAvariant=shift;
	my $GCvariant=shift;
	my $GGvariant=shift;
	my $GTvariant=shift;
	my $TAvariant=shift;
	my $TCvariant=shift;
	my $TGvariant=shift;
	my $TTvariant=shift;

	my $kmerranksHref=shift; #ref to ranks
	my $collapsedranksHref=shift;

	undef %{$collapsedranksHref};
	my @observedranks;
	my $i; my $j; my $key;
	my $tester;

	for($i=0; $i<=$#{$kmerranksHref->{$AAvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$AAvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$AAvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$ACvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$ACvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$ACvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$AGvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$AGvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$AGvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$ATvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$ATvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$ATvariant}[$i];}	
	}	

	for($i=0; $i<=$#{$kmerranksHref->{$CAvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$CAvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$CAvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$CCvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$CCvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$CCvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$CGvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$CGvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$CGvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$CTvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$CTvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$CTvariant}[$i];}	
	}	

	for($i=0; $i<=$#{$kmerranksHref->{$GAvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$GAvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$GAvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$GCvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$GCvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$GCvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$GGvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$GGvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$GGvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$GTvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$GTvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$GTvariant}[$i];}	
	}	

	for($i=0; $i<=$#{$kmerranksHref->{$TAvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$TAvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$TAvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$TCvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$TCvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$TCvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$TGvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$TGvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$TGvariant}[$i];}	
	}	
	for($i=0; $i<=$#{$kmerranksHref->{$TTvariant}}; $i++){
		$tester=1;
		for($j=0; $j<=$#observedranks; $j++){
			if($kmerranksHref->{$TTvariant}[$i] == $observedranks[$j]){$tester=0; last;}
		}
		if($tester){push @observedranks,$kmerranksHref->{$TTvariant}[$i];}	
	}	

	@observedranks = sort {$a<=>$b} @observedranks;

	for($i=0; $i<=$#observedranks; $i++){
		$collapsedranksHref->{$observedranks[$i]} = $i+1;
	}

	return($#observedranks+1);
}
