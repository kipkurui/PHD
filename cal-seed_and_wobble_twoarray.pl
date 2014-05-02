#!/usr/bin/perl
#
# (c)2006-2008 The Brigham and Women's Hospital, Inc and President and Fellows of Harvard College
#
use strict;
use warnings;

use lib './';
use seed_and_wobble_modules;
use Data::Dumper;

##################################################################
### seed_and_wobble_twoarray.pl
###
###
### Takes as input TWO lists of all "combinatorial" probe intensities
###   and sequences: one from design v1 and one from design v2.
###   MUST BE PRE-ORDERED FROM BRIGHTEST TO DIMMEST.
### Also takes file of gapped patterns (e.g., 11111.111) to
###   consider as candidate seeds.
### Also takes file of all gapped patterns (e.g., 111..11.1.11) that
###   are covered evenly on the universal array design.  This may be
###   the same as or longer than the above file of candidate seeds.
###
### Outputs list of all possible k-mers and their corresponding
###   average enrichment score (separate file for each seed).
### Also outputs single integrated file for all seeds with an average
###   E-score above a cutoff.
### Also outputs file of top "N" position weight matrices from *combined*
###   data.
###
### M. Berger 04/04/08
##################################################################

if($#ARGV !=6){die "Usage is:\n\tPBM data (sorted by intensities) FOR DESIGN v1\n\tPBM data (sorted by intensities) FOR DESIGN v2\n\twidth of k-mers to inspect\n\tlist of candidate seed patterns (e.g., 11111.111)\n\tlist of all covered gapped patterns (e.g., 111..11.1.11)\n\toutput file prefix\nFor example, \"perl seed_and_wobble_twoarray.pl TF_v1_combinatorial.txt TF_v2_combinatorial.txt 8 query_patterns.txt all_patterns.txt output_prefix\"\n\n";}

my $intensity_file_v1=shift; #PBM data file FOR DESIGN v1, two columns (intensity, sequence), pre-sorted from brightest to dimmest
my $intensity_file_v2=shift; #PBM data file FOR DESIGN v2, two columns (intensity, sequence), pre-sorted from brightest to dimmest
my $order=shift; #k-mer order to use for seed patterns
my $seed_pattern_file=shift; #spaced seeds covered on array for seed determination
my $total_pattern_file=shift; #all spaced seeds covered on array
my $output_prefix=shift; #prefix for output files
my $DNasefile=shift; #contain DNase file to count

##############################################################################
### PARAMETERS
my $spotlength=36; # total number of nucleotides after common primer
my $startposition=2; # position from end of strand to consider, starting at 1
my $Escore_cutoff=0.25; # to store for integrated list of top k-mers from all seeds
my $print_seeds = "yes"; # "yes" if output file for each seed is desired
my $print_top_seeds = "yes"; # "yes" if output file for top seeds (above cutoff) is desired
my $topN = 3; # number of top k-mer seeds to use for seed-and-wobble PWM construction

###########################################################################
### Read in list of gapped patterns to consider for seeds
###########################################################################

my @seed_patterns;

open (FH1, "<$seed_pattern_file") or die "Cannot open seed patterns file.\n";

while (my $text = <FH1>) {
	my $numinfopos=0;
	chomp $text;
	push @seed_patterns, $text;	
	my @characters = split//,$text;
	for (my $p=0; $p<=$#characters; $p++) {
		if ($characters[$p] ne ".") {$numinfopos++;}
	}
	if ($numinfopos != $order) {
		die "Number of positions in seed $text in $seed_pattern_file does not agree with order = $order\n";
	}
}	

close FH1;

#################################################################
### Read in intensities and create arrays of sequences
#################################################################

my @data_matrix_v1;
my @data_matrix_v2;

open (FHA, "<$intensity_file_v1") or die "Cannot open intensity file #1.\n";
my $spot = 1;
while (my $text = <FHA>) {
    chomp $text;
    my @line = split ("\t", $text);
    $data_matrix_v1[$spot][0] = $line[0]; ### spot intensity
    $data_matrix_v1[$spot][1] = $line[1]; ### spot sequence
    if ($spot > 1) {
	if ($data_matrix_v1[$spot][0] > $data_matrix_v1[$spot-1][0]) {
	    die "Probes in $intensity_file_v1 are not sorted from brightest to dimmest.\n";
	}
    }
    $spot++;
}
close FHA;

open (FHB, "<$intensity_file_v2") or die "Cannot open intensity file #2.\n";
$spot = 1;
while (my $text = <FHB>) {
    chomp $text;
    my @line = split ("\t", $text);
    $data_matrix_v2[$spot][0] = $line[0]; ### spot intensity
    $data_matrix_v2[$spot][1] = $line[1]; ### spot sequence
    if ($spot > 1) {
        if ($data_matrix_v2[$spot][0] > $data_matrix_v2[$spot-1][0]) {
            die "Probes in $intensity_file_v2 are not sorted from brightest to dimmest.\n";
        }
    }
    $spot++;
}
close FHB;

###############################################################################
### Calculate median intensity and enrichment score and Z score for each 8-mer
###############################################################################

my $numberspotsarray_v1=$#data_matrix_v1; #number of spots on array 1
my $numberspotsarray_v2=$#data_matrix_v2; #number of spots on array 2
my $tally=0;
my %top_kmerareas_combined = ();
my %kmerranks_v1;
my %kmerranks_v2;
my $outputfile1;
my $keepfraction=0.5;

while ($tally<=$#seed_patterns) {

	if ($print_seeds eq "yes" || $tally ==0) {
	    if ($print_seeds eq "yes") {
		$outputfile1=$output_prefix."_".$order."mers_".$seed_patterns[$tally]."_combined.txt";
	    }
	    else {$outputfile1=$output_prefix."_".$order."mers_combined.txt";}
	    open(OUTPUT1,">$outputfile1") || die "Cannot create k-mers output file.\n";
	    my $toplabel = $order."-mer";
	    print OUTPUT1 "$toplabel\t$toplabel\tE-score\n";
	}

# array design v1

	%kmerranks_v1 = (); #data structure that stores ranks of every k-mer
	my %kmerintensities_v1 = (); #data structure that stores intensities of every k-mer
	my %kmerareas_v1 = (); #stores area for each gapped k-mer

	my %observedpatterns = ();
	my $key;
	my $spaced_seed = $seed_patterns[$tally];
	my @split_seed = split//,$spaced_seed;
	my $rev_spaced_seed = reverse $spaced_seed;
	my @split_rev_seed = split//,$rev_spaced_seed;
	my $spaced_seed_width = length $spaced_seed;

	print "Currently on $intensity_file_v1, spaced seed: $spaced_seed\n";

	for (my $spotnumber=1; $spotnumber<=$#data_matrix_v1; $spotnumber++) {
	    for (my $i=$startposition; $i<=$spotlength-$spaced_seed_width+1; $i++) {
		    my $currentstring = substr($data_matrix_v1[$spotnumber][1], $i-1, $spaced_seed_width);
		    my @splitstring = split//,$currentstring;
		    my $fwd_element="";
		    my $rev_element="";
		    my $counter;
		    for ($counter=0; $counter<$spaced_seed_width; $counter++) {
			if ($split_seed[$counter] eq "1") {
				$fwd_element = $fwd_element.$splitstring[$counter];
			}
			else {
				$fwd_element = $fwd_element.".";
			}
		    }

		    my $rc_fwd_element = rc($fwd_element);
		    if ($fwd_element lt $rc_fwd_element) {$observedpatterns{$fwd_element}=0;}
		    else {$observedpatterns{$rc_fwd_element}=0;}

		    if ($spaced_seed ne $rev_spaced_seed) {
			for ($counter=0; $counter<$spaced_seed_width; $counter++) {
				if ($split_rev_seed[$counter] eq "1") {
					$rev_element = $rev_element.$splitstring[$counter];
				}
				else {
					$rev_element = $rev_element.".";
				}
			}

			my $rc_rev_element = rc($rev_element);
			if ($rev_element lt $rc_rev_element) {$observedpatterns{$rev_element}=0;}
			else {$observedpatterns{$rc_rev_element}=0;}
		    }			
	    }
	    foreach $key (keys %observedpatterns) {
		    push @{$kmerranks_v1{$key}},$spotnumber;
		    push @{$kmerintensities_v1{$key}},$data_matrix_v1[$spotnumber][0];
	    }
	    undef %observedpatterns;
	    %observedpatterns=();
	    if ($spotnumber%1000==0) {print "$spotnumber\n";}
	}

	get_kmer_truncatedarea(\%kmerranks_v1,\%kmerareas_v1,$numberspotsarray_v1,$keepfraction);

# array design v2

	%kmerranks_v2 = (); #data structure that stores ranks of every k-mer
	my %kmerintensities_v2 = (); #data structure that stores intensities of every k-mer
	my %kmerareas_v2 = (); #stores area for each gapped k-mer

	print "Currently on $intensity_file_v2, spaced seed: $spaced_seed\n";

	for (my $spotnumber=1; $spotnumber<=$#data_matrix_v2; $spotnumber++) {
	    for (my $i=$startposition; $i<=$spotlength-$spaced_seed_width+1; $i++) {
		    my $currentstring = substr($data_matrix_v2[$spotnumber][1], $i-1, $spaced_seed_width);
		    my @splitstring = split//,$currentstring;
		    my $fwd_element="";
		    my $rev_element="";
		    my $counter;
		    for ($counter=0; $counter<$spaced_seed_width; $counter++) {
			if ($split_seed[$counter] eq "1") {
				$fwd_element = $fwd_element.$splitstring[$counter];
			}
			else {
				$fwd_element = $fwd_element.".";
			}
		    }

		    my $rc_fwd_element = rc($fwd_element);
		    if ($fwd_element lt $rc_fwd_element) {$observedpatterns{$fwd_element}=0;}
		    else {$observedpatterns{$rc_fwd_element}=0;}

		    if ($spaced_seed ne $rev_spaced_seed) {
			for ($counter=0; $counter<$spaced_seed_width; $counter++) {
				if ($split_rev_seed[$counter] eq "1") {
					$rev_element = $rev_element.$splitstring[$counter];
				}
				else {
					$rev_element = $rev_element.".";
				}
			}

			my $rc_rev_element = rc($rev_element);
			if ($rev_element lt $rc_rev_element) {$observedpatterns{$rev_element}=0;}
			else {$observedpatterns{$rc_rev_element}=0;}
		    }			
	    }
	    foreach $key (keys %observedpatterns) {
		    push @{$kmerranks_v2{$key}},$spotnumber;
		    push @{$kmerintensities_v2{$key}},$data_matrix_v2[$spotnumber][0];
	    }
	    undef %observedpatterns;
	    %observedpatterns=();
	    if ($spotnumber%1000==0) {print "$spotnumber\n";}
	}

	get_kmer_truncatedarea(\%kmerranks_v2,\%kmerareas_v2,$numberspotsarray_v2,$keepfraction);

	my $bit_seed = $spaced_seed;
	$bit_seed =~ tr/./0/;
	$bit_seed = reverse $bit_seed;
	my $decimal_bit_seed = oct( "0b$bit_seed" );

	for (my $k=0; $k<4**($order); $k++) {
		my $word = gapped_convert_to_letters ($k, $order, $decimal_bit_seed, $spaced_seed_width);
		my $revcomp = rc($word);
		if ($word le $revcomp || $spaced_seed ne $rev_spaced_seed) {
			if ($print_seeds eq "yes" || $tally == 0) {print OUTPUT1 "$word\t$revcomp\t";}
			if ($word gt $revcomp) {$word = $revcomp;}
			if ($kmerareas_v1{$word} && $kmerareas_v2{$word}) {
			    my $avg_escore = ($kmerareas_v1{$word}+$kmerareas_v2{$word})/2;
			    if ($print_seeds eq "yes" || $tally == 0) {
				printf OUTPUT1 "%.5f\n", $avg_escore;
			    }
			    if ($avg_escore > $Escore_cutoff) {
				$top_kmerareas_combined{$word}=$avg_escore;
			    }
			}
			else {
			    if ($print_seeds eq "yes" || $tally == 0) {
				print OUTPUT1 "NA\n";
			    }
			}
		}
	}
	undef %kmerranks_v1;
	undef %kmerranks_v2;
	undef %kmerareas_v1;
	undef %kmerareas_v2;

	close (OUTPUT1);

	$tally++;
}
################################################################################
# Read the DNAse data and count for the occurance of the kmers\
###############################################################################
#my %kmerareas=();

#open (FH, "<$kmerfile") or die "Cannot open DNase file.\n";

#while (my $text = <FH>) {
#    chomp $text;
#    my @line = split ("\t", $text);
#    $kmerareas{$line[0]} = $line[2]; ### spot intensity
#	}
	
#close FH;
#print Dumper(\%kmerareas);
#####################################################
#Count the kmer frequency
#####################################################
my %counts=();
my $count=0;
my $Track=0;
open (FH1, "<$DNasefile") or die "Cannot open DNAse file.\n";
#while ($Track <1000) {
foreach my $key (keys %top_kmerareas_combined) {
	if ($Track==1000) {last}
	seek(FH1, $DNasefile, 0);
	while (my $line = <FH1>) {
		chomp $line;
		while ($line =~ /$key/g) { $count++ }
		}
	$counts{$key}=$count;
	$count=0;
	$Track++;
	print "$Track\n";
	#print Dumper(\%counts);
	}
close FH1;

#####################################################
#Calculate the frequency of each kmer
####################################################
my $total=0;
foreach my $value (values %counts){
	$total+=$value;
	print "$total \n";
}
#print $total;
my %newEscore=();
foreach my $key (keys %counts){
	#get frequency of each kmer
	$newEscore{$key}=($top_kmerareas_combined{$key})*(($counts{$key})/($total)); 
}
print Dumper(\%newEscore);

####################################################
#Sort the new Escore file then use as the top_kmer areas
#######################################################
foreach my $name (sort { $newEscore{$a} <=> $newEscore{$b} or $a cmp $b } keys %newEscore) {
    #printf "%-8s %s\n", $name, $newEscore{$name};
}
#################################################
# Find top N seeds (adapted from A. Philippakis)
#################################################

print "Finding top $topN seeds.\n";

my @top_N_elements = ();
my @elementvalues;
my $key;

my $topNcounter=0;
foreach $key (keys %newEscore) {
    push @elementvalues, {
	element=>$key,
	value=>$newEscore{$key},
    };
    $topNcounter++;
}
if ($topNcounter < $topN) {$topN = $topNcounter;}

@elementvalues = sort {$b->{value} <=> $a->{value}} @elementvalues;

for (my $N=0; $N<$topN; $N++) {
    push @top_N_elements, {
	element=>$elementvalues[$N]{element},
	value=>$elementvalues[$N]{value}
    };
}

if ($print_top_seeds eq "yes") {
    my $outputfile2 = $output_prefix."_".$order."mers_top_enrichment_combined.txt";
    open(OUTPUT2,">$outputfile2") || die "Cannot create top k-mers output file.\n";
    my $toplabel = $order."-mer";
    print OUTPUT2 "$toplabel\t$toplabel\tE-score\n";
    my $N=0;
    while ($elementvalues[$N]) {
	my $rc_element = rc($elementvalues[$N]{element});
	printf OUTPUT2 "$elementvalues[$N]{element}\t$rc_element\t%.5f\n", $elementvalues[$N]{value};
	$N++;
    }
    close (OUTPUT2);
}


#########################################################################################
### Read in list of all gapped patterns covered in universal PBM design (for wobble step)
#########################################################################################

my @total_patterns;

open (FH2, "<$total_pattern_file") or die "Cannot open total patterns file.\n";

while (my $text = <FH2>) {
	my $numinfopos=0;
	chomp $text;
	push @total_patterns, $text;	
	my @characters = split//,$text;
	for (my $p=0; $p<=$#characters; $p++) {
		if ($characters[$p] ne ".") {$numinfopos++;}
	}
	if ($numinfopos != $order) {
		die "Number of positions in seed $text in $total_pattern_file does not agree with order = $order\n";
	}
}	

close FH2;


###########################################################################
### Seed-and-Wobble PWM construction
###########################################################################

my %areapwm_v1;
my %areapwm_v2;
my $array_v1=$intensity_file_v1;
my $array_v2=$intensity_file_v2;

my $outputfile3 = $output_prefix."_".$order."mers_pwm_combined.txt";
open(OUTPUT3,">$outputfile3") || die "Cannot create pwm output file.\n";

for (my $z=0; $z<=$#top_N_elements; $z++) {
	my $ranking = $z+1;
	print "Currently on element ranked:\t$ranking\n";
	my $seed = $top_N_elements[$z]{element};
	print "$seed\t";
	my $topescore = $top_N_elements[$z]{value};
	print "$topescore\n";
	print OUTPUT3 "$ranking\t$seed\t$newEscore{$seed}\n\n";
	$seed = ".......".$seed.".......";

# array design v1

	wobble_seed_rerank(\%kmerranks_v1,$seed,\%{$areapwm_v1{$seed}},$array_v1,$spotlength,$startposition,1);
	my $minimuminfopos_v1=find_minimum_info_pos(\%{$areapwm_v1{$seed}},$seed,log(10000));
	extend_seed_allpatterns_rerank(\%kmerranks_v1,$seed,$minimuminfopos_v1,\%{$areapwm_v1{$seed}},$array_v1,$spotlength,$startposition,1,\@total_patterns);

#array design v2

        wobble_seed_rerank(\%kmerranks_v2,$seed,\%{$areapwm_v2{$seed}},$array_v2,$spotlength,$startposition,1);
        my $minimuminfopos_v2=find_minimum_info_pos(\%{$areapwm_v2{$seed}},$seed,log(10000));
        extend_seed_allpatterns_rerank(\%kmerranks_v2,$seed,$minimuminfopos_v2,\%{$areapwm_v2{$seed}},$array_v2,$spotlength,$startposition,1,\@total_patterns);

#must have data for both PWM-v1 and PWM-v2 before averaging

      	while (($areapwm_v1{$seed}{A}[0]==0) || ($areapwm_v2{$seed}{A}[0]==0)) {
	    shift @{$areapwm_v1{$seed}{A}}; shift @{$areapwm_v1{$seed}{C}}; shift @{$areapwm_v1{$seed}{G}}; shift @{$areapwm_v1{$seed}{T}};
            shift @{$areapwm_v2{$seed}{A}}; shift @{$areapwm_v2{$seed}{C}}; shift @{$areapwm_v2{$seed}{G}}; shift @{$areapwm_v2{$seed}{T}};
	}

	while (($areapwm_v1{$seed}{A}[-1]==0) || ($areapwm_v2{$seed}{A}[-1]==0)) {
	    pop @{$areapwm_v1{$seed}{A}}; pop @{$areapwm_v1{$seed}{C}}; pop @{$areapwm_v1{$seed}{G}}; pop @{$areapwm_v1{$seed}{T}};
            pop @{$areapwm_v2{$seed}{A}}; pop @{$areapwm_v2{$seed}{C}}; pop @{$areapwm_v2{$seed}{G}}; pop @{$areapwm_v2{$seed}{T}};
	}

#average matrix elements

	my %areapwm_combined;

	print OUTPUT3 "Enrichment score matrix\n\n";
	foreach $key (sort keys %{$areapwm_v1{$seed}}) {
		print OUTPUT3 "$key:";
		for (my $y=0; $y<=$#{$areapwm_v1{$seed}{$key}}; $y++) {
		    $areapwm_combined{$seed}{$key}[$y] = ($areapwm_v1{$seed}{$key}[$y]+$areapwm_v2{$seed}{$key}[$y])/2;
		    print OUTPUT3 "\t$areapwm_combined{$seed}{$key}[$y]";
		}
		print OUTPUT3 "\n";
	}
	print OUTPUT3 "\nEnergy matrix for enoLOGOS\n\n";
	print OUTPUT3 "PO";
	for (my $counter=1; $counter<=($#{$areapwm_combined{$seed}{A}}+1); $counter++) {
		print OUTPUT3 "\t$counter";
	}
	print OUTPUT3 "\n";
	foreach $key (sort keys %{$areapwm_combined{$seed}}) {
		print OUTPUT3 "$key:";
		for (my $y=0; $y<=$#{$areapwm_combined{$seed}{$key}}; $y++) {
		    my $logscaled = $areapwm_combined{$seed}{$key}[$y]*(-log(10000));
		    print OUTPUT3 "\t$logscaled";
		}
		print OUTPUT3 "\n";
	}
#	print OUTPUT3 "\nReverse complement matrix for enoLOGOS\n\n";
#	print OUTPUT3 "PO";
#	for (my $counter=1; $counter<=($#{$areapwm{$seed}{A}}+1); $counter++) {
#		print OUTPUT3 "\t$counter";
#	}
#	print OUTPUT3 "\n";
#	foreach $key (sort keys %{$areapwm{$seed}}) {
#		my $compkey;
#		if ($key eq "A") {$compkey="T";}
#		if ($key eq "C") {$compkey="G";}
#		if ($key eq "G") {$compkey="C";}
#		if ($key eq "T") {$compkey="A";}
#		print OUTPUT3 "$compkey:";
#		for (my $y=$#{$areapwm{$seed}{$key}}; $y>=0; $y--) {
#		    my $logscaled = $areapwm{$seed}{$key}[$y]*(-log(10000));
#		    print OUTPUT3 "\t$logscaled";
#		}
#		print OUTPUT3 "\n";
#	}
	print OUTPUT3 "\nProbability matrix\n\n";
	foreach $key (sort keys %{$areapwm_combined{$seed}}) {
		print OUTPUT3 "$key:";
		for (my $y=0; $y<=$#{$areapwm_combined{$seed}{$key}}; $y++) {
			my $numerator = exp(log(10000)*$areapwm_combined{$seed}{$key}[$y]);
			my $denominator = exp(log(10000)*$areapwm_combined{$seed}{A}[$y]) + exp(log(10000)*$areapwm_combined{$seed}{C}[$y]) + exp(log(10000)*$areapwm_combined{$seed}{G}[$y]) + exp(log(10000)*$areapwm_combined{$seed}{T}[$y]);
			my $probability = $numerator/$denominator;
			print OUTPUT3 "\t$probability";
		}
		print OUTPUT3 "\n";
	}
	print OUTPUT3 "\n\n\n";
}

close (OUTPUT3);
