#!/bin/bash
#A collection of perl scripts to generate PWM.

for kmer in  $(cat Max_3863.txt);
do
	perl wobble_single_seed_twoarray.pl Max_3863.1_v1_deBruijn.txt Max_3863.1_v1_deBruijn.txt $kmer patterns_4x44k_all_8mer.txt >>Max-$kmer-0.467-269-2.uniprobe
	
done

