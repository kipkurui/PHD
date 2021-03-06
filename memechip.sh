#!/bin/bash
# Author: Caleb Kibet
#Date: 24.08.2013
#Script makes fasta files from bed files
#Use the for loop to execute all TFss

function memechip {
	downloads=$HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
	derived=$HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Derived
	results=$HOME/Project/Results/Chipseq/ENCODE/$lab/$cl/$tf
	hg=$HOME/Project/Data/genomes/hg/hg19
	db=$HOME/meme/db/motif_databases/All-max.meme #JASPAR_CORE_2009.meme
	meme-chip -oc $results/$cl-$tf-Maxtest-trans -db $db -nmeme 500 -meme-nmotifs 5 -dreme-m 5 $derived/$tf-$cl.fa
	#centrimo --neg $neg/$tf-$pk.fa --verbosity 1 -oc $results/$tf-$pk/GM78-K562 -bgfile $results/$tf-$pk/background -score 5 -ethresh 10 $derived/$tf-$pk.fa $db $results/$tf-$pk/meme_out/meme.xml $results/$tf-$pk/dreme_out/dreme.xml
}
PATH=/home/kipkurui/meme/lib/perl/:$PATH
export PATH
file1=$HOME/Desktop/
for tf in Max # $(cat  $file1/tennew);
do
	for lab in Sydh Haib
	do 
		for cl in K562 Helas3 T47d Shsy5y Imr90 A549 Nb4 Hct116 Hek293 Panc1 Mcf7 Gm12878 H1hesc Hepg2 
		do
   			echo "doing $tf-$cl..."
			memechip
		#else
   			#echo "File $downloads/wgEncodeAwgTfbs$lab$cl$tf*narrowPeak.gz does not exists"
		#fi
			
		done
	done
done

