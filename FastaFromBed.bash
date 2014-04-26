#!/bin/bash
#Author: Caleb Kibet
#Date: 16.02.2013
#Script makes fasta files from bed files
#Use the for loop to execute all TFs

function frombed {
	#initialize folder names
	downloads=$HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
	derived=$HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Derived
	hg=$HOME/Project/Data/genomes/hg19
	gunzip $downloads/wgEncodeAwgTfbs$lab$cl$tf*narrowPeak.gz
	cut -f 1-3 $downloads/wgEncodeAwgTfbs$lab$cl$tf*narrowPeak | bed-widen -width 500 > /tmp/$tf-$cl.bed
	sort -u -o $derived/$tf-$cl.bed /tmp/$tf-$cl.bed #remove duplicates
	fastaFromBed -fi "$hg/hg19.fa" -bed $derived/$tf-$cl.bed -fo $derived/$tf-$cl.fa
}
file1=$HOME/Desktop/
for lab in Sydh Haib
do 
	for cl in Gm12878 H1hesc Hepg2 K562 Helas3 T47d Shsy5y Imr90 A549 Nb4 Hct116 Hek293 Panc1 Mcf7
	do
		#echo $cl
		for tf in $(cat  $file1/newtfs);
		do
		#if [ -f ]; #$downloads/wgEncodeAwgTfbs$lab$cl$tf*narrowPeak.gz 
		#then
   			echo "doing $tf-$cl..."
			frombed
		#else
   			#echo "File $downloads/wgEncodeAwgTfbs$lab$cl$tf*narrowPeak.gz does not exists"
		#fi
			
		done
	done
done

