#!/bin/bash
#Author: Caleb Kibet
results=$HOME/Project/Results/Chipseq/ENCODE/$lab/$cl/$tf/$cl-$tf
file1=$HOME/Desktop
output=$HOME/Project/Results
#head -2 $output/centrimo.txt | tail -1 |cut -f 5,3,7,10,12,2 | awk -F"\t" '{print $2,$3,$4,$5,$6,$1}' >>$output/parsedcentrimo.txt	
for lab in Sydh Haib
do 
	for cl in Gm12878 H1hesc Hepg2 K562 Helas3 T47d Shsy5y Imr90 A549 Nb4 Hct116 Hek293 Panc1 Mcf7
	do
		for tf in $(cat  $file1/tennew);
		do
			ls 
		done
	done
done
