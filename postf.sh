#!/bin/bash
#Author: Caleb Kibet
#Date: 14.05.2013
#Sorts and parses centrimo test file output
function parsecentrimo {
	results=$HOME/Project/Results/Chipseq/ENCODE/$lab/$cl/$tf/$cl-$tf-Maxtest
	if [ -d "$results/centrimo_out" ]; then 
		#echo "$lab-$cl-$tf" >>$output/Maxparsed.txt  #save the details of the TF being parsed
		#tail -n +3 $results/centrimo_out/centrimo.txt |cut -f 5,3,7,10,12,2 | awk -F"\t" '{print $2,$3,$4,$5,$6,$1}'|sort -g -k2 | head -5 >>$output/parsedcentrimo.txt
		#echo "$tf-$lab-$cl-${mots}"  >>$output/Maxparsed.txt
		tail -n +3 $results/centrimo_out/centrimo.txt |cut -f 5,3,7,10,12,2 | awk -F"\t" '{print $2,$3,$4,$5,$6,$1}'|sort -g -k2 | grep -n "$mots"|cut -f 1 -d ":" >>$output/Maxparsed2.txt
		#tail -n +3 $results/centrimo_out/centrimo.txt |cut -f 5,3,7,10,12,2 | awk -F"\t" '{print $2,$3,$4,$5,$6,$1}'|sort -g -k2 | grep -n "$tf"|cut -f 1 -d ":" >>$output/mymotifspos.txt
	fi
}

file1=$HOME/Desktop
output=$HOME/Project/Max
#head -2 $output/centrimo.txt | tail -1 |cut -f 5,3,7,10,12,2 | awk -F"\t" '{print $2,$3,$4,$5,$6,$1}' >>$output/parsedcentrimo.txt	
for lab in Sydh Haib
do 
	for cl in Gm12878 H1hesc Hepg2 K562 Helas3 T47d Shsy5y Imr90 A549 Nb4 Hct116 Hek293 Panc1 Mcf7
	do
		echo "$lab-$cl" >>$output/Maxparsed2.txt
		for tf in Max #$(cat  $file1/tennew);
		do
			for mots in $(cat $output/Max-uni.txt);
			do
				parsecentrimo
			done
		done
	done
done
