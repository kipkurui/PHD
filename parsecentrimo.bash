#!/bin/bash
#Author: Caleb Kibet
#Date: 14.05.2013
#Sorts and parses centrimo test file output
centri=/home/Caleb/Dropbox/Projects/results/Chipseq/ENCODE/HAIB/K562/CTCF/CTCF-Rep2/centrimo_out
output=$HOME/Dropbox/Projects/results/Centrimo-output-2.txt
head -2 $centri/centrimo.txt | tail -1 |cut -f 5,3,7,10,12,2 | awk -F"\t" '{print $2,$3,$4,$5,$6,$1}' >>$output
for lab in HAIB #SYDH OpenChromUTA
do
for cl in GM78 K562 #hESC 
do
for tf in ATF3 CREB1 EGR1 YY1 #CTCF MYC SP1 MAX  JUND NRF1 TBP CEPB
do
for pk in Rep1 Rep2
do
results=$HOME/Dropbox/Projects/results/Chipseq/ENCODE/$lab/$cl/$tf/$tf-$pk
if [ -d "$rsults"]; then 
	echo "$lab-$cl-$tf-$pk" >>$output  #save the details of the TF being parsed
fi
tail -n +3 $results/centrimo_out/centrimo.txt |cut -f 5,3,7,10,12,2 | awk -F"\t" '{print $2,$3,$4,$5,$6,$1}'|sort -g -k2 | head -5 >>$output

done
done
done
done

