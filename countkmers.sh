#!/bin/bash

#Caleb Kibet

#Count the numbert of occurance of a kmer in the DNase sites

file=/home/kipkurui/PHD-project/Data/k-mers/Human/Combined #/home/kipkurui/PHD-project/kmerprobs

file2=/home/kipkurui/PHD-project
file3=/home/kipkurui/PHD-project/Results/Human

#ls $file >$file2/listTFs.txt
for line in $(cat  $file2/listTFs.txt);
do
	cut -f 1 $file/$line >$file2/Zscan4-2.txt #these can go to temp later	
	for word in $(cat  $file2/Zscan4-2.txt); 
	do
		grep -c $word /home/kipkurui/PHD-project/Data/DNase/Derived/DNaseclustured.fa >> $file2/count.txt;
#count for the occurance of forward motif, assuming reverse has a similar score 
	done

	paste -d '\t'  $file/$line $file2/count.txt > $file3/$line
#combine both files

	rm $file2/count.txt

done

#for line in $(cat $file2/counts.txt);
#do
	
#awk -F"\t" '{print $2,$3,$4,$5,$6,$1}'
#for word in  $(cat $file2/counts.txt);
#do
#	awk '{print $0, $word }' $file/Zscan4.txt;
#done
#for line in $(cat $f1); do for l2 in $(cat $f2); do grep -o '\($line\|$l2\)' DNaseclustured.fa|wc -l; done done

#for line in $(cat  $file);
#do
	#awk '/`cut -f 1 $file`/ || /`cut -f 2 $file`/ {count++} END {print count}' /home/kipkurui/PHD-project/Data/DNase/Derived/DNaseclustured.fa
#done
