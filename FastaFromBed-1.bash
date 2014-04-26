#!/bin/bash
#Author: Caleb Kibet
#Date: 16.02.2013
#Script makes fasta files from bed files
#Use the for loop to execute both TFss
lab="HAIB"
#for lab in HAIB # SYDH
#do
tf="EGR1"
for cl in GM78 hESC K562
do
#for tf in CTCF NRF1 TBP JUND MAX MYC CEPB 
#do
#rep="Rep1"
for rep in Rep1 Rep2
do
results=$HOME/Dropbox/Projects/results/Chipseq/ENCODE/$lab/$cl/$tf
downloads=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
derived=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Derived
hg=$HOME/Projects/Data/genomes/hg/hg19
#cut -f 1-3  $downloads/wgEncodeHaibTfbs*Peak|bed-widen -width 500 > $derived/$tf-$rep.bed
fastaFromBed -fi "chrM.fa" -bed "test.bed" -fo /tmp/$tf-$rep.fa #$derived/$tf-$rep.bed
tr acgt N < /tmp/$tf-$rep.fa | sed 's/Nhr/chr/g' > "$rep.fa"
#done 
done
done

