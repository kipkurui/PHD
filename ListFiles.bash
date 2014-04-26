#!/bin/bash
#Author: Caleb Kibet
#Date: 16.02.2013
#Script makes fasta files from bed files
#Use the for loop to execute both TFss
for lab in HAIB SYDH OpenChromUTA
do
for cl in GM78 hESC K562
do
for tf in ATF3 CTCF MYC CREB1 EGR1 YY1 SP1 MAX  JUND NRF1 TBP CEPB
do

results=$HOME/Dropbox/Projects/results/Chipseq/ENCODE/$lab/$cl/$tf
downloads=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
derived=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Derived
hg=$HOME/Projects/Data/genomes/hg/hg19

ls $downloads/wgEncode*Peak >>TFS2.txt
done
done
done
