#!/bin/bash
# Author: Caleb Kibet
lab="HAIB"
#for lab in  #SYDH
#do
for cl in GM78 #K562
do
#tf="MYC"
for tf in ATF3
do
#pk="Rep1"
for pk in Rep1 Rep2
do
results=$HOME/Dropbox/Projects/results/Chipseq/ENCODE/$lab/$cl/$tf
downloads=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
derived=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Derived
db=$HOME/meme/db/motif_databases/JASPAR_CORE_2009.meme

ame --oc $results/$tf-$pk/ame_out --scoring totalhits --bgformat 0 --bgfile $derived/$tf-$pk.fa $db $results/$tf-$pk/meme_out/meme.xml $results/$tf-$pk/dreme_out/dreme.xml
done
done 
done
#done

