#!/bin/bash
# Author: Caleb Kibet
lab="SYDH"
#for lab in  #SYDH
#do
for cl in NB4
do
#tf="MYC"
for tf in MAX MYC POL2 #CEPB CREB1 CTCF EGR1 JUND MAX MYC NRF1 SP1 TBP YY1
do
pk="Rep1"
#for pk in Rep1 Rep2
#do
echo "Doing $tf..."
results=$HOME/Dropbox/Projects/results/Chipseq/ENCODE/$lab/$cl/$tf
downloads=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
derived=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Derived
db=$HOME/meme/db/motif_databases/JASPAR_CORE_2009.meme
neg=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/GM78/$tf/Derived
meme-chip -oc $results/$tf-$pk -db $db -nmeme 500 -meme-nmotifs 5 -dreme-m 5 $derived/$tf-$pk.fa 

centrimo --neg $neg/$tf-$pk.fa --verbosity 1 -oc $results/$tf-$pk/GM78-K562 -bgfile $results/$tf-$pk/background -score 5 -ethresh 10 $derived/$tf-$pk.fa $db $results/$tf-$pk/meme_out/meme.xml $results/$tf-$pk/dreme_out/dreme.xml
done
done 
#done
#done
