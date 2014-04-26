#!/bin/bash
# Author: Caleb Kibet
lab="HAIB"
#for lab in  SYDH
#do
#for cl in K562
#do
cl="hESC"
ne="GM78"
#for ne in GM78
#do
#tf="MYC"
for tf in ATF3 TBP MAX POL2 CEPB CREB1 CTCF EGR1 JUND NRF1 SP1 TBP YY1 MYC ATF3
do

echo "Doing $tf..."
pk="Rep1"
#for pk in Rep1
#do
results=$HOME/Dropbox/Projects/results/Chipseq/ENCODE/$lab/$cl/$tf
downloads=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
derived=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$cl/$tf/Derived
db=$HOME/meme/db/motif_databases/JASPAR_CORE_2009.meme
neg=$HOME/Dropbox/Projects/Data/Chipseq/ENCODE/$lab/$ne/$tf/Derived

centrimo --neg $neg/$tf-$pk.fa --verbosity 1 -oc $results/$tf-$pk/hESC-GM78 -bgfile $results/$tf-$pk/background -score 5 -ethresh 10 $derived/$tf-$pk.fa $db $results/$tf-$pk/meme_out/meme.xml $results/$tf-$pk/dreme_out/dreme.xml
#done
done 
#done
#done

