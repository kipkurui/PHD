#!/bin/bash
#Caleb Kibet
#Make the directoried where data and results are stored
file1=$HOME/Desktop/
for lab in Sydh Haib
do 
	for cl in Gm12878 H1hesc Hepg2 K562 Helas3 T47d Shsy5y Imr90 A549 Nb4 Hct116 Hek293 Panc1 Mcf7
	do
		for tf in $(cat  $file1/newtfs);
		do
			mkdir -p $HOME/Project/Results/Chipseq/ENCODE/$lab/$cl/$tf
			mkdir -p $HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
			mkdir -p $HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Derived
		done
	done
done
