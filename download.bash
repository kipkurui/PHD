#!/bin/bash
#Caleb Kibet
#05-05-2013
#DownloadChipseq data
file1=$HOME/Desktop/
for lab in Sydh Haib
do 
	for cl in Gm12878 H1hesc Hepg2 K562 Helas3 T47d Shsy5y Imr90 A549 Nb4 Hct116 Hek293 Panc1 Mcf7
	do
		for tf in $(cat  $file1/newtfs);
		do
			echo "downloading $tf ..."
				#uncomment below if directories had not been created
				#mkdir -p $HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
				cd $HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
				wget -nd  -nc -a,  --output-file=$HOME/Project/Data/Chipseq/ENCODE/downlog3.txt http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbs"$lab""$cl""$tf"UniPk.narrowPeak.gz
		done
	done
done
for lab in #enter the labs eg Haib Uta
do 
		for cl in # cell lines eg Gm12878 K562
		do
			for tf in #list of TFs to download
			do
				echo "downloading $tf ..."
				#uncomment below if directories had not been created
				#mkdir -p $HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
				cd $HOME/Project/Data/Chipseq/ENCODE/$lab/$cl/$tf/Downloads
				wget -nd  -nc -a,  --output-file=$HOME/Project/Data/Chipseq/ENCODE/downlog3.txt http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbs"$lab""$cl""$tf"UniPk.narrowPeak.gz
		done 
	done
done
