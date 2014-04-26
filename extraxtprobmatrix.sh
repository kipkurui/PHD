#1/bin/bash
#remove spaces from file names


paths="/home/kipkurui/Project/Data/PBM/All_deBruijn/SCI09-PWM"

cd $paths
ls  > /tmp/Tfs.txt
for tf in $(cat /tmp/Tfs.txt);
do
	cd $paths/$tf
	ls > /tmp/mots.txt
	for motif in $(cat /tmp/mots.txt);
	do
		mkdir -p "/home/kipkurui/Project/Results/PBM/PWM/SCI09/$tf"
		results="/home/kipkurui/Project/Results/PBM/PWM/SCI09/$tf"
		echo $motif | sed 's/.uniprobe//g' > $results/$motif
		tail -6 $motif >> $results/$motif
		#rename "s/\s+//g" *
	done
done
