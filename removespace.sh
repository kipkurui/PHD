#1/bin/bash
#remove spaces from file names

paths="/home/kipkurui/Project/Data/PBM/All_deBruijn/SCI09-pwm4"
cd $paths
ls  > ../Tfs.txt

for tf in $(cat ../Tfs.txt);
	do
		cd $paths/$tf
		rename "s/\s+//g" *
	done
