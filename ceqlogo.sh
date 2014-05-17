#!/bin/bash

save=$HOME/PHD-project/Results/logos
db=$HOME/meme/db/motif_databases/uniprobe_mouse.meme

#for i in {1..`grep -c "MOTIF" uniprobe_mouse.meme`};
grep MOTIF $db | cut -d " " -f 3 >$save/names.txt

for MOTIF in $(cat  $save/names.txt);
do
	echo "Now on $MOTIF"
	#ceqlogo -i1 /home/kipkurui/Project/Results/PBM/SCI09-PWM/Max/Max1_pwm.meme -Y -N -y IC -t Max-1 > Max-1.eps
	#ceqlogo -i`grep -n $MOTIF | cut -f 1 -d ":"` $db/ -t "$MOTIF" -b 2 -h 6 -o $save/$MOTIF.eps
	ceqlogo `grep -n $MOTIF /home/kipkurui/meme/db/motif_databases/uniprobe_mouse.meme |cut -f 1 -d ":"` $db -t "$MOTIF" -b 2 -h 6 -o $save/$MOTIF.eps

done

/home/kipkurui/meme/db/motif_databases/uniprobe_mouse.meme
