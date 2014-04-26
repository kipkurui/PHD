#1/bin/bash
#remove spaces from file names


paths="/home/kipkurui/Project/Results/PBM/PWM/SCI09"
results="/home/kipkurui/Project/Results/PBM/PWM"

cd $paths
ls  > /tmp/Tfs.txt

for tf in $(cat /home/kipkurui/Desktop/tennew);
do
	cd $paths/$tf
	#ls > /tmp/mots.txt
	#for motif in $(cat /tmp/mots.txt);
	#do	
		#touch $results/PWM/$tf.txt
		ls >>/home/kipkurui/Desktop/mymotifs.txt
		#cat *.uniprobe >$results/PWM/$tf.txt
		#cat $results/All-seednwobble-pwm-2.uniprobe $motif >$results/hold.txt
		#cat $results/hold.txt $results/final.txt >>$results/All-seednwobble-pwm-2.uniprobe
	done
#done

#uniprobe2meme -logodds  < calall.uniprobe >! Calnew9.meme
