import os, sys
import subprocess

#Need to put more effort in naming variables

#s1=open("./kmers.txt", "r")
#DNase="/home/kipkurui/Project/Data/DNase/DnaseClusteredV2.fa"
DNase="/home/kipkurui/Project/Data/Chipseq/ENCODE/Derived/RegTfbsClusteredV3.fa"
patterns="/home/kipkurui/Project/Data/PBM/Patterns"
debruin="/home/kipkurui/Project/Data/PBM/All_deBruijn"
results="/home/kipkurui/Project/Results/PBM"
script="/home/kipkurui/Project/Scripts/PHD"
	
os.chdir(patterns) #get the patterns

os.system("ls /home/kipkurui/Project/Data/PBM/All_deBruijn/SCI09 >SCI09_deBruijn_list.txt") #save the cimbinatorial file

f1=open("./SCI09_deBruijn_list.txt", "r")
debru=f1.readlines()
tflist=open("./SCI09-tfs.txt", "r")
for tfs in tflist:
    #os.chdir("/home/kipkurui/PHD-project/Data/k-mers/All_Contig8mers")
    tf=tfs.strip()[:-16]
    for n in range(1,6):
        if tf+".%s_v1_deBruijn.txt\n" % (str(n)) in debru:
            tfs=tfs.strip()
            os.system("mkdir -p %s/SCI09-PWM/%s" % (results,tf[:-5]))
            print "Now on: %s"% tfs.strip()[:-21]+"-"+str(n)
            os.chdir(patterns)
            errf=open("./log.txt", "a")
            os.chdir("%s/SCI09-PWM/%s" % (results,tf[:-5]))
            #with open ("%s/SCI09-PWM/%s/%s-%s.uniprobe" % (results,tfs.strip()[:-21],tf,str(n)), "w") as outfile:
            subprocess.call("perl %s/cal-seed_and_wobble_twoarray.pl %s/SCI09/%s.%s_v1_deBruijn.txt %s/SCI09/%s.%s_v2_deBruijn.txt 8 %s/patterns_8of10.txt %s/patterns_4x44k_all_8mer.txt %s %s" 
                                   % (patterns,debruin,tf,str(n),debruin,tf,str(n),patterns,patterns,tfs.strip()[:-21]+str(n),DNase), shell=True, stderr=errf)
        else:
            continue
f1.close()
tflist.close()
errf.close()

#perl $script/cal-seed_and_wobble_twoarray.pl $debru/$kmer*v1_deBruijn.txt $debru/$kmer*v2_deBruijn.txt 8 $script/patterns_8of10.txt $script/patterns_4x44k_all_8mer.txt $kmer $
