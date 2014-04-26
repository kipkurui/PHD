import os, sys
import subprocess

os.chdir("/home/kipkurui/PHD-project/Data/k-mers/All_Contig8mers")

#create directory to put the seeds
os.system("mkdir SC09-seeds")
os.system("ls ./SC09-freq >kmers.txt")
files=open("kmers.txt","r")
#os.system("rename "s/ *//g" *.txt")

#extract the top 4 and bottom 2 as the seeds and save using the seed name
for line in files:
    #print line.strip()
    line=line.strip()
    os.system("sort -r -g -k5 ./SC09-freq/%s | tail -4 | cut -f 1,3,4 | sed 's/\t/-/g' >./SC09-seeds/%s "% (line,line))
    os.system("sort -r -g -k5 ./SC09-freq/%s | head -2 | cut -f 1,3,4 | sed 's/\t/-/g' >>./SC09-seeds/%s "% (line,line))
path="/home/kipkurui/PHD-project/Data/k-mers/All_Contig8mers/SC09-seeds"
os.system("ls ./SC09-seeds >kmers.txt")
#files=open("kmers.txt","r")

#Need to put more effort in naming variables

s1=open("./kmers.txt", "r")
#os.system("ls /home/kipkurui/Project/Data/PBM/All_deBruijn/SCI09 >debru.txt")
s=open("/home/kipkurui/Project/Data/PBM/All_deBruijn/debru.txt", "r")
k=s.readlines()
for tfs in s1:
    os.chdir("/home/kipkurui/PHD-project/Data/k-mers/All_Contig8mers")
    tf=tfs.strip()[:-16]
    for n in range(1,6):
        if tf+".%s_v1_deBruijn.txt\n" % (str(n)) in k:
            tfs=tfs.strip()
            files=open("/home/kipkurui/PHD-project/Data/k-mers/All_Contig8mers/SC09-seeds/%s" %(tfs),"r")
            os.chdir("/home/kipkurui/Project/Data/PBM/All_deBruijn")
            for seed in files:
                os.system("mkdir -p ./SCI09-PWM/%s" % (tf[:-5]))
                #print tf[:-5]
                errf=open("./log.txt", "a")
                with open ("./SCI09-PWM/%s/%s-%s.uniprobe" % (tfs.strip()[:-21],tfs.strip()[:-21]+"-"+seed.strip()[:14]+seed.strip()[26:],str(n)), "w") as outfile:
                    subprocess.call("perl ./wobble_single_seed_twoarray.pl ./SCI09/%s.%s_v1_deBruijn.txt ./SCI09/%s.%s_v2_deBruijn.txt %s ./patterns_4x44k_all_8mer.txt" \
                                   % (tf,str(n),tf,str(n),seed[:8]), shell=True, stdout=outfile, stderr=errf)
        else:
            continue
files.close()
