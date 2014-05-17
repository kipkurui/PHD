import os, sys
import subprocess
import os.path
subprocess.call("export LC_ALL=C", shell=True)
os.chdir("/home/kipkurui/Project/Results/PBM/SCI09-PWM")
#os.system("ls >/home/kipkurui/Project/Results/PBM/SCI09.txt")
tflist=open("/home/kipkurui/Project/Results/PBM/SCI09.txt", "r")
for tfs in tflist:
    tf=tfs.strip()
    for n in range(1,6):
        path="/home/kipkurui/Project/Results/PBM/SCI09-PWM/%s/%s%s_8mers_pwm_combined.txt" %(tf,tf,n)
        #print path
        if os.path.isfile(path):
            print n
            path2="/home/kipkurui/Project/Results/PBM/SCI09-PWM/%s/%s%s_pwm.uniprobe" % (tf,tf,n)
            with open(path) as f1:
                with open(path2, 'w') as f2:
                    Lines=f1.readlines()
                    for i, line in enumerate(Lines):
                            print i
                            if line.startswith("1") or line.startswith("2") or line.startswith("3"):
                                    f2.write("%s-%s\n" % (tf,line))
                            if line.startswith("Probability"):
                                    f2.write(Lines[i+2])
                                    f2.write(Lines[i+3])
                                    f2.write(Lines[i+4])
                                    f2.write(Lines[i+5]+"\n")
            res="/home/kipkurui/Project/Results/PBM/SCI09-PWM/%s/%s%s_pwm.uniprobe >/home/kipkurui/Project/Results/PBM/SCI09-PWM/%s/%s%s_pwm.meme" %(tf,tf,n,tf,tf,n)
            #print res
            err=open("../log.txt","a")
            subprocess.call("uniprobe2meme -logodds %s" % res,shell=True,stderr=err)
        else:
            #os.system("rm /home/kipkurui/Project/Results/PBM/SCI09-PWM/%s/%s%s_pwm.meme" %(tf,tf,n))
            continue
tflist.close()
err.close()
