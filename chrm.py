import os
labs=('HAIB','SYDH','OpenChromUTA')
cells=('GM78','hESC','K562')
tfs=('ATF3','CTCF','MYC','CREB1','EGR1','YY1','SP1','MAX','JUND','NRF1','TBP','CEPB')
reps=('Rep1','Rep2')
writefile=open('chrM.txt','w')
chrom=''
for lab in labs:
    for cl in cells:
        for tf in tfs:
            for rep in reps:
                path=("/home/Caleb/Dropbox/Projects/Data/Chipseq/ENCODE/%s/%s/%s/Derived/%s-%s.bed"%(lab,cl,tf,tf,rep))
                if os.path.exists(path):
                    bedfile=open(path,'r')
                    print bedfile
                    for line in bedfile:
                        if 'chrM' in line:
                            w=("%s/%s/%s/Derived/%s-%s.bed \n"%(lab,cl,tf,tf,rep))
                            chrom+=w
                            print w
                            #writefile.writelines(w)
                        else:
                            continue
                else:
                    continue
for i in chrom:
    writefile.write(i)
writefile.close()
