f = open("TFS.txt", "r")
tf=["Atf3","Ctcf","Cmyc","Creb1","Egr1","Yy1","Sp1","Max","Jund","Nrf1","Tbp","Cepb"]
Cline= ["K562","Gm12878","H1hesc"]
Prot = ["Iggmus","Iggrab","V0", "Std", "Pcr1x"]
Rep = ["Rep1", "Rep2"]
typ = ["broadPeak","narrowPeak"]
writefile = open("Output.txt", "a")
writefile.write("Group \t Lab \t Cell-line \t Tf \t Protocol \t Replication \t File-Type \n")

for line in f:
    rep = "single"
    group = line[2:8]
    lab = line[8:12]
    for i in Cline:
        if i in line:
            cline=i
    for j in tf:
        if j in line:
            Tf = j
    for k in Prot:
        if k in line:
            prot = k
    for l in Rep:
        if l in line:
            rep = l  
    for m in typ:
        if m in line:
            Typ = m
    writefile.write("%s \t %s \t %s \t %s \t %s \t %s \t %s \n" % (group, lab, cline, Tf, prot, rep, Typ))
writefile.close()
f.close()

    

    
