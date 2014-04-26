#Auestion now is how to perform this on multiple files.
#read the text file
files <- list.files(path="/home/kipkurui/PHD-project/Data/k-mers/All_Contig8mers/SCI09-counts", pattern=".txt", all.files=T, full.names=T)
for (file in files) {
  kmers = read.table(file)
  #create an extra column after calculations
  kmers[5]=(kmers[4]/sum(kmers[4]))*kmers[3]
  #Write the output to a new file 
  l=nchar(file)
  d=69
  kmername=substr(file,d,l)
  #trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  #kmername <- trim(kmername)
  write.table(kmers, file=paste("/home/kipkurui/PHD-project/Data/k-mers/All_Contig8mers/SC09-freq/",kmername),row.names=F,col.names=F, sep="\t", quote=FALSE)
  #write.table(kmers,file="output.txt",row.names=F,col.names=F,sep="\t", quote = FALSE)
}