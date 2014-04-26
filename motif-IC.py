#!/usr/bin/python
from __future__ import print_function
import sys
from math import log

found = 0;
row = 0
nrows = 0
entropy = 0
totalentropy = 0
motifs = 0
motiffile=raw_input("Enter motif file: ")
outfile=raw_input("Enter motif file: ")
outfile=open(outfile, "w")
motiffile=open(motiffile, "r")
for line in motiffile: #sys.stdin.readlines():
  words = line.split()
  if found == 0 :
    if line.startswith ("MOTIF"):
      # allow for motifs without an alternative name
      if len(words) < 3:
        words.append("")
      outfile.write("%s \t %s \t " % (words[1],words[2]))
      found = 1
      motifs = motifs + 1
      entropy = 0
      continue
  if found == 1:
    if line.startswith ("letter-probability"):
      nrows = int((line.split("w="))[1].split()[0])
      found = 2
    continue
  if found == 2:
    #outfile.write("< %s > %s \n" % (row,words))
    for val in words:
        if float(val) > 0:
           entropy += float(val) * log(float(val))/log(2.0)
    row += 1
    if row >= nrows:
      v=2*nrows+entropy
      outfile.write(str(v))
      outfile.write("\n")
      found = 0
      row = 0
      totalentropy = totalentropy + 2*nrows+entropy

if row > 0:
  outfile.write(str(v))
  outfile.write("\n")
  totalentropy = totalentropy + 2*nrows+entropy
    
if motifs > 0:
  outfile.write("%s motifs mean IC= %s \n \n" % (motifs,totalentropy/motifs))
outfile.close()
motiffile.close()
