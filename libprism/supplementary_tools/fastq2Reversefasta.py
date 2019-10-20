#########################################################################
# File Name: fastq2Reversefasta.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Mon 29 Jul 2019 17:02:56 AEST
#########################################################################
#!/bin/bash
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tools
inFile = sys.argv[1]
outFile = sys.argv[2]
#records = []
nosiy = ['W', 'N', 'R', 'M', 'K', 'Y', 'S']
fout = open(outFile, "w")
for record in SeqIO.parse(inFile, "fastq"):
    #print(record.id)
    s = 0
    for ele in nosiy:
        s += record.seq.count(ele)
    if s > 0:
        continue
    rec1 = SeqRecord(record.seq, id=record.id+"_1")
    rec2 = SeqRecord(tools.reverse(record.seq), id=record.id+"_2")
    #records.append(rec1)
    #records.append(rec2)
    fout.write(">%s\n" % rec1.id)
    fout.write("%s\n" % rec1.seq)
    fout.write(">%s\n" % rec2.id)
    fout.write("%s\n" % rec2.seq)
#SeqIO.write(records, sys.argv[2], "fasta")
fout.close()
