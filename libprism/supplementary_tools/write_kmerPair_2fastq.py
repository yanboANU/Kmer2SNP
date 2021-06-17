#########################################################################
# File Name: write_kmerPair_2fastq.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Mon 07 Jun 2021 09:06:48 PM AEST
#########################################################################
#!/bin/bash
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from pbcore.io import fastqWriter

if __name__ == '__main__':   

    if len(sys.argv) < 2:
        print ("k_31*snp kmerPair.fastq/extendKmerPair.fastq")
        sys.exit()
    cnt = 0
    records = []
    with open(sys.argv[1], "r") as f:
        for line in f:
            words = line.strip().split()
            cnt += 1
            r0 = SeqRecord( Seq(words[0]) , id= str(cnt) +"_0" )
            r0.letter_annotations["phred_quality"] = [10]*len(words[0])
            r1 = SeqRecord( Seq(words[1]) , id= str(cnt) +"_1" )
            r1.letter_annotations["phred_quality"] = [10]*len(words[1])
            records.append(r0)
            records.append(r1)
    #SeqIO.write(records, "kmerPair.fastq", "fastq")
    SeqIO.write(records, sys.argv[2], "fastq")
