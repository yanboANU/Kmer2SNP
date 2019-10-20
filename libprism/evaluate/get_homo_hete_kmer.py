#########################################################################
# File Name: get_homo_hete_kmer.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 01 Oct 2019 10:47:49 AEST
#########################################################################
#!/bin/bash

import os
import sys
#dsk first
#/home/yulin/software/dsk/build/bin/dsk -nb-cores 10 -file chr22_2_haplotypes.fa -out chr22_k31 -kmer-size 31 -abundance-min 1
#/home/yulin/software/dsk/build/bin/dsk2ascii -nb-cores 10 -file chr22_k31 -out chr22_k31.txt


m={}
m['A'] = 'T'
m['T'] = 'A'
m['C'] = 'G'
m['G'] = 'C'

def reverse(s): # AATC
    news = ""
    for ele in s:
       news = m[ele] + news #CGATT
    return news   


homoFile = open("chr22_homo_k31.txt", "w")
heteFile = open("chr22_hete_k31.txt", "w")
repeatFile = open("chr22_repeat_k31.txt", "w")
with open("chr22_k31.txt") as f:
    for line in f:
        words = line.split(' ')
        cov = int(words[1])
        kmer = words[0]
        Rkmer = reverse(kmer)
        if kmer > Rkmer:
            kmer = Rkmer
        if cov == 2:
            homoFile.write("%s\n" % kmer)
        elif cov == 1:
            heteFile.write("%s\n" % kmer)
        elif cov > 2:
            repeatFile.write("%s\n" % kmer)
homoFile.close()
heteFile.close()
repeatFile.close()
#sort three files
