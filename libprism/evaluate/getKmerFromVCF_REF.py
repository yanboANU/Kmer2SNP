#########################################################################
# File Name: getKmerFromVCF_REF.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 10:45:06 AEST
#########################################################################
#!/bin/bash

import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys
#from libprismv2.local import tools
import tools
#from libprism.evaluate import read
import read
#import read

def write_pair_kmer(outFile, kmers):

    sortedKmers = sorted(kmers)
    with open(outFile, "w") as f:
        for eles in sortedKmers:
            #f.write("%s %s %s %s %s %s %s %s\n" % (ele[0], ele[1], ele[2], ele[3], ele[4], ele[5], tools.reverse(ele[0]), tools.reverse(ele[1]) ) )
            for e in eles:
                f.write("%s " % (e) )
            f.write("\n")    

def get_snp_pair_kmer(vcfFilename):

    snps = read.read_vcf(vcfFilename)
    sortedSNP = sorted(snps)
    kmerFilename="chr" + sys.argv[1] + ".snp.real." + sys.argv[2] + "mer"
    nonFilename="chr" + sys.argv[1] + ".nonIsolated.snp.real." + sys.argv[2] + "mer"
    # two snp close to each other, use one kmer pair (put one snp in the middle) represent
    nonSepFilename="chr" + sys.argv[1] + ".non.sep.snp.real." + sys.argv[2] + "mer"
    kmers = []
    nons = []
    nonSeps=[]
    #for key in snps:  
    mid = int(k/2)
    i = 0
    while i < len(sortedSNP):
        key = sortedSNP[i]
        if i+1 < len(sortedSNP):
            after = sortedSNP[i+1]
        if i==len(sortedSNP)-1 or  after-key > mid:
            assert seq[key-1] == snps[key][0] or  seq[key-1] == snps[key][1]
            h1 = seq[key-mid-1 : key-1] + snps[key][0] + seq[key : key+mid ] # 0
            h2 = seq[key-mid-1 : key-1] + snps[key][1] + seq[key : key+mid ] # 1
            if h1.count('N') > 0 or h2.count('N') > 0:
                continue
            smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
            kmers.append( (smallerH1, smallerH2, key) )
        else:
            #assert seq[key-1] == snps[key][0] and seq[after-1] == snps[after][0]
            assert seq[key-1] == snps[key][0] or  seq[key-1] == snps[key][1]
            assert seq[after-1] == snps[after][0] or  seq[after-1] == snps[after][1]
            h1 = seq[key-mid-1 : key-1] + snps[key][0] + seq[key : after-1 ] + snps[after][0] + seq[after : after+mid ]# 0
            h2 = seq[key-mid-1 : key-1] + snps[key][1] + seq[key : after-1 ] + snps[after][1] + seq[after : after+mid ]  # 1
            if h1.count('N') > 0 or h2.count('N') > 0:
                continue
            smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
            nons.append( (smallerH1, smallerH2, key, after) )
            
            h1Pre, h2Pre, h1Suf, h2Suf = h1[:k], h2[:k], h1[-k:], h2[-k:]
            smallerH1P, smallerH2P = tools.get_smaller_pair_kmer(h1Pre, h2Pre)
            smallerH1S, smallerH2S = tools.get_smaller_pair_kmer(h1Suf, h2Suf)
            #if smallerH1P <= smallerH1S:
            nonSeps.append( (smallerH1P, smallerH2P, key, after) )
            #else:
            nonSeps.append( (smallerH1S, smallerH2S, key, after) )


            i += 1
        i += 1    

    write_pair_kmer(kmerFilename, kmers)
    write_pair_kmer(nonFilename, nons)
    write_pair_kmer(nonSepFilename, nonSeps)


def get_indel_pair_kmer(vcfFilename):

    indels = read.read_vcf(vcfFilename)
    kmerFilename="chr" + sys.argv[1] + ".indel.real." + sys.argv[2] + "mer"
    kmers = []
    indel_length1_cnt = 0
    for key in indels:
        s1, s2, ID = indels[key]
        lenS1, lenS2 = len(s1), len(s2)
        if lenS1 + lenS2 > 3:
            continue

        indel_length1_cnt += 1
        assert lenS1 + lenS2 >= 2
        if len(s1) == 1 and len(s2) == 2:
            assert seq[key-1] == s1
            if s2[0] == s2[1]:
                print (seq[key-1], s1, s2)
                continue
            assert s2[0] != s2[1]
            #while s2[1] == seq[key-1]: # delete content is s2[1]
                #key = key-1            # delete happen at "AAA" region, always think delete first poisition                  
            h1 = seq[key-int(k/2) : key+int(k/2)] # k-1 
            h2 = seq[key-int(k/2) : key-1] + s2 + seq[key : key+int(k/2)] # 1 # len: k
            assert len(h1) == k-1 and len(h2) == k
            h1, h2 = h2, h1 # h1 always is longer one
            initialH1 = h1
            if h1.count('N') > 0 or h2.count('N') > 0:
                continue
            #print key, "11", h1
            smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
            kmers.append( (smallerH1, smallerH2, key) )
            # delete happen at multipe "AAAA" region, more pair kmer happen
            l = len(h1) 
            mid = int(l/2)
            i=1
            while mid+i<l and initialH1[mid+i] == initialH1[mid]:
                h1 = seq[key-int(k/2)+i : key+int(k/2)+i] # move right i
                h2 = seq[key-int(k/2)+i : key-1] + s2 + seq[key : key+int(k/2)+i] # move right i
                h1, h2 = h2, h1 # h1 always is longer one
                smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
                #print key, "aa"
                kmers.append( (smallerH1, smallerH2, key) )
                i+=1
            i=1
            while mid-i>=0 and initialH1[mid-i] == initialH1[mid]:
                h1 = seq[key-int(k/2)-i : key+int(k/2)-i] # move left i
                h2 = seq[key-int(k/2)-i : key-1] + s2 + seq[key : key+int(k/2)-i] # move right i
                h1, h2 = h2, h1 # h1 always is longer one
                smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
                #print key, "bb"
                kmers.append( (smallerH1, smallerH2, key) )
                i+=1

            ''' # for test can grouth-truth can always keep min strand delete first
            if h1 > tools.reverse(h1):
                print "aa"
                print h1, h2
                print tools.reverse(h1), tools.reverse(h2)
                while s2[1] == seq[key]:
                    key+=1
                h1 = seq[key-int(k/2) : key+int(k/2)] # k-1 
                h2 = seq[key-int(k/2) : key] + s2[1] + seq[key : key+int(k/2)] # 1 # len: k
                h1, h2 = h2, h1
                print h1, h2
                print tools.reverse(h1), tools.reverse(h2)
            '''    
        elif len(s1) == 2 and len(s2) == 1: 
            if seq[key-1:key+1] != s1 or s1[0] == s1[1]:
                print (seq[key-1:key+1], s1, s2)
                continue
            assert seq[key-1:key+1] == s1
            assert s1[0] != s1[1]            
            h1 = seq[key-int(k/2) : key+int(k/2)+1] # k 
            h2 = seq[key-int(k/2) : key] + seq[key+1 : key+int(k/2)+1] # 1 # len: k-1
            assert len(h1) == k and len(h2) == k-1
            initialH1 = h1
            if h1.count('N') > 0 or h2.count('N') > 0:
                continue
            smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
        
            #print key, "22"
            kmers.append( (smallerH1, smallerH2, key) )

            l = len(h1) 
            mid = int(l/2)
            i=1
            while mid+i<l and initialH1[mid+i] == initialH1[mid]:
                h1 = seq[key-int(k/2)+i : key+int(k/2)+1+i] # k 
                h2 = seq[key-int(k/2)+i : key] + seq[key+1 : key+int(k/2)+1+i] # 1 # len: k-1
                smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
                #print key, "cc"
                kmers.append( (smallerH1, smallerH2, key) )
                i+=1
            i=1
            while mid-i>=0 and initialH1[mid-i] == initialH1[mid]:
                h1 = seq[key-int(k/2)-i : key+int(k/2)+1-i] # k 
                h2 = seq[key-int(k/2)-i : key] + seq[key+1 : key+int(k/2)+1-i] # 1 # len: k-1
                smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
                #print key, "dd"
                kmers.append( (smallerH1, smallerH2, key) )
                i+=1
    print ("there are ", indel_length1_cnt, "indels, create ", len(kmers), "indel pair kmer")

    write_pair_kmer(kmerFilename, kmers)


'''
allFile = "chr" + sys.argv[1] + ".all." + sys.argv[2] + "mer"
foutAll = open(allFile, "w")
for i in range(seqLen-21):
    mer = seq[i:i+k]
    if mer.count('N') > 0:
        continue
    Rmer = tools.reverse(mer)
    if Rmer < mer:
        mer = Rmer
    foutAll.write("%s %s\n" % (mer, i))    
foutAll.close()       
'''


# this simulate data is based on hg18

if len(sys.argv) < 2:
    print ("chrID kmer size")
    sys.exit()

#path=/home/yulin/bio/VariationCalling/data/
#path=/home/yanbo/bio/Kmer2SNP/data/HG002

#NA12878
#refFilename="/home/yulin/bio/Data/reference/NCBI36_hg18/chr"+ sys.argv[1] +".fa"
#snpVCFFile="/home/yulin/bio/VariationCalling/data/NA12878/VCF/NA12878_hg18_snp_VCFs/chr" +sys.argv[1] + ".vcf"
#indelVCFFile="/home/yulin/bio/VariationCalling/data/NA12878/VCF/NA12878_hg18_indel_VCFs/chr"+ sys.argv[1] +".vcf"
       
#HG002
#refFilename="/home/yulin/bio/Data/reference/GRCh38_hg38/chr"+ sys.argv[1] +".fa"
refFilename="/home/yanbo/bio/Kmer2SNP/data/HG002/GRCh38_reference/chr"+ sys.argv[1] +".fa"
snpVCFFile="/home/yanbo/bio/Kmer2SNP/data/HG002/HG002_GRCh38_snp_VCFs/chr" +sys.argv[1] + ".vcf"
indelVCFFile="/home/yanbo/bio/Kmer2SNP/data/HG002/HG002_GRCh38_indel_VCFs/chr"+ sys.argv[1] +".vcf"
        

record = SeqIO.read(open(refFilename), "fasta")
print (record.id)
seq = str(record.seq).upper()
seqLen = len(seq)
k=int(sys.argv[2])
print(snpVCFFile)
get_snp_pair_kmer(snpVCFFile)
print(indelVCFFile)
get_indel_pair_kmer(indelVCFFile)

