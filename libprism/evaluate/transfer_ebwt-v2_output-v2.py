#########################################################################
# File Name: transfer_discoSNP_output.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Wed Jul  3 19:44:14 2019
#########################################################################
#!/bin/bash

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import os
import sys
import tools
#from tools import tools
#######################################
#input: .snp 
#output: snp_pair_kmer indel_pair_kmer
#######################################
#vcfFile = "/home/yulin/software/HapCUT2/reproduce_hapcut2_paper/run_hapcut2_fosmid/data/NA12878_hg19_VCFs/chr15.hg19.vcf"


def updataIso(pos, kmer1, kmer2, cov0, cov1):
    if pos+mid+1 > len(kmer1) and pos -mid >= 0:
        temp1 = kmer1[pos-mid, pos+mid+1]
        temp2 = kmer2[pos-mid, pos+mid+1]
        k1, k2 = tools.get_smaller_pair_kmer(temp1, temp2)
        isol.add( (k1, k2) )#, cov0, cov1) )
    return

    
if len(sys.argv) < 3:
    print ("output.snp 6 31")
    sys.exit()

inFile = sys.argv[1]
filterCov = sys.argv[2]
k = int(sys.argv[3])
#inFile = "output.5.snp"
mid = int(k/2)

nonSepSnpFile = "nonSep_snp_"+ str(k) +"mer_cov" + filterCov + "_pair_kmer"
nonIsoSnpFile = "nonIso_snp_"+ str(k) +"mer_cov" + filterCov + "_pair_kmer"
snpFile = "snp_" + str(k) + "mer_cov" + filterCov + "_pair_kmer"
indelFile = "indel_"+ str(k) + "mer_cov"+ filterCov +"_pair_kmer"

print ("input: %s" % (inFile) )
# python3 sorted snp pair kmer correct
# python2 sorted snp pair kmer fail

# main
snpPairKmer, indelPairKmer = {}, {}
state = 0
snpNumber = 0
with open(inFile, "r") as f:
    line = f.readline()
    while line:
        if state == 0:
            assert line.startswith(">")
            words = line.split("_")
            ID = words[0][1:].split(":")[1]
            cov = words[3].split(":")[1]
            #print (ID)
            mutation = words[5]
            if mutation == "INDEL":
                indelLen = len(line.strip().split(":")[-1]) - 1
            state = 1 
            line = f.readline()
            continue
        elif state == 1:
            kmer = line.strip()
            if mutation == "SNP":
                if ID not in snpPairKmer:
                    snpPairKmer[ID] = []
                #snpPairKmer[ID].append(kmer[mid:-mid])
                snpPairKmer[ID].append((kmer, cov))
            elif mutation == "INDEL" and indelLen == 1:
                if ID not in indelPairKmer:
                    indelPairKmer[ID] = []
                indelPairKmer[ID].append((kmer) )
            state = 0
            line = f.readline()
            continue
        else:
            print ("error")

#print snpPairKmer
#print indelPairKmer
print ( "total snp number:", len(snpPairKmer) )
nonSep, nonIsol, isol = set(), set(), set()

for key in snpPairKmer:
    m = snpPairKmer[key]
    if len(m) == 2:
        #assert tools.hamming_distance(m[0], m[1]) == 1
        s0, cov0 = m[0]
        s1, cov1 = m[1]
        if tools.hamming_distance(s0, s1) == 1:
            cLen = len(s0) 
            cmid = int(cLen/2)
            kmer1 = s0[cmid-mid:cmid+mid+1]
            kmer2 = s1[cmid-mid:cmid+mid+1]
            assert len(kmer1) + len(kmer2) == 2*k
            k1, k2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
            isol.add( (k1, k2))#, cov0, cov1) )
        elif tools.hamming_distance(s0, s1) == 2:
            diffPos = tools.hamming_distance2(s0, s1)
            pos1, pos2 = diffPos[0], diffPos[1]
            if pos2-pos1 > mid:
                updataIso(pos1, s0, s1, cov0, cov1)
                updataIso(pos2, s0, s1, cov0, cov1)
            else:
                if pos1-mid>=0 and pos2+mid+1 < len(s0):
                    kmer1 = s0[pos1-mid:pos2+mid+1]
                    kmer2 = s1[pos1-mid:pos2+mid+1] 
                    #assert len(kmer1) + len(kmer2) == 2*k
                    k1, k2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
                    nonIsol.add( (k1, k2))#, cov0, cov1) )

                    h1Pre, h2Pre, h1Suf, h2Suf = kmer1[:k], kmer2[:k], kmer1[-k:], kmer2[-k:]
                    smallerH1P, smallerH2P = tools.get_smaller_pair_kmer(h1Pre, h2Pre)
                    smallerH1S, smallerH2S = tools.get_smaller_pair_kmer(h1Suf, h2Suf)
                    nonSep.add( (smallerH1P, smallerH2P) )
                    nonSep.add( (smallerH1S, smallerH2S) )
    #else:
        #print ("a cluter more than 2 choice", len(m))
        '''
        if len(m) >= 3:
            print (key)
            print ("dis:", tools.hamming_distance(m[0], m[1]) )
            print ("dis:", tools.hamming_distance(m[0], m[2]) )
            print ("dis:", tools.hamming_distance(m[1], m[2]) )
        '''

print ( "isolated snp number:", len(isol) )
sortedSnpPairKmer = sorted(list(isol))
snpOUT = open(snpFile, "w")
#for (kmer1, kmer2, cov0, cov1) in sortedSnpPairKmer:
   #snpOUT.write("%s %s %s %s\n" % (kmer1, kmer2, cov0, cov1) )
for (kmer1, kmer2) in sortedSnpPairKmer:
    snpOUT.write("%s %s\n" % (kmer1, kmer2) )
snpOUT.close()


print ( "non-isolated snp number:", len(nonIsol) )
sortedKmer = sorted(list(nonIsol))
nonOUT = open(nonIsoSnpFile, "w")
for (kmer1, kmer2) in sortedKmer:
    nonOUT.write("%s %s\n" % (kmer1, kmer2) )
nonOUT.close()



print ( "non-isolated sep snp number:", len(nonIsol) )
sortedKmer = sorted(list(nonSep))
nonSepOUT = open(nonSepSnpFile, "w")
for (kmer1, kmer2) in sortedKmer:
    nonSepOUT.write("%s %s\n" % (kmer1, kmer2) )
nonSepOUT.close()


print ( "indel number", len(indelPairKmer) )

indels = set()
countf, countb = 0, 0 
for key in indelPairKmer:
    m = indelPairKmer[key]
    if len(m) == 2:
        #if kmer1[:cmid-1] != kmer2[:cmid-1]:
            #print (kmer1, kmer2, m[0], m[1])
        cLen = len(kmer1)
        cmid = int(cLen/2)
        if m[0][cmid-mid:cmid-1]  == m[1][cmid-mid+1:cmid]: # a pair of indel, m[0] or m[1] which one be deleted
            kmer1 = m[0]
            kmer2 = m[1][1:]
            countf += 1
        elif m[1][cmid-mid:cmid-1]  == m[0][cmid-mid+1:cmid]:
            kmer1 = m[1]
            kmer2 = m[0][1:]
            countb += 1
        else:
            #print (m[0], m[1])
            #break
            continue
        
        temp1 = kmer1[cmid-mid : cmid+mid+1]
        temp2 = kmer2[cmid-mid : cmid+mid]
        check = temp1[:mid+1] + temp1[mid+2:]
        if temp2 == check:
            #print (kmer1, kmer2, temp1, temp2)
            #print (kmer1, kmer2 temp1, temp2)
            assert len(temp1) + len(temp2) == 2*k-1
            small1, small2 = tools.get_smaller_pair_kmer(temp1, temp2)
            indels.add( (small1, small2) )
         
        initialH1 = temp1
        i=1
        while mid+i<k and initialH1[mid+i] == initialH1[mid] and cmid+mid+1+i < cLen:
            temp1 = kmer1[cmid-mid+i : cmid+mid+1+i]
            temp2 = kmer2[cmid-mid+i : cmid+mid+i]
            check = temp1[:mid+1] + temp1[mid+2:]
            if check == temp2:
                assert len(temp1) + len(temp2) == 2*k-1
                small1, small2 = tools.get_smaller_pair_kmer(temp1, temp2)
                indels.add( (small1, small2) )
            i+=1
        i=1
        while mid-i>=0 and initialH1[mid-i] == initialH1[mid] and cmid+mid+1+i < cLen:
            temp1 = kmer1[cmid-mid-i : cmid+mid+1-i]
            temp2 = kmer2[cmid-mid-i : cmid+mid-i]
            check = temp1[:mid+1] + temp1[mid+2:]
            if check == temp2:
                assert len(temp1) + len(temp2) == 2*k-1
                small1, small2 = tools.get_smaller_pair_kmer(temp1, temp2)
                indels.add( (small1, small2) )
            i+=1
        
    #else:
        #print ("a cluter more than 2 choice, indel", len(m))

print ("delete in second sequence", countb, "delete in first sequence", countf)
print ( "indel pair number", len(indels) )
sortedIndelPairKmer = sorted(list(indels))
indelOUT = open(indelFile, "w")
for (kmer1, kmer2) in sortedIndelPairKmer:
    indelOUT.write("%s %s\n" % (kmer1, kmer2) )
indelOUT.close()


