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
#######################################
#input: .vcf and .fa
#output: snp_pair_kmer indel_pair_kmer
#######################################

if len(sys.argv) < 3:
    print ("python3 *vcf *fa k=31")
    sys.exit()


vcfFile = sys.argv[1]
faFile = sys.argv[2] 
k = int(sys.argv[3])
mid = int(k/2)

nonSnpFile = "non_snp_" + str(k) + "mer_pair_kmer"
snpFile = "snp_" + str(k) + "mer_pair_kmer"
indelFile = "indel_" + str(k)  + "mer_pair_kmer"

print ("input: %s %s %s" % (vcfFile, faFile, k) )

def read_vcf(filename):
    mutations = {}
    with open(filename, "r")as f:
        for line in f:
            if line.startswith('#'):
                continue
            words = line.split()
            contigID = words[0] 
            pos = int(words[1])
            s1 = words[3].strip()
            s2 = words[4].strip()
            homo = words[9].split(':')[0]
            #print homo
            if homo == "0/0" or homo == "0|0" or homo == "1|1" or homo == "1/1":
                continue
            assert homo == "1/0" or homo == "1|0" or homo == "0|1" or homo == "0/1"
            s1Len = len(s1)
            s2Len = len(s2)
            if contigID not in mutations:
                mutations[contigID] =[]
            mutations[contigID].append( ( pos, s1, s2, homo) )
    return mutations    

mutations = read_vcf(vcfFile)
print ( "mutations number:", len(mutations) )

def updata_one_heterozygous_pair(s1, s2, pos, seq, snpPairKmer, indelPairKmer):
    
    if len(s1) + len(s2) >3:
        return
    if len(s1) == 1 and len(s2) == 1:
        if seq[pos] != s1.upper():
            print ("snp error pos range", ID, pos, seq[pos-2:pos+3], s1, s2)
            print (seq)
            sys.exit()
        assert seq[pos] == s1.upper()
        kmer1 = seq[pos-mid:pos+mid+1]
        kmer2 = seq[pos-mid:pos] + s2.upper() + seq[pos+1:pos+mid+1]
        assert len(kmer1) + len(kmer2) == 2*k
        smallerKmer1, smallerKmer2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
        snpPairKmer.add( (smallerKmer1, smallerKmer2) )
    elif len(s1) == 1 and len(s2) == 2:
        kmer1 = seq[pos-(mid-1):pos+(mid+1)]
        kmer2 = seq[pos-(mid-1):pos] + s2.upper() + seq[pos+1:pos+(mid+1)]
        assert len(kmer1) + len(kmer2) == 2*k - 1
        smallerKmer1, smallerKmer2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
        indelPairKmer.add( (smallerKmer1, smallerKmer2) )
        initialH1 = kmer2
        i=1
        while mid+i<k and initialH1[mid+i] == initialH1[mid]:
            kmer1 = seq[pos-(mid-1)+i:pos+(mid+1)+i]
            kmer2 = seq[pos-(mid-1)+i:pos] + s2.upper() + seq[pos+1:pos+(mid+1)+i]
            assert len(kmer1) + len(kmer2) == 2*k - 1
            smallerH1, smallerH2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
            indelPairKmer.add( (smallerH1, smallerH2) )
            i+=1
        i=1
        while mid-i>=0 and initialH1[mid-i] == initialH1[mid]:
            kmer1 = seq[pos-(mid-1)-i:pos+(mid+1)-i]
            kmer2 = seq[pos-(mid-1)-i:pos] + s2.upper() + seq[pos+1:pos+(mid+1)-i]
            assert len(kmer1) + len(kmer2) == 2*k - 1
            smallerH1, smallerH2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
            indelPairKmer.add( (smallerH1, smallerH2) )
            i+=1

    elif len(s1) == 2 and len(s2) == 1:   # vcf pos record correct but s1 and s2 not correct
        kmer1 = seq[pos-(mid-1):pos+(mid+2)]
        kmer2 = seq[pos-(mid-1):pos+1] + seq[pos+2:pos+(mid+2)]
        assert len(kmer1) + len(kmer2) == 2*k - 1
        smallerKmer1, smallerKmer2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
        indelPairKmer.add( (smallerKmer1, smallerKmer2) )
        initialH1 = kmer1
        i=1
        while mid+i<k and initialH1[mid+i] == initialH1[mid]:
            kmer1 = seq[pos-(mid-1)+i:pos+(mid+2)+i]
            kmer2 = seq[pos-(mid-1)+i:pos+1] + seq[pos+2:pos+(mid+2)+i]
            assert len(kmer1) + len(kmer2) == 2*k - 1
            smallerH1, smallerH2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
            indelPairKmer.add( (smallerH1, smallerH2) )
            i+=1
        i=1
        while mid-i>=0 and initialH1[mid-i] == initialH1[mid]:
            kmer1 = seq[pos-(mid-1)-i:pos+(mid+2)-i]
            kmer2 = seq[pos-(mid-1)-i:pos+1] + seq[pos+2:pos+(mid+2)-i]
            assert len(kmer1) + len(kmer2) == 2*k - 1
            smallerH1, smallerH2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
            indelPairKmer.add( (smallerH1, smallerH2) )
            i+=1
    return    


def updata_non(s1_1, s2_1, pos_1, s1_2, s2_2, pos_2, seq, nonPairKmer):
    
    assert pos_1 < pos_2
    assert (seq[pos_1] == s1_1.upper() and 
            seq[pos_2] == s1_2.upper() )
    kmer1 = seq[pos_1-mid : pos_2+mid+1]
    kmer2 = seq[pos_1-mid:pos_1] + s2_1.upper() + seq[pos_1+1:pos_2] + s2_2.upper() + seq[pos_2+1:pos_2+mid+1]
    smallerKmer1, smallerKmer2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
    nonPairKmer.add( (smallerKmer1, smallerKmer2) )
    return

snpPairKmer = set()
indelPairKmer = set()
nonPairKmer = set()

for record in SeqIO.parse(faFile, "fasta"):
    ID = record.id.split('|')[0] 
    if ID in mutations:
        seq = record.seq.upper()
        mLen = len(mutations[ID])
        mutations[ID] = sorted(mutations[ID])
        if mLen == 1:   
            pos, s1, s2, homo = mutations[ID][0]
            updata_one_heterozygous_pair(s1, s2, pos, seq, snpPairKmer,  
                    indelPairKmer)
        elif mLen == 2:   
            pos_1, s1_1, s2_1, homo_1 = mutations[ID][0]
            pos_2, s1_2, s2_2, homo_2 = mutations[ID][1]
            pos_1, pos_2 = pos_1 - 1, pos_2 - 1
            # one mutation, vcf position starts from 0
            # more than on mutation, vcf position starts from 1 
            assert (len(s1_1) == 1 and len(s2_1) == 1)
            assert (len(s1_2) == 1 and len(s2_2) == 1)
            if abs(pos_2-pos_1) > mid:
                updata_one_heterozygous_pair(s1_1, s2_1, pos_1, 
                        seq, snpPairKmer, indelPairKmer)
                updata_one_heterozygous_pair(s1_2, s2_2, pos_2, 
                        seq, snpPairKmer, indelPairKmer)
            else:
                # updata non isolated
                updata_non(s1_1, s2_1, pos_1,
                           s1_2, s2_2, pos_2, seq, nonPairKmer)
        elif mLen==3:
            pos_1, s1_1, s2_1, homo_1 = mutations[ID][0]
            pos_2, s1_2, s2_2, homo_2 = mutations[ID][1]
            pos_3, s1_3, s2_3, homo_3 = mutations[ID][2]
            pos_1, pos_2, pos_3 = pos_1 - 1, pos_2 - 1, pos_3-1
            #print (ID, pos_1, pos_2, pos_3)
            # one mutation, vcf position starts from 0
            # more than on mutation, vcf position starts from 1 
            assert (len(s1_1) == 1 and len(s2_1) == 1)
            assert (len(s1_2) == 1 and len(s2_2) == 1)
            assert (len(s1_3) == 1 and len(s2_3) == 1)
            if abs(pos_2-pos_1) > mid:
                updata_one_heterozygous_pair(s1_1, s2_1, pos_1, 
                        seq, snpPairKmer, indelPairKmer)
            else:
                updata_non(s1_1, s2_1, pos_1,
                           s1_2, s2_2, pos_2, seq, nonPairKmer)
            if abs(pos_3-pos_2) > mid:    
                updata_one_heterozygous_pair(s1_2, s2_2, pos_2, 
                        seq, snpPairKmer, indelPairKmer)
                updata_one_heterozygous_pair(s1_3, s2_3, pos_3, 
                        seq, snpPairKmer, indelPairKmer)   
            else:
                # updata non isolated
                updata_non(s1_2, s2_2, pos_2,
                           s1_3, s2_3, pos_3, seq, nonPairKmer)
        else:    
            print ("need change code, mLen: ", mLen)
            sys.exit()


print ( "snp number:", len(snpPairKmer) )
sortedSnpPairKmer = sorted(list(snpPairKmer))
print ( "sorted snp number", len(sortedSnpPairKmer) )
snpOUT = open(snpFile, "w")
for (kmer1, kmer2) in sortedSnpPairKmer:
    snpOUT.write("%s %s\n" % (kmer1, kmer2) )
snpOUT.close()

sortedIndelPairKmer = sorted(indelPairKmer)
indelOUT = open(indelFile, "w")
for (kmer1, kmer2) in sortedIndelPairKmer:
    indelOUT.write("%s %s\n" % (kmer1, kmer2) )
indelOUT.close()

sortedNonPairKmer = sorted(nonPairKmer)
nonOUT = open(nonSnpFile, "w")
for (kmer1, kmer2) in sortedNonPairKmer:
    nonOUT.write("%s %s\n" % (kmer1, kmer2) )
nonOUT.close()

