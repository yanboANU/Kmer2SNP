#########################################################################
# File Name: getKmerPair_FromVCF_REF.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Fri 18 Jun 2021 06:19:47 PM AEST
#input: .vcf and .fa
# for those unphased vcf file, like 0/1 not 0|1
#output: xx_k_31_pair.snp xx_k_31_pair.non.sep xx_k_31_pair.non.more 
#########################################################################
#!/bin/bash
import os
import sys

import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys
import tools
import time

# faster than transfer_hybrid.py
# consider combination of non-Isolated snp
# A/T 0/1 C/G 0/1
# A----C        A-----G
#          and             both include
# T----G        T-----C

def read_vcf(filename):
    snps = {}
    mutationNum = {}
    with open(filename, "r") as f: 
        for line in f:
            if line.startswith('#'):
                continue
            words = line.split()
            if len(words) < 10:
                continue
            ID = words[0]
            pos = int(words[1]) - 1
            s1 = words[3].strip()
            s2 = words[4].strip().split(',')[0]
            homo = words[9].split(':')[0]
            if homo == "0/0" or homo == "0|0" or homo == "1|1" or homo == "1/1" or homo == "./.":
                continue
            if homo == "1/2":
                continue
            assert homo == "1/0" or homo == "1|0" or homo == "0|1" or homo == "0/1"

            if ID not in snps:
                snps[ID] = []

            if s1 != '.' and s2 != '.':
                if ID not in mutationNum:
                    mutationNum[ID] = 1
                else:
                    mutationNum[ID] += 1
            if len(s1) == 1 and len(s2) == 1:
                snps[ID].append( (pos, s1, s2, homo) ) # single mutaion/ins/del
    return snps, mutationNum      

def get_first_last_pos(snps, seq, mid):

    firstPos = snps[0][0]
    lastPos = snps[-1][0]
    while firstPos - mid < 0:
        snps = snps[1 : ]
        if len(snps) == 0:
            return pairKmer
        firstPos = snps[0][0]
        print ("seq head not enough for left half k")
        #return pairKmer
    seqLen = len(seq)    
    while lastPos + mid + 1 > seqLen:
        snps = snps[ :-1]
        if len(snps) == 0:
            return pairKmer
        lastPos = snps[-1][0]
        print ("seq tail not enough for right half k")

    return firstPos, lastPos  


def get_pairKmer_from2seqs(kmer1, kmer2, mid, pairKmers, ID, ID2):

    dis = tools.hamming_distance2(kmer1, kmer2)
    #assert len(dis) == 1
    for pos in dis:
        h1, h2 = kmer1[pos-mid : pos+mid+1], kmer2[pos-mid : pos+mid+1] 
        smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
        if len(smallerH1) > 0 and  len(smallerH2) > 0:
            pairKmers.append( (smallerH1, smallerH2, ID, ID2) )
        else:
            print ("something wrong")
    #print ("pairKmers num", len(pairKmers))
    return     

def get_pairKmer_from_oneSNV(seg, seq, mid, pairKmers, ID):

    ID2 = ""
    for (pos, s1, s2, homo) in seg:
        if ID2 == "":
            ID2 = str(pos)
        else:
            ID2 = ID2 + "_" + str(pos)

        if len(s1) == 1 and len(s2) == 1:
            if seq[pos] != s1.upper():
                print ("snp error pos range", ID, pos, seq[pos-2:pos+3], seq[pos], s1, s2)
                print (seq)
                sys.exit()
            assert seq[pos] == s1.upper()
            seqLen = len(seq)
            if pos-mid >= 0 and pos+mid+1 < seqLen:
                kmer1 = seq[pos-mid:pos+mid+1]
                kmer2 = seq[pos-mid:pos] + s2.upper() + seq[pos+1:pos+mid+1]
                assert len(kmer1) + len(kmer2) == 2*k
                #print (kmer1)       
                #print (kmer2)        
                get_pairKmer_from2seqs(kmer1, kmer2, mid, pairKmers, ID, ID2)
    return


def get_pairKmer_from_muliSNV(seg, seq, mid, pairKmers, ID):

    ID2 = ""
    for (pos, s1, s2, homo) in seg:
        if ID2 == "":
            ID2 = str(pos)
        else:
            ID2 = ID2 + "_" + str(pos)

    if len(seg) == 2:
        pos_1, s1_1, s2_1, homo1 = seg[0]
        pos_2, s1_2, s2_2, homo2 = seg[1]

        assert pos_1 < pos_2
        assert (seq[pos_1] == s1_1.upper() and 
                seq[pos_2] == s1_2.upper() )
        if homo1 == homo2 and (homo1 == "0|1" or homo1 == "1|0"):
            kmer1 = seq[pos_1-mid : pos_2+mid+1]
            kmer2 = ( seq[pos_1-mid:pos_1] + s2_1.upper() + 
                    seq[pos_1+1:pos_2] + s2_2.upper() + seq[pos_2+1:pos_2+mid+1] ) 
            get_pairKmer_from2seqs(kmer1, kmer2, mid, pairKmers, ID, ID2)
            return
        elif ( (homo1 == "0|1" and homo2 == "1|0") or (homo2 == "0|1" and homo1 == "1|0") ): 
            kmer1 = seq[pos_1-mid : pos_2] + s2_2.upper() + seq[pos_2+1:pos_2+mid+1]
            kmer2 = seq[pos_1-mid : pos_1] + s2_1.upper() + seq[pos_1+1:pos_2+mid+1]
            get_pairKmer_from2seqs(kmer1, kmer2, mid, pairKmers, ID, ID2)
            return
        else:  
            kmer1 = seq[pos_1-mid : pos_2+mid+1]
            kmer2 = ( seq[pos_1-mid:pos_1] + s2_1.upper() + 
                    seq[pos_1+1:pos_2] + s2_2.upper() + seq[pos_2+1:pos_2+mid+1] ) 
            get_pairKmer_from2seqs(kmer1, kmer2, mid, pairKmers, ID, ID2)
            kmer1 = seq[pos_1-mid : pos_2] + s2_2.upper() + seq[pos_2+1:pos_2+mid+1]
            kmer2 = seq[pos_1-mid : pos_1] + s2_1.upper() + seq[pos_1+1:pos_2+mid+1]
            get_pairKmer_from2seqs(kmer1, kmer2, mid, pairKmers, ID, ID2)
            return
    elif len(seg) > 2:
        #print (seg)
         
        firstPos, lastPos = get_first_last_pos(seg, seq, mid)
        lenSNPs = len(seg)
        secondPos = seg[1][0]
        totalPair = 2**(lenSNPs-1)
        for i in range(totalPair):
            pos, s1, s2, homo = seg[0]      
            kmer1 = seq[firstPos - mid : secondPos ]
            kmer2 = seq[firstPos - mid : firstPos ] + s2.upper() + seq[firstPos+1: secondPos]  
            x = str(bin(i)[2:])
            x = '0' * (lenSNPs - 1 - len(x)) + x
            #print (x)
            for j in range(1, lenSNPs):
                pos, s1, s2, homo = seg[j]      
                #if seq[pos] != s1.upper():
                    #print (i, pos, seq[pos], s1.upper() )
                assert ( seq[pos] == s1.upper() ) 
                if j + 1 < lenSNPs:
                    pos_2, s1_2, s2_2, homo2 = seg[j+1]
                    if x[j-1] == '0':
                        kmer1 = kmer1 + s1.upper() + seq[pos+1: pos_2]
                        kmer2 = kmer2 + s2.upper() + seq[pos+1: pos_2]
                    else:    
                        kmer1 = kmer1 + s2.upper() + seq[pos+1: pos_2]
                        kmer2 = kmer2 + s1.upper() + seq[pos+1: pos_2]

                else:
                    if x[j-1] == '0':
                        kmer1 = kmer1 + s1.upper() + seq[pos+1: pos+mid+1]
                        kmer2 = kmer2 + s2.upper() + seq[pos+1: pos+mid+1]
                    else:
                        kmer1 = kmer1 + s2.upper() + seq[pos+1: pos+mid+1]
                        kmer2 = kmer2 + s1.upper() + seq[pos+1: pos+mid+1]

            #print (kmer1)       
            #print (kmer2)        
            get_pairKmer_from2seqs(kmer1, kmer2, mid, pairKmers, ID, ID2)
    return


def updata_one_pair(snps, seq, mid, ID):
    pairKmer = []

    firstPos, lastPos = get_first_last_pos(snps, seq, mid)
    #print ("lastPos", lastPos)
    lenSNPs = len(snps)
    prePos = firstPos
    seg = []
    seg.append(snps[0])
    i = 1
    snvNumber = 0
    snvPair = 0
    pairKmers = []
    while i < lenSNPs:
        pos, s1, s2, homo = snps[i]   
        while pos <= prePos + mid: # <=
            seg.append(snps[i])
            prePos = pos
            i += 1
            if i >= lenSNPs:
                break
            pos, s1, s2, homo = snps[i]
        if len(seg) == 1:
            get_pairKmer_from_oneSNV( sorted(seg) , seq , mid, pairKmers, ID)
            snvNumber += len(seg)
            snvPair += 1
            seg = []
        if i >= lenSNPs:
            break
        seg = []
        seg.append(snps[i])
        prePos = pos
        i += 1
    #assert len(seg) <= 1
    
    if len(seg) == 1:
        #print (seg)
        snvNumber += len(seg)
        snvPair += 1
        get_pairKmer_from_oneSNV( sorted(seg) , seq , mid, pairKmers, ID)
    
    #print ("isolated SNV number", snvNumber)
    #print ("isolated SNV pair", snvPair)
    #print ("Possible isolated SNV pair", len(pairKmers) )
    return pairKmers

def updata_more(snps, seq, mid, ID):
    pairKmer = []

    firstPos, lastPos = get_first_last_pos(snps, seq, mid)
    #print ("lastPos", lastPos)
    lenSNPs = len(snps)
    prePos = firstPos
    seg = []
    seg.append(snps[0])
    i = 1
    snvNumber = 0
    snvPair = 0
    pairKmers = []
    while i < lenSNPs:
        pos, s1, s2, homo = snps[i]   
        while pos <= prePos + mid:
            seg.append(snps[i])
            prePos = pos
            i += 1
            if i >= lenSNPs:
                break
            pos, s1, s2, homo = snps[i]   
        if len(seg) > 1:
            get_pairKmer_from_muliSNV( sorted(seg) , seq , mid, pairKmers, ID)
            snvNumber += len(seg)
            snvPair += 1
            seg = []
        if i >= lenSNPs:
            break
        seg = []
        seg.append(snps[i])
        prePos = pos
        i += 1
    #assert len(seg) <= 1
    
    if len(seg) > 1:
        #print (seg)
        snvNumber += len(seg)
        snvPair += 1
        get_pairKmer_from_muliSNV( sorted(seg) , seq , mid, pairKmers, ID)
    
    #print ("non-isolated SNV number", snvNumber)
    #print ("non-isolated SNV pair", snvPair)
    #print ("Possible isolated SNV pair", len(pairKmers) )
    return pairKmers


def get_isolated_pair(faFile, mutations, realMNum, k):
    
    snpPairKmer = [] 
    mid = int(k/2)
    for record in SeqIO.parse(faFile, "fasta"):
        ID = record.id.split('|')[0] 
        if ID in mutations:
            seq = record.seq.upper()
            mLen = len(mutations[ID])
            #print ("in contig/ref %s, mutations(include indel) number: %s real mutaion number %s " % (ID, mLen, realMNum[ID] ))
            mutations[ID] = sorted(mutations[ID])
            if mLen == 0:
                continue
            pairKmer = updata_one_pair(mutations[ID], seq, mid, ID) 
            snpPairKmer.extend( pairKmer )
    return snpPairKmer

def get_nonIsolated_pair(faFile, mutations, realMNum, k):
    
    snpPairKmer = [] 
    mid = int(k/2)
    for record in SeqIO.parse(faFile, "fasta"):
        ID = record.id.split('|')[0] 
        if ID in mutations:
            seq = record.seq.upper()
            mLen = len(mutations[ID])
            #print ("in contig/ref %s, mutations(include indel) number: %s real mutaion number %s " % (ID, mLen, realMNum[ID] ))
            mutations[ID] = sorted(mutations[ID])
            if mLen == 0:
                continue
            pairKmer = updata_more(mutations[ID], seq, mid, ID) 
            snpPairKmer.extend( pairKmer )
    return snpPairKmer

def write_isolated_snp_pair_kmer(snpPairKmer, k, preFilename):

    snpFile = preFilename +  "_k_" + str(k) + "_pair.snp"

    isol, nonSep, more = set(), set(), set()
    for (k1, k2, ID, pos) in snpPairKmer:
        dis = tools.hamming_distance(k1, k2)
        #print (k1, k2, dis)
        if dis == 1:
            isol.add( (k1, k2, ID, pos) )
        elif dis == 2:
            nonSep.add( (k1, k2, ID, pos) )
        elif dis > 2:
            more.add( (k1, k2, ID, pos) )
        else:
            print ("error", k1, k2, dis)
        
    assert ( len(nonSep)==0 and len(more)==0 )
    #print ("3 type snp size:", len(isol), len(nonSep), len(more) )
    
    print ( "isolated snp number:", len(isol) )
    sortedSnpPairKmer = sorted(list(isol))
    snpOUT = open(snpFile, "w")
    for (kmer1, kmer2, ID, pos) in sortedSnpPairKmer:
        snpOUT.write("%s %s %s %s\n" % (kmer1, kmer2, ID, pos) )
    snpOUT.close()
    return

def write_nonIsolated_snp_pair_kmer(snpPairKmer, k, preFilename):

    nonSepSnpFile = preFilename+ "_k_"+ str(k) + "_pair.non.sep"
    moreFile = preFilename + "_k_" + str(k) + "_pair.more"

    #print ( "snp pair kmer:", len(snpPairKmer) )
    #print ( "set snp pair kmer:", len(set(snpPairKmer) ) ) 
    # len(set(snpPairKmer) ) < len(snpPairKmer) is normal
    isol, nonSep, more = set(), set(), set()
    for (k1, k2, ID, pos) in snpPairKmer:
        #print (k1, k2, ID, pos)
        dis = tools.hamming_distance(k1, k2)
        #print (k1, k2, dis)
        if dis == 1:
            isol.add( (k1, k2, ID, pos) )
        elif dis == 2:
            nonSep.add( (k1, k2, ID, pos) )
        elif dis > 2:
            more.add( (k1, k2, ID, pos) )
        else:
            print ("error", k1, k2, dis)
        
    print ( "non-isolated sep snp, including all possible combination, number:", len(nonSep) )
    sortedKmer = sorted(list(nonSep))
    nonSepOUT = open(nonSepSnpFile, "w")
    for (kmer1, kmer2, ID, pos) in sortedKmer:
        nonSepOUT.write("%s %s %s %s\n" % (kmer1, kmer2, ID, pos) )
    nonSepOUT.close()

    print ( "more than 2 snps,including all possible combination, number:", len(more) )
    morePairKmer = sorted(list(more))
    moreOUT = open(moreFile, "w")
    for (kmer1, kmer2, ID, pos) in morePairKmer:
        moreOUT.write("%s %s %s %s\n" % (kmer1, kmer2, ID, pos) )
    moreOUT.close()
    return


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print ("ref.fa, kmer size, *.vcf")
        sys.exit()

    faFile = sys.argv[1]
    k=int(sys.argv[2])
    vcfFile = sys.argv[3]
    print ("input: %s %s %s" % (vcfFile, faFile, k) )
    t1 =time.time()
    mutations, realMNum = read_vcf(vcfFile)

    t2 =time.time()
    print ("read vcf cost %s s" % (t2-t1))
    print ( "chr number:", len(mutations) )
    isolatedPairKmer = get_isolated_pair(faFile, mutations, realMNum, k)

    t3 =time.time()
    print ("get isolated snp pair cost %s s" % (t3-t2))
    prefilename = vcfFile.split('.')[0]
    write_isolated_snp_pair_kmer(isolatedPairKmer, k, prefilename)


    nonPairKmer = get_nonIsolated_pair(faFile, mutations, realMNum, k)

    t4 =time.time()
    print ("get non-isolated snp pair cost %s s" % (t4-t3))
    write_nonIsolated_snp_pair_kmer(nonPairKmer, k, prefilename)
    
