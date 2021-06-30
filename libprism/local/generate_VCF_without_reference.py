#########################################################################
# File Name: count_each_heterozygous_kmer.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 08 Jun 2021 12:24:01 PM AEST
#########################################################################
#!/bin/bash
import os
import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


m={}
m['A'] = 'T'
m['T'] = 'A'
m['C'] = 'G'
m['G'] = 'C'

def reverse(s): # AAGTC
    news = ""
    for ele in s:
        news = m[ele] + news #GACTT
    return news   

def get_smaller_kmer(kmer):

    RKmer = reverse(kmer)
    if kmer <= RKmer:
        return kmer
    else:
        return RKmer

def read_kmer_freq(filename):
    kmerFreq = []
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split()
            cnt = int(words[1])
            if cnt <= 2:
                continue
            kmer  = words[0]
            kmer = get_smaller_kmer(kmer)
            kmerFreq.append( (kmer, cnt) )
    return sorted(kmerFreq)

def read_hete_kmer(filename, sampleID):
    kmerID = []
    cnt = 1
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            k1 = words[0]
            k2 = words[1]
            weight = words[2]
            ID1  = sampleID + "_" + str(cnt) +  "_w" + weight + "_1"
            ID2  = sampleID + "_" + str(cnt) +  "_w" + weight + "_2" 
            kmerID.append( (k1, ID1) )
            kmerID.append( (k2, ID2) )
            cnt += 1
    return sorted(kmerID) 


def read_edgeWeight(filename):
    ew = {} #edge weight
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            k1 = words[0]
            k2 = words[1]
            weight = words[2]
            #assert k1 < k2
            ew[(k1,k2)] = weight
    return ew


def compare_two_lists(list2, list4):

    res=[]
    index2, index4=0,0
    len2, len4 = len(list2), len(list4)

    while index2 < len2 and index4 < len4:
        if (list2[index2][0:1] < list4[index4][0:1]):
            #print (list2[index2][1], 0)
            res.append( (list2[index2][1], list2[index2][0], 0  ) )
            index2 +=1
        elif list2[index2][0:1] > list4[index4][0:1]:   
            index4 +=1
        elif list2[index2][0:1] == list4[index4][0:1]: 
            #print (list2[index2])
            #print (list4[index4])
            #print (list2[index2][1], list4[index4][1]) 
            res.append( (list2[index2][1], list2[index2][0], list4[index4][1] )  )
            index2 +=1
            #index4 +=1 #list2 have repeat kmer
        else: 
            print ("unknow")
            print (list2[index2])
            print (list4[index4])
            index2 +=1
            index4 +=1

    return sorted(res)

def generate_high_kmer(res, sampleID, vcfMap, prefix, edgeWeight):
    
    kmers =set()
    kLen = len(res[0][1])
    mid = int(kLen/2)
    for i in range(0, len(res), 2):
        ID1, kmer1, cnt1 = res[i]
        ID2, kmer2, cnt2 = res[i+1]
        if kmer1 > kmer2:
            ID1, kmer1, cnt1 = res[i+1]
            ID2, kmer2, cnt2 = res[i]
  
        if ID1.startswith(sampleID):
            #key = kmer1+"/"+kmer2
            w = edgeWeight[ (kmer1, kmer2) ]
            allele1, allele2  = kmer1[ mid ], kmer2[ mid ] 
            key = kmer1[ : mid ] + '[' + allele1 + "/" + allele2 + ']' + kmer1[mid+1 : ]
            if key not in vcfMap:
                vcfMap[key] = prefix
            #GT:AD:EW (edge weight)      
            vcfMap[key] = vcfMap[key] + "0/1:" + str(cnt1) +","+ str(cnt2) + ":"+ w +";"
            kmers.add(kmer1)
            kmers.add(kmer2)
          
    return kmers, mid


def generate_vcfMap(res, sampleID, homoCov, vcfMap, prefix, ew):
    

    kmers, mid = generate_high_kmer(res, sampleID, vcfMap, prefix, ew)
    heteCov = homoCov/2
    
    for i in range(0, len(res), 2):
        ID1, kmer1, cnt1 = res[i]
        ID2, kmer2, cnt2 = res[i+1]

        if kmer1 > kmer2:
            ID1, kmer1, cnt1 = res[i+1]
            ID2, kmer2, cnt2 = res[i]

        if cnt1 + cnt2 <=5:
            continue
        allele1, allele2  = kmer1[ mid ], kmer2[ mid ] 
        key = kmer1[ : mid ] + '[' + allele1 + "/" + allele2 + ']' + kmer1[mid+1 : ]
        if kmer1 not in kmers and kmer2 not in kmers: 
            if ( cnt1/(cnt1+cnt2) > 0.95 ):    
                if key not in vcfMap:
                    vcfMap[key] = prefix
                vcfMap[key] = vcfMap[key] + "0/0:" +str(cnt1) +","+str(cnt2) + ":-1;"
                kmers.add(kmer1)
                kmers.add(kmer2)
                continue
            elif ( cnt2/(cnt1+cnt2) > 0.95 ): 
                if key not in vcfMap:
                    vcfMap[key] = prefix
                vcfMap[key] = vcfMap[key] + "1/1:" + str(cnt1) +","+ str(cnt2) + ":-1;"    
                kmers.add(kmer1)
                kmers.add(kmer2)
                continue
            elif (cnt1 > 2 and cnt2 > 2):
                if key not in vcfMap:
                    vcfMap[key] = prefix
                vcfMap[key] = vcfMap[key] + "0/1:" + str(cnt1) +","+ str(cnt2) + ":-1;"    
                kmers.add(kmer1)
                kmers.add(kmer2)
                continue
            else:
                print (ID1, ID2, kmer1, kmer2, cnt1, cnt2)

    return

def read_coverage(filename):
    covMap = {}
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split(' ')
            sampleID = words[0]
            cov = int(words[1])
            covMap[sampleID] = cov         
    return covMap

def check_vcfMap(vcfMap, cnt):

    for key in vcfMap:
        if vcfMap[key].count(';') == cnt:
            continue
        elif vcfMap[key].count(';') == cnt-1:
            vcfMap[key] = vcfMap[key] + "-;"
            continue
        else:
            print ("error in vcfMap")
    return        


def run(filename):

    covMap = read_coverage(filename)
    kmerID = []
    for sampleID in covMap:
        temp = read_hete_kmer(sampleID + "/k_31_pair.snp", sampleID)
        kmerID.extend( temp )
    kmerID = sorted(kmerID) # have repeat kmer, different ID 
    vcfMap = {}
    prefix = ""
    sampleCnt = 0
    for sampleID in covMap:
        
        edgeWeight = read_edgeWeight(sampleID + "/k_31_pair.snp")
        filename = sampleID + "/chr_k31.txt"
        kmerFreq =  read_kmer_freq(filename)
        homoCov = covMap[ sampleID] 
        res = compare_two_lists(kmerID, kmerFreq)
        #print (res)
        generate_vcfMap(res, sampleID, homoCov, vcfMap, prefix, edgeWeight)
        sampleCnt += 1
        if sampleCnt > 2:
            check_vcfMap(vcfMap, sampleCnt)
        prefix = prefix + "-;"

    fout = open("all_samples.vcf", "w")
    fout.write("kmer\t")
    for sampleID in covMap:
        fout.write("%s\t" % (sampleID))
    fout.write("\n")
    for key in vcfMap:
        print (key, vcfMap[key])
        fout.write("%s\t" % key)
        words = vcfMap[key][:-1].split(";")
        for w in words:
            fout.write("%s\t" % (w))
        fout.write("\n")

    fout.close()    

'''
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print ("covList")
        sys.exit()
   
    run(sys.argv[1])
'''    


