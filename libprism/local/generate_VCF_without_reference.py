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

def generate_high_kmer(res, sampleID, vcfMap, prefix):
    
    kmers =set()
    for i in range(0, len(res), 2):
        ID1, kmer1, cnt1 = res[i]
        ID2, kmer2, cnt2 = res[i+1]
        if kmer1 > kmer2:
            ID1, kmer1, cnt1 = res[i+1]
            ID2, kmer2, cnt2 = res[i]
 
        if ID1.startswith(sampleID):
            key = kmer1+"/"+kmer2
            if key not in vcfMap:
                vcfMap[key] = prefix
            vcfMap[key] = vcfMap[key] + "GT=1:AD=" + str(cnt1) +","+ str(cnt2) + ";"
            kmers.add(kmer1)
            kmers.add(kmer2)
          
    return kmers


def generate_vcfMap(res, sampleID, homoCov, vcfMap, prefix):
    

    kmers = generate_high_kmer(res, sampleID, vcfMap, prefix)
    heteCov = homoCov/2
    
    for i in range(0, len(res), 2):
        ID1, kmer1, cnt1 = res[i]
        ID2, kmer2, cnt2 = res[i+1]

        if kmer1 > kmer2:
            ID1, kmer1, cnt1 = res[i+1]
            ID2, kmer2, cnt2 = res[i]

        if cnt1 == 0 and cnt2 == 0:
            continue
        if kmer1 not in kmers and kmer2 not in kmers:
            if (cnt1==0 and cnt2 >5) or (cnt2 > heteCov*1.6 and cnt1 < 5):
                
                key = kmer1+"/"+kmer2
                if key not in vcfMap:
                    vcfMap[key] = prefix
                vcfMap[key] = vcfMap[key] + "GT=0:AD=" +str(cnt1) +","+str(cnt2) + ";"
                continue

            elif (cnt2==0 and cnt1 >5) or (cnt1 > heteCov*1.6 and cnt2 < 5): 
                key = kmer1+"/"+kmer2
                if key not in vcfMap:
                    vcfMap[key] = prefix
                vcfMap[key] = vcfMap[key] + "GT=0:AD=" + str(cnt1) +","+ str(cnt2) + ";"     
                continue
            elif (cnt1 > 2 and cnt2 > 2):
                key = kmer1+"/"+kmer2
                if key not in vcfMap:
                    vcfMap[key] = prefix
                vcfMap[key] = vcfMap[key] + "GT=1:AD=" + str(cnt1) +","+ str(cnt2) + ";"     
                continue
            else:
                print (ID1, ID2, kmer1, kmer2, cnt1, cnt2)

    return






def write_fastq(res, sampleID, homoCov, filename):
    

    kmers, records = write_hete_kmer(res, sampleID)
    if len(kmers) == 0:
        return 
    heteCov = homoCov/2
    index = int( (len(kmers)/2) + 1)
    for i in range(0, len(res), 2):
        ID1, kmer1, cnt1 = res[i]
        ID2, kmer2, cnt2 = res[i+1]
        if cnt1 == 0 and cnt2 == 0:
            continue
        if kmer1 not in kmers and kmer2 not in kmers:
            #if cnt1 < cnt2 and cnt1 < 5 and (cnt2 > heteCov*1.6 or (cnt1==0 and cnt2 >heteCov*1.2) ):
            if (cnt1==0 and cnt2 >5) or (cnt2 > heteCov*1.6 and cnt1 < 5):
                ID1 = sampleID + ":" + str(index) + ":homo:0:" + str(cnt2) +":1"
                ID2 = sampleID + ":" + str(index) + ":homo:0:" + str(cnt2) +":2"
                r0 = SeqRecord( Seq(kmer2) , id= ID1 )
                r0.letter_annotations["phred_quality"] = [43]*len(kmer2)
                r1 = SeqRecord( Seq(kmer2) , id= ID2 )
                r1.letter_annotations["phred_quality"] = [43]*len(kmer2)
                records.append(r0)
                records.append(r1) 
                kmers.add(kmer2)
                index += 1
                continue
            #elif cnt1 > cnt2 and cnt2 < 5 and ( cnt1 >heteCov*1.6 or (cnt2==0 and cnt1 >heteCov*1.2) ) :
            elif (cnt2==0 and cnt1 >5) or (cnt1 > heteCov*1.6 and cnt2 < 5): 
                ID1 = sampleID + ":" + str(index) + ":homo:0:" + str(cnt1) +":1"
                ID2 = sampleID + ":" + str(index) + ":homo:0:" + str(cnt1) +":2"
                r0 = SeqRecord( Seq(kmer1) , id= ID1 )
                r0.letter_annotations["phred_quality"] = [43]*len(kmer1)
                r1 = SeqRecord( Seq(kmer1) , id= ID2 )
                r1.letter_annotations["phred_quality"] = [43]*len(kmer1)
                records.append(r0)
                records.append(r1)
                kmers.add(kmer1)
                index += 1
                continue
            #elif cnt1 > heteCov*0.6 and cnt1 < heteCov*1.2 and cnt2 > heteCov*0.6 and cnt1 < heteCov*1.2:
            elif (cnt1 > 2 and cnt2 > 2):
                ID1 = sampleID + ":" + str(index) + ":het_other:0:"+ str(cnt1) + ID1[-2:]
                ID2 = sampleID + ":" + str(index) + ":het_other:0:"+ str(cnt2) + ID2[-2:]
                r0 = SeqRecord( Seq(kmer1) , id= ID1 )
                r0.letter_annotations["phred_quality"] = [43]*len(kmer1)
                r1 = SeqRecord( Seq(kmer2) , id= ID2 )
                r1.letter_annotations["phred_quality"] = [43]*len(kmer2)
                records.append(r0)
                records.append(r1)
                kmers.add(kmer1)
                kmers.add(kmer2)
                index += 1
                continue
            else:
                print (ID1, ID2, kmer1, kmer2, cnt1, cnt2)

    print ("total kmer number", len(records))        
    SeqIO.write(records, filename, "fastq")
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

def run(filename):

    covMap = read_coverage(filename)
    kmerID = []
    for sampleID in covMap:
        temp = read_hete_kmer(sampleID + "/k_31_pair.snp", sampleID)
        kmerID.extend( temp )
    kmerID = sorted(kmerID) 
    vcfMap = {}
    prefix = ""
    for sampleID in covMap:
        filename = sampleID + "/chr_k31.txt"
        kmerFreq =  read_kmer_freq(filename)
        homoCov = covMap[ sampleID] 
        res = compare_two_lists(kmerID, kmerFreq)
        #print (res)
        generate_vcfMap(res, sampleID, homoCov, vcfMap, prefix)
        prefix = prefix + "-;"

        #write_fastq(res, sampleID, homoCov, sampleID + "_kmerPair.fastq")

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


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print ("covList")
        sys.exit()
   
    covMap = read_coverage(sys.argv[1])
    kmerID = []
    for sampleID in covMap:
        temp = read_hete_kmer(sampleID + "/k_31_pair.snp", sampleID)
        kmerID.extend( temp )
    kmerID = sorted(kmerID) 
    vcfMap = {}
    prefix = ""
    for sampleID in covMap:
        filename = sampleID + "/chr_k31.txt"
        kmerFreq =  read_kmer_freq(filename)
        homoCov = covMap[ sampleID] 
        res = compare_two_lists(kmerID, kmerFreq)
        #print (res)
        generate_vcfMap(res, sampleID, homoCov, vcfMap, prefix)
        prefix = prefix + "-;"

        #write_fastq(res, sampleID, homoCov, sampleID + "_kmerPair.fastq")

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
