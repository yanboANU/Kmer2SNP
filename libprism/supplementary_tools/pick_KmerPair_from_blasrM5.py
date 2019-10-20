#/*************************************************************************
#	> File Name: pick_KmerPair_from_blasrM5.py
#	> Author: Yanbo Li
#	> Mail: liyanbotop@163.com 
#	> Created Time: Sat 19 Oct 2019 10:47:14 AEDT
# ************************************************************************/
import os
import sys
import tools

def write_pair_kmer(outFile, kmers):

    sortedKmers = sorted(kmers)
    with open(outFile, "w") as f:
        for eles in sortedKmers:
            #f.write("%s %s %s %s %s %s %s %s\n" % (ele[0], ele[1], ele[2], ele[3], ele[4], ele[5], tools.reverse(ele[0]), tools.reverse(ele[1]) ) )
            for e in eles:
                f.write("%s " % (e) )
            f.write("\n")    


infile = ['chr1A_B.blasr.m5', 'chr3A_B.blasr.m5', 'chr5A_B.blasr.m5', 'chr7A_B.blasr.m5',
          'chr2A_B.blasr.m5', 'chr4A_B.blasr.m5', 'chr6A_B.blasr.m5', 'chrRA_B.blasr.m5']

k = int(sys.argv[1])
mid = int(k/2)

kmers, nons, mores = [], [], []
snpNumIn1kmer={} 
for onefile in infile:
    with open(onefile, "r") as f:
        chrID = onefile[0:4]
        for line in f:
            if line.startswith('[INFO]'):
                continue
            words = line.split()
            print chrID , "mismatch number:", words[-7]
            qAlignedSeq = words[-3]
            matchPattern = words[-2]
            tAlignedSeq = words[-1]
            pairNum = len(matchPattern)
            #print matchPattern.count('|') # only two pattern * and |
            #print "check mismatch num", len(matchPattern) - matchPattern.count('|')
            #print "check mismatch num", matchPattern.count('*')
            misPos = []
            for i in range(pairNum):
                if matchPattern[i] !=  '|':
                    assert matchPattern[i] == '*'
                    #misPos[i] = matchPattern[i]
                    misPos.append(i)
            i = 0
            misNum = len(misPos)
            while i < misNum:
                key = misPos[i]
                h1 = qAlignedSeq[key-mid-1 : key+mid ] # 0
                h2 = tAlignedSeq[key-mid-1 : key+mid ] # 1
                if h1.count('N') > 0 or h2.count('N') > 0 or h1.count('-')>0 or h2.count('-')>0:
                    i += 1
                    continue
                smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
                dis = tools.hamming_distance(smallerH1, smallerH2)
                if dis == 1:
                    kmers.append( (smallerH1, smallerH2, chrID, key) )    
                elif dis == 2:
                    nons.append( (smallerH1, smallerH2, chrID, key) )
                else:
                    if dis not in snpNumIn1kmer:
                        snpNumIn1kmer[dis] = 1
                    else:
                        snpNumIn1kmer[dis] += 1
                    mores.append( (smallerH1, smallerH2, chrID, key) )
                i += 1    

print "snp number differ in one kmer", sorted(snpNumIn1kmer.items())
write_pair_kmer("kmer_pair_k" + str(k), kmers)
write_pair_kmer("non_pair_k" + str(k), nons)
write_pair_kmer("more_pair_k" + str(k), mores)





