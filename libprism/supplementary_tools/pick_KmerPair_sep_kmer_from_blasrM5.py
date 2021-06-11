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
            for e in eles:
                f.write("%s " % (e) )
            f.write("\n")    

def write_sep_kmer(h1Kmers, h2Kmers , k):

    fout1 = open("h1_" + str(k) + "mers" ,"w")
    fout2 = open("h2_" + str(k) + "mers" ,"w")
    h1Kmers = sorted(h1Kmers.items())
    h2Kmers = sorted(h2Kmers.items())
    
    for (kmer, IDs) in h1Kmers:
        if IDs.count(',') > 0:
            continue
        fout1.write("%s %s\n" %  (kmer, IDs))

    for (kmer, IDs) in h2Kmers:
        if IDs.count(',') > 0:
            continue
        fout2.write("%s %s\n" %  (kmer, IDs))

    fout1.close()    
    fout2.close()    


def read_one_blasr_m5_line(line):
    
    words = line.split()

    qname = words[0].split('_')[0]
    qlen = int(words[1])
    qStart = int(words[2])
    qEnd = int(words[3])
    tname = words[5].split('_')[0]
    if qname[:-1] != tname[:-1] or (qEnd-qStart)/float(qlen)<0.9:
        return "","",[], ""
    chrID = qname[:-1]
    print (qname, tname)
    print ("match number", words[-8])
    mis = int(words[-7]) 
    ins = int(words[-6])
    dele = int(words[-5])
    print (chrID , "mismatch/insert/delete number:", mis, ins, dele)
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
    return qAlignedSeq, tAlignedSeq, misPos, chrID 


# in this version, blasr.m5 file can include multiple chromosome alignment

if __name__ == '__main__':
    
    if len(sys.argv) < 2:
        print ("ksize, *blasr.m5")
        ##############
        #reference name must in chr1A_xxxx/1A_xxx and chr1B_xxx/1B_xxx format
        ########
        sys.exit()

    k = int( sys.argv[1] )
    mid = int(k/2)
    filename = sys.argv[2]

    kmers, nons, mores = [], [], []
    snpNumIn1kmer = {}
    mis, ins, dele = 0, 0, 0
    h1Kmers, h2Kmers = {}, {}

    with open(filename, "r") as f:
     
        for line in f:
            if line.startswith('[INFO]') or len(line.strip())<1:
                continue
            qAlignedSeq, tAlignedSeq, misPos, chrID = read_one_blasr_m5_line(line)
            if len(qAlignedSeq) == 0:
                continue
            i = 0
            misNum = len(misPos)
            while i < misNum:
                key = misPos[i]
                h1 = qAlignedSeq[key-mid : key+mid+1 ] # 0
                h2 = tAlignedSeq[key-mid : key+mid+1 ] # 1
                if h1.count('N') > 0 or h2.count('N') > 0 or h1.count('-')>0 or h2.count('-')>0:
                    i += 1
                    continue                 
                k1, k2 = tools.get_smaller_pair_kmer_keep_order(h1, h2)
                k1Key = chrID + "_" + str(key) + "_0"
                if k1 not in h1Kmers:
                    h1Kmers[k1] = k1Key 
                else:
                    h1Kmers[k1] = h1Kmers[k1] +"," + k1Key
                
                k2Key = chrID + "_" + str(key) + "_1"
                if k2 not in h2Kmers:
                    h2Kmers[k2] = k2Key 
                else:
                    h2Kmers[k2] = h2Kmers[k2] +"," + k2Key
                
                #fout1.write("%s %s\n" % (k1, chrID + "_" + str(key) + "_0") )
                #fout2.write("%s %s\n" % (k2, chrID + "_" + str(key) + "_1") )

                h1, h2 = tools.get_smaller_pair_kmer(h1, h2)
                dis = tools.hamming_distance(h1, h2)
                if dis == 1:
                    kmers.append( (h1, h2, chrID, key) )    
                elif dis == 2:
                    nons.append( (h1, h2, chrID, key) )
                else:
                    if dis not in snpNumIn1kmer:
                        snpNumIn1kmer[dis] = 1
                    else:
                        snpNumIn1kmer[dis] += 1
                    mores.append( (h1, h2, chrID, key) )
                i += 1    

    print ("snp number differ in one kmer", sorted(snpNumIn1kmer.items()))
    write_pair_kmer("k_" + str(k)+"_pair.snp", kmers)
    write_pair_kmer("k_" + str(k)+"_pair.non.sep", nons)
    write_pair_kmer("k_" + str(k)+"_pair.more", mores)
    write_sep_kmer(h1Kmers, h2Kmers , k)

'''
# remove this function
def write_sep_kmer(kmers, nons, mores, chrID, ksize):

    group1, group2 = [] , []
    posS = []
    for (k1, k2, name, pos) in kmers:
        if tools.reverse(k1) < k1:
            k1 = tools.reverse(k1)
        if tools.reverse(k2) < k2:
            k2 = tools.reverse(k2)
        group1.append( (k1, pos) )
        group2.append( (k2, pos) )
        posS.append(pos)


    for (k1, k2, name, pos) in nons:
        if tools.reverse(k1) < k1:
            k1 = tools.reverse(k1)
        if tools.reverse(k2) < k2:
            k2 = tools.reverse(k2)
        
        group1.append( (k1, pos) )
        group2.append( (k2, pos) )
    sortedGroup1 = sorted(group1)
    sortedGroup2 = sorted(group2)
    
    lenPos = len(posS)
    sortedPos = sorted(posS)
    disDistribution= {}
    for i in range(lenPos - 1):
        dis = int( (sortedPos[i+1] - sortedPos[i])/100 )
        if dis in disDistribution:
            disDistribution[dis] += 1
        else:
            disDistribution[dis] = 1
        
    #print ("distance distribution")
    #sortedDisD = sorted(disDistribution.items())
    #for (dis, fre) in sortedDisD:
        #print ("[%s %s] %s" % ( dis*100, (dis+1)*100 - 1, fre ) )
    

    fout = open("h1_" + str(ksize) + "mers" ,"w")
    for (k, pos) in sortedGroup1:
        fout.write("%s %s\n" % (k, chrID + "_" + str(pos) + "_0") )
    fout.close()    

    fout = open("h2_" + str(ksize) + "mers" ,"w")
    for (k, pos) in sortedGroup2:
        fout.write("%s %s\n" % (k, chrID + "_" + str(pos) + "_1") )
    fout.close()    
    return
'''
