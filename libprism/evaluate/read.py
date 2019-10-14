#########################################################################
# File Name: read.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 30 May 2019 15:34:53 AEST
#########################################################################
#!/bin/bash
from libprism.local import tools

def read_snp_sam(filename):
    #bothHappen, oneHappen, zeroHappen = 0, 0, 0
    bothHappen, oneHappen, zeroHappen = [], [], []
    count = 0
    align_0, align_1 = "", "" #alt_align_0, alt_align_1, = "","","",""
    with open(filename, "r") as f:
        line = f.readline()
        state = 0
        while line:
            if line.startswith("@"):
                line = f.readline()
                continue
            assert line.startswith("kmer")
            words = line.split()
            if state == 0:
                count +=1
                state = 1
                ID = words[0].split('_')[1]
                kmer_0 = words[9]
                refPos_0 = words[3]
                if words[-1].startswith("MD:Z:"):
                    align_0 = words[-1].split(":")[-1]  
                elif words[-1].startswith("XA:Z:"):
                    align_0 = words[-2].split(":")[-1]
                    alt_align_0 = words[-1].split(":")[-1]
                line = f.readline()
            elif state == 1:
                state = 2
                kmer_1 = words[9]
                refPos_1 = words[3]
                if words[-1].startswith("MD:Z:"):
                    align_1 = words[-1].split(":")[-1]
                elif words[-1].startswith("XA:Z:"):
                    align_1 = words[-2].split(":")[-1]
                    alt_align_1 = words[-1].split(":")[-1]
                line = f.readline()
            if state == 2:
                state = 0
                if align_0 == "31" and align_1 == "31":
                    print ("both exactly happen:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1)
                    bothHappen.append( (ID, refPos_0, refPos_1, kmer_0, kmer_1) )
                elif align_0 == "31":
                    print ("one exactly happen:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1)
                    oneHappend.append( (ID, refPos_0, kmer_0, kmer_1) )
                elif align_1 == "31":
                    print ("one exactly happen:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1)
                    oneHappend.append( (ID, refPos_1, kmer_0, kmer_1) )
                else:
                    print ("none exactly happen:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1)
                    zeroHappend.append( (ID, refPos_0, refPos_1, kmer_0, kmer_1) )

                '''
                if refPos_0 == refPos_1:
                    if ( (align_0 == "31" and ( align_1.startswith("15") and align_1.endswith("15") ) ) or
                         (align_1 == "31" and ( align_0.startswith("15") and align_0.endswith("15") ) ) ):
                        TP.append((ID, kmer_0, kmer_1))
                    else:
                        print "align to same position, but not 31/15*15:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1
                        print "alt align", alt_align_0
                        print "alt align", alt_align_1
                else:        
                    if align_0 == "31" and align_1 == "31":
                        print "both exactly happen:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1
                        print "alt align", alt_align_0
                        print "alt align", alt_align_1
                        bothHappen += 1
                    elif align_0 == "31" or align_1 == "31":
                        print "one exactly happen:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1
                        print "alt align", alt_align_0
                        print "alt align", alt_align_1
                        oneHappen += 1
                    elif align_0 != "31" and align_1 != "31":
                        print "none exactly happen:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1
                        print "alt align", alt_align_0
                        print "alt align", alt_align_1
                        zeroHappen += 1
                    else:    
                        print ID, align_0, align_1
                '''    
                    
    print (len(TP), zeroHappen, oneHappen, bothHappen)
    print (len(TP) + zeroHappen + oneHappen + bothHappen, count)

        
def read_indel_sam(filename):
    bothHappen, oneHappen, zeroHappen = 0, 0, 0
    TP = []
    TP_refPos = []
    count = 0
    indelCount = 0
    align_0, align_1, alt_align_0, alt_align_1, = "","","",""
    with open(filename, "r") as f:
        line = f.readline()
        state = 0
        while line:
            if line.startswith("@"):
                line = f.readline()
                continue
            assert line.startswith("kmer")
            words = line.split()
            if state == 0:
                count +=1
                state = 1
                ID = words[0].split('_')[1]
                kmer_0 = words[9]
                refPos_0 = words[3]
                if words[-1].startswith("MD:Z:"):
                    align_0 = words[-1].split(":")[-1]  
                elif words[-1].startswith("XA:Z:"):
                    align_0 = words[-2].split(":")[-1]
                    alt_align_0 = words[-1].split(":")[-1]
                line = f.readline()
            elif state == 1:
                state = 2
                kmer_1 = words[9]
                refPos_1 = words[3]
                if words[-1].startswith("MD:Z:"):
                    align_1 = words[-1].split(":")[-1]
                elif words[-1].startswith("XA:Z:"):
                    align_1 = words[-2].split(":")[-1]
                    alt_align_1 = words[-1].split(":")[-1]
                line = f.readline()
            if state == 2:
                state = 0
                if refPos_0 == refPos_1:
                        TP.append((ID, kmer_0, kmer_1))
                        TP_refPos.append(int(refPos_0))
                        assert len(kmer_0) == 31 
                        midCnt = tools.count_mid_same(kmer_0)
                        indelCount += 1.0/midCnt  
                else:        
                    print ("FP:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1)
                    print ("alt align", alt_align_0)
                    print ("alt align", alt_align_1)
                    '''
                    bothHappen += 1
                    elif align_0 == "31" or align_1 == "31":
                        print "one exactly happen:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1
                        print "alt align", alt_align_0
                        print "alt align", alt_align_1
                        oneHappen += 1
                    elif align_0 != "31" and align_1 != "31":
                        print "none exactly happen:", ID, align_0, align_1, "refpos:", refPos_0, refPos_1
                        print "alt align", alt_align_0
                        print "alt align", alt_align_1
                        zeroHappen += 1
                    else:    
                        print ID, align_0, align_1
                    '''    
                    
    print ("TP number:", len(TP))
    print ("FP number:", count-len(TP))
    print ("TP indel number:", indelCount)


def read_vcf(filename):
    f = open(filename, "r")
    mutations = {}
    ID = 1
    for line in f:
        if line.startswith('#'):
            continue
        words = line.split()
        pos = int(words[1])
        s1 = words[3].strip()
        s2 = words[4].strip()
        homo = words[9].split(':')[0]
        assert homo == "1/0" or homo == "1|0" or homo == "0|1" or homo == "0/1"
        s1Len = len(s1)
        s2Len = len(s2)
        #assert s1Len == 1 and s2Len == 1
        if s1 == '.' or s2 == '.':
            continue
        mutations[pos] = (s1, s2, ID)
        ID += 1

    f.close()     
    return mutations    


def read_group_kmer(filename):
    pair_kmer = []
    pair_group_length = []
    state = 0
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            if line.startswith('#'):
                continue
            if state == 0 and line.startswith("group"):
                pair_kmer.append(words[1:])
                state = 1
            elif state == 1:
                a = len(words)
                state = 2
            elif state == 2:
                b = len(words)
                state = 0
                pair_group_length.append( (a, b ) )

    return pair_kmer, pair_group_length      

def read_multip_columns(filename):
    pair_kmer = []
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            w = []
            for ele in words:
                if ele.isdigit():
                    w.append(int(ele))
                else:
                    w.append(ele)
            pair_kmer.append(w)
    print (pair_kmer[0])        
    return pair_kmer      


def read_2_int_columns(filename):
    ID2ID = {}
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split()
            ID2ID[ int(words[0]) ] = int(words[1])
    
    return ID2ID       


def read_2_columns(filename):
    kmers = {}
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split()
            kmers[ words[0] ] = int(words[1])
    
    return kmers    

def read1column(filename, pos):

    #ID = set()
    ID = list()
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            words = line.split()
            #ID.add( int( words[pos] ) )
            ID.append( int( words[pos] ))
    return ID

def get_diff(filename, pos):
    nums = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            words = line.split()
            nums.append( int( words[pos] ) )
    numsLen = len(nums)
    diff = []
    for i in range(numsLen-1):
        diff.append(nums[i+1]-nums[i])
    diffSorted = sorted(diff)    
    for ele in diffSorted:
        print (ele)



def read1column_str(filename):
    ID = set()
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            ID.add( words[0] )
    return ID


def read_matrix(filename):
    groupID = {}
    readID = 1
    reads = {}
    with open(filename, "r") as f:
        for line in f:
            if readID % 10000 == 0:
                print ("deal read number", readID)
            words = line.split()
            num = int(words[0])
            reads [ readID ] = {}
            for i in range(num):
                ID = int(words[2*i+2])
                reads[readID][ID] = int(words[2*i+3])
                if ID not in groupID:
                    groupID[ID] = set()
                groupID[ID].add(readID)
            readID += 1
    return groupID, reads     
            
