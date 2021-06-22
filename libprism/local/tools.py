#########################################################################
# File Name: tools.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 13:12:45 AEST
#########################################################################
#!/bin/bash

import sys
import os
import logging
#import subprocess

m={}
m['A'] = 'T'
m['T'] = 'A'
m['C'] = 'G'
m['G'] = 'C'

binTran={}
binTran['A'] = '00'
binTran['C'] = '01'
binTran['G'] = '10'
binTran['T'] = '11'

m2={}
m2['00'] = 'A'
m2['01'] = 'C'
m2['10'] = 'G'
m2['11'] = 'T'

RbinaryEle={}
RbinaryEle['00'] = '11'
RbinaryEle['01'] = '10'
RbinaryEle['10'] = '01'
RbinaryEle['11'] = '00'

def transfer_kmer_int(s):
    news=""
    for ele in s:
        news = news +  binTran[ele]  
    val = int(news,2)
    return val 

def transfer_int_kmer(val, k):
    s = bin(val)[2:].zfill(k*2)
    #lenS = len(s)
    return transfer_binary_string_2_kmer(s)   

def transfer_binary_string_2_kmer(s): # no exact kmer, nucleotide string

    kmer = ""
    lenS = len(s)
    for i in range(0, lenS, 2):
        kmer += m2[s[i:i+2]]
    return kmer
def reverse(s): # AAGTC
    news = ""
    for ele in s:
        news = m[ele] + news #GACTT
    return news   

def reverse_int(v, k):
    s = bin(v)[2:].zfill(2*k)
    news = reverse_binary_string(s)

    return int(news, 2)

def reverse_binary_string(s): #00 00 10 11 01
   
    news = ""
    lenS = len(s)
    if lenS %2 != 0:
        print ("binary string %s isnot correct len: %s" % (s, lenS) )
        logging.info("binary string %s isnot correct, len: %s " % (s, lenS) )
        sys.exit()
    for i in range(0, lenS, 2):
        news = RbinaryEle[ s[i:i+2] ] + news #10 00 01  11 11
    return news   


def reverse_ward(ward):
    if ward == 'f':
        return 'b'
    elif ward == 'b':
        return 'f'
    else:
        print ("reverse forward or backward error")
        sys.exit()

        
def calc_ATCG_zero_num(k1):
    l=[]
    l.append( k1.count('A') )
    l.append( k1.count('C') )
    l.append( k1.count('G') )
    l.append( k1.count('T') )

    return l.count(0)

def calc_extendLen(meanCov, heteRate):
    if heteRate >= 0.005:
        return 2
    else:
        return 1

def calc_weightThreshold(heteRate, k):
    if heteRate >= 0.005:
        return 4 #max( (k-1)/2 - heteRate-0.005, 4)
    else:
        return (k-3)/2


'''
s="ATCG"
print s
print reverse(s)
'''
def print_list(l):
    for ele in l:
        print (ele,)
    print ("")

def file_lines(filename):
    command = "wc -l " + filename + " >file_lines"
    os.system(command)
    with open("file_lines", "r") as f:
        for line in f:
            words = line.split()
            return int(words[0])


def hamming_distance(s1, s2):
    count = 0
    lenS = len(s1)
    #print s1, s2
    assert lenS == len(s2)
    for i in range(lenS):
        if s1[i] == -1 or s2[i] == -1:
            continue
        if s1[i] != s2[i]:
            count +=1
    return count

def hamming_distance2(s1, s2): # return pos
    pos = []
    lenS = len(s1)
    assert lenS == len(s2)
    for i in range(lenS):
        if s1[i] == -1 or s2[i] == -1:
            continue
        if s1[i] != s2[i]:
            pos.append(i)
    return pos
def get_smaller_kmer(kmer):

    RKmer = reverse(kmer)
    if kmer <= RKmer:
        return kmer
    else:
        return RKmer

def get_smaller_pair_kmer_keep_order(k1, k2):
    

    Rk1 = reverse(k1)
    Rk2 = reverse(k2)
    if (Rk1 < k1 and Rk1 < k2) : # the mini is Rk1 or Rk2
        k1,k2 = Rk1, Rk2    
    if (Rk2 < k1 and Rk2 < k2) : # the mini is Rk1 or Rk2
        k1,k2 = Rk1, Rk2
    return k1, k2



def get_smaller_pair_kmer(kmer1, kmer2):
   
    lenKmer1 = len(kmer1)
    lenKmer2 = len(kmer2) 
    if lenKmer1 < lenKmer2:
        kmer1, kmer2 = kmer2, kmer1
    RKmer1 = reverse(kmer1)
    RKmer2 = reverse(kmer2)
    if lenKmer1 == lenKmer2:
        if kmer2 < kmer1:
            kmer1, kmer2 = kmer2, kmer1
        if RKmer2 < RKmer1:
            RKmer1, RKmer2 = RKmer2, RKmer1
    if RKmer1 < kmer1:
        kmer1, kmer2 = RKmer1, RKmer2
    return kmer1, kmer2               

def get_smaller_pair_int(v1, v2, k):
   
    Rv1 = reverse_int(v1, k)
    Rv2 = reverse_int(v2, k)
    if v2 < v1:
        v1, v2 =v2, v1
    if Rv2 < Rv1:
        Rv1, Rv2 = Rv2, Rv1
    if Rv1 < v1:
        v1, v2 = Rv1, Rv2
    return v1, v2               

def get_smaller_pair_binary(v1, v2):
   
    Rv1 = reverse_binary_string(v1)
    Rv2 = reverse_binary_string(v2)
    if v2 < v1:
        v1, v2 =v2, v1
    if Rv2 < Rv1:
        Rv1, Rv2 = Rv2, Rv1
    if Rv1 < v1:
        v1, v2 = Rv1, Rv2
    return v1, v2               


def count_mid_same(kmer):
    ans = 1 
    l = len(kmer)
    mid = l/2
    i=1
    while kmer[mid+i] == kmer[mid]:
        ans += 1
        i+=1
    i=1
    while kmer[mid-i] == kmer[mid]:
        ans += 1
        i+=1
    return ans  


def min_edit_distance(word1, word2): 
    if not word1:
        return len(word2 or '') or 0

    if not word2:
        return len(word1 or '') or 0

    size1 = len(word1)
    size2 = len(word2)

    last = 0
    tmp = range(size2 + 1)
    value = None

    for i in range(size1):
        tmp[0] = i + 1
        last = i
        # print word1[i], last, tmp
        for j in range(size2):
            if word1[i] == word2[j]:
                value = last
            else:
                value = 1 + min(last, tmp[j], tmp[j + 1])
                # print(last, tmp[j], tmp[j + 1], value)
            last = tmp[j+1]
            tmp[j+1] = value
    #print value
    return value


def read_findGSE_result(filename):
    
    cnt = 0
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split()
            cnt += 1
            if cnt == 1:
                r = float(words[3])
            if cnt == 2:
                lowC = int(words[3])
            if cnt == 3:
                highC = int(words[3])

    extendLen = calc_extendLen((lowC+highC)/2, r) 
    lowC = lowC - extendLen
    highC = highC + extendLen
    return lowC, highC, r
