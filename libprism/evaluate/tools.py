#########################################################################
# File Name: tools.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 13:12:45 AEST
#########################################################################
#!/bin/bash

import sys
import os
#import subprocess

m={}
m['A'] = 'T'
m['T'] = 'A'
m['C'] = 'G'
m['G'] = 'C'

def reverse(s): # AATC
    news = ""
    for ele in s:
       news = m[ele] + news #CGATT
    return news   

def reverse_ward(ward):
    if ward == 'f':
        return 'b'
    elif ward == 'b':
        return 'f'
    else:
        print ("reverse forward or backward error")
        sys.exit()

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
