#########################################################################
# File Name: kmercalling.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Fri 09 Aug 2019 11:40:22 AEST
#########################################################################
#!/bin/bash

import os
import sys
import time
from libprism.local import tools
import logging
#from tools import *

def pick_smaller_unique_kmer(input_filename, low, high):
    uniqKmer = {}
    with open(input_filename, "r") as f:
        for line in f:
            words = line.strip().split()
            coverage = int(words[1])
            if coverage < low or coverage > high:
                continue 
            kmer = words[0]
            newkmer = tools.reverse(kmer)
            if kmer > newkmer:    
                kmer = newkmer
            val = tools.transfer_kmer_int(kmer) 
            uniqKmer[ val ] = coverage
    logging.debug("uniqKmer memory %s bytes" % (sys.getsizeof(uniqKmer)) )
    return uniqKmer


def build_index(uniqKmer, k):
    
    left_index , right_index = {}, {}
    m_index = {} # remove middle index 
    count = 0
    filterLKey = set()
    filterRKey = set()

    for intkmer in iter(uniqKmer):
        binarykmer = bin(intkmer)[2:].zfill(k*2) 
        
        key = int(binarykmer[:k-1] + binarykmer[k+1:], 2)
        if key not in m_index:
            m_index[key] = []
        m_index[key].append( intkmer )
    
        key = int(binarykmer[:-2], 2)
        if key not in left_index:
            left_index[ key  ] = int(binarykmer[-2:], 2) #intkmer #int(binarykmer[-2:]) # intkmer
        else:
            filterLKey.add(key)
        
        key = int(binarykmer[2:], 2)
        if key not in right_index:
            right_index[ key ] = int(binarykmer[:2], 2) #intkmer #int(binarykmer[:2]) # intkmer
        else:
            filterRKey.add(key)

    for key in iter(filterLKey):
        left_index.pop(key)
    for key in iter(filterRKey):
        right_index.pop(key)

    logging.debug( "left_index Memory %s bytes " % ( sys.getsizeof(left_index) ) )
    logging.debug( "right_index Memory %s bytes " % ( sys.getsizeof(right_index) ) )
    logging.debug( "m_index Memory %s bytes " % (sys.getsizeof(m_index) ) )
    return left_index, right_index, m_index


'''
def build_map_merge(left, k):
    mapMerge= {}
    mid = int(k/2)
    hisMap = {}
    for key in left:
        groupSize = len(left[key])
        for i in range(0, groupSize-1):
            for j in range(i+1, groupSize):
                k1 = left[key][i]
                k2 = left[key][j]
                if k1[mid] == k2[mid]: #or (k1 in mappedKmer) or (k2 in mappedKmer):
                    continue
                dis = 0
                cnt = mid + 1
                diffPos = 0
                while cnt < k:
                    if k1[cnt] != k2[cnt]:
                        dis += 1
                        diffPos = cnt
                    cnt += 1    
                    if dis >= 2:
                        break
                if cnt == k and dis == 1:
                    key1, key2 = k1[ diffPos - mid : ] , k2[ diffPos - mid : ]
                    mink1, mink2 = tools.get_smaller_pair_kmer(k1, k2)
                    #fout.write("%s %s %s %s\n" % (mink1, mink2, cov1, cov2) )
                    #candidateNonPair.append((mink1, mink2))
                    if mink1 not in hisMap:
                        hisMap[mink1] = 0
                    if mink2 not in hisMap:
                        hisMap[mink2] = 0
                    hisMap[mink1] += 1
                    hisMap[mink2] += 1
                    update_map_merge(left[key][i], left[key][j], key1, key2, mapMerge)
                    #mappedKmer.add(k1)
                    #mappedKmer.add(k2)
                    #break #3 lines add 22 Aug. a kmer only allow one kmer hamming distance equal to 2
    
    highRepeat = set()
    for key in hisMap:
        if hisMap[key] > 2:
            highRepeat.add(key)
    print ("high Repeat kmer number", len(highRepeat) )
    return mapMerge, highRepeat

def merge_pair(mapMerge, highRepeat, k, left_index, right_index, uniqKmers, kmerCov):
    pairSet, nonPair = set(), []
    usedKmers = []
    for (key1, key2) in mapMerge:
        mlen = len( mapMerge [ (key1, key2) ] )
        if mlen > 2 or mlen==1: # one key only allow a pair
            continue
        Left1, Left2 = mapMerge[ (key1, key2) ][0]
        Right1, Right2 = mapMerge[ (key1, key2) ][1]
        l = len(key1)

        Lmin1 = min(Left1, tools.reverse(Left1))
        Lmin2 = min(Left2, tools.reverse(Left2))
        Rmin1 = min(Right1, tools.reverse(Right1))
        Rmin2 = min(Right2, tools.reverse(Right2))
        covL1, covL2 = kmerCov[Lmin1], kmerCov[Lmin2]
        covR1, covR2 = kmerCov[Rmin1], kmerCov[Rmin2]
        if (Lmin1 in highRepeat or Lmin2 in highRepeat or
                Rmin1 in highRepeat or Rmin2 in highRepeat):
            continue
        if covL1 > 1.5*covR1 or 1.5*covL1 < covR1:
            continue
        if covL2 > 1.5*covR2 or 1.5*covL2 < covR2:
            continue
        if Right1[-l:] == Left1[:l] and Right2[-l:] == Left2[:l]:
            merge1 = Right1 + Left1[l:]
            merge2 = Right2 + Left2[l:]
        elif Left1[-l:] == Right1[:l] and Left2[-l:] == Right2[:l]:
            merge1 = Left1 + Right1[l:]
            merge2 = Left2 + Right2[l:]
        else:
            continue #print ("one side")
         
        #TTTTTTTTTTTTTTTCAAAAAAAAAAAAAAAA # also useful SNP and indel 
        #TTTTTTTTTTTTTTTTCAAAAAAAAAAAAAAA 
        if merge1[1:] == merge2[:-1] or merge1[:-1]==merge2[1:]:
            print (merge1, merge2)
            continue
        small1, small2 = tools.get_smaller_pair_kmer(merge1, merge2) 
        ek1, ek2, flag, group1, group2 = extend_one_pair(small1, small2, left_index, right_index, k, 0)
        #print "group size for non snp", len(group1), len(group2)
        #print group1, group2
        if flag == False:
            continue
        #leftKmer = uniqKmer
        if (small1, small2) not in pairSet:
            usedKmers.extend( group1 )
            usedKmers.extend( group2 )
            pairSet.add( (small1, small2) )
            nonPair.append( (small1, small2, ek1, ek2, covL1, covL2, covR1, covR2) )
    print ( "for non snp, used kmer number", len(usedKmers) )
    return nonPair, uniqKmers-set(usedKmers)


def find_non_pair_kmer(uniqKmer, kmerCov, k, left_index, right_index):

    mid = int(k/2)
    left = {}
    for kmer in uniqKmer:
        leftKey = kmer[:mid] 
        if leftKey not in left:
            left[ leftKey ] = []  
        left[leftKey].append( (kmer) )
        Rkmer = tools.reverse(kmer)
        leftKey = Rkmer[:mid] 
        if leftKey not in left:
            left[ leftKey ] = []  
        left[leftKey].append( (Rkmer) )
    
    #build map: overlap is key
    print ("left size", len(left) )
    mapMerge, highRepeat = build_map_merge(left, k)
    print ("map Merge size", len(mapMerge) )
    nonPair, leftKmer = merge_pair(mapMerge, highRepeat, k, left_index, right_index, uniqKmer, kmerCov) 

    print ("non pair size", len(nonPair) )
    fout = open("non_pair", "w")
    sortedNon = sorted(nonPair)
    print ("non pair size", len(sortedNon) )
    for (k1, k2, ek1, ek2, c1, c2, c3, c4) in sortedNon:     
        fout.write("%s %s %s %s %s %s %s %s\n" % (k1,k2,ek1,ek2,c1,c2,c3,c4) )
    fout.close()        
    
    return nonPair, highRepeat, leftKmer



def extend_to_left(h1, left_index, right_index, k, group):
    
    mid = int(k/2)
    key = h1[ : (k-1)]
    Rkey = tools.reverse(key)
    temp, Rtemp = h1, tools.reverse(h1)
    add, Radd = "", ""
    #group.add(h1)
    for i in range(0, mid):
        flag, flagR = False, False
        if key in right_index :
            temp = right_index[key][0] + temp
            Rtemp = tools.reverse(temp)
            flag = True

        if Rkey in left_index:
            Rtemp = Rtemp + left_index[Rkey][0]
            temp = tools.reverse(Rtemp)
            flagR = True

        if flag == True and flagR == False:
            add = right_index[key][0] + add
            Radd = tools.reverse(add) 
            key = temp[: (k-1) ]
            group.append(temp[:k])
            Rkey = tools.reverse(key)
        elif flag == False and flagR == True: 
            Radd = Radd + left_index[Rkey][0]
            add = tools.reverse(Radd)
            Rkey = Rtemp[-(k-1):]
            group.append(Rtemp[-k:])
            key = tools.reverse(Rkey)
        elif flag == flagR:
            if flag == True:
                temp = temp[1:]
                Rtemp = Rtemp[:-1]
            break
    return temp, add

def extend_to_right(h1, left_index, right_index, k, group):

    mid = int(k/2)
    key = h1[-(k-1):]
    Rkey = tools.reverse(key)
    temp, Rtemp = h1, tools.reverse(h1)
    add, Radd = "", ""
    for i in range(0, mid):
        flag, flagR = False, False
        if key in left_index:
            temp = temp + left_index[key][0]
            Rtemp = tools.reverse(temp)
            flag = True

        if Rkey in right_index:
            Rtemp = right_index[Rkey][0] + Rtemp
            temp = tools.reverse(Rtemp)
            flagR = True

        if flag == True and flagR == False:
            add = add + left_index[key][0]
            Radd = tools.reverse(add)
            key = temp[-(k-1):]
            group.append(temp[-k:])
            Rkey = tools.reverse(key)
        elif flag == False and flagR == True: 
            Radd = right_index[Rkey][0] + Radd
            add = tools.reverse(Radd)
            Rkey = Rtemp[:(k-1)]
            group.append(Rtemp[:k])
            key = tools.reverse(Rkey)
        elif flag == flagR:
            if flag == True:
                temp = temp[:-1]
                Rtemp = Rtemp[1:]
            break
    return temp, add
   
def extend_one_pair(h1,h2,left_index, right_index,k, threshold):

    flag = True
    group1, group2 = [], []
    if len(h1) == k: 
        group1.append(h1)
    if len(h2) == k:
        group2.append(h2)
    if len(h1)>k:       
        lenh1 = len(h1)
        for i in range(lenh1-k+1):
            k1temp =h1[i:i+k]
            group1.append(min(k1temp, tools.reverse(k1temp) ) )
    if len(h2)>k:
        lenh2 = len(h2)
        for i in range(lenh2-k+1):
            k2temp = h2[i:i+k]
            group2.append(min(k2temp, tools.reverse(k2temp) ) )

    temp1, add1 = extend_to_left(h1, left_index, right_index, k, group1)
    temp2, add2 = extend_to_left(h2, left_index, right_index, k, group2)
    minL = min (len(add1), len(add2) )
    # some TP add content shift one position equal 
    if minL!= 0 and add1[-minL:] != add2[-minL:]:
        flag = False    
    if minL == int(k/2) and tools.hamming_distance(add1, add2) == 1:
    #if minL == int(k/2) and tools.min_edit_distance(add1, add2) <= 2:
        flag = True    
    ekmer1, add1R = extend_to_right(temp1, left_index, right_index, k, group1)
    ekmer2, add2R = extend_to_right(temp2, left_index, right_index, k, group2) 
    minR = min (len(add1R), len(add2R) )
    if minR!=0 and add1R[0:minR] != add2R[0:minR]:
        flag = False
    # the distance of two snps larger than k/2 smaller than k    
    if minR == int(k/2) and tools.hamming_distance(add1R, add2R) == 1:
    #if minR == int(k/2) and tools.min_edit_distance(add1R, add2R) <= 2:
        flag = True
    if max(minL, minR) <= 2: 
        flag = False
    if min(minL, minR) <= threshold:
        flag = False
    return ekmer1, ekmer2, flag, group1, group2

def extend_pair(hetePairs, left_index, right_index, k):
    extendPair, extendSet = [], set()
    for (h1, h2, ID) in hetePairs:
        ekmer1, ekmer2 = extend_one_pair(h1, h2, left_index, right_index,k)
        small1, small2 = tools.get_smaller_pair_kmer(ekmer1, ekmer2)
        if (small1, small2) not in extendSet:
            extendSet.add( (small1, small2) )
            extendPair.append( (small1, small2, ID) )
    return extendPair   

def check_unique_next(kmer, left_index, right_index):
    Rkmer = tools.reverse(kmer)
    if Rkmer < kmer:
        kmer = Rkmer
    leftKey = kmer[1:]
    rightKey = kmer[:-1]
    #if leftKey in left_index and rightKey in right_index:

    if (leftKey in left_index) or (rightKey in right_index):
        #assert len(left_index[leftKey]) == 1
        #assert len(right_index[rightKey]) == 1
        return True

    return False
''' 
