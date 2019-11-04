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
#from tools import *

def pick_smaller_unique_kmer(input_filename, low, high):
    uniqKmer = {}
    #fout = open("uniq_kmer" ,"w") 
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
            uniqKmer[kmer] = coverage
            #fout.write("%s %s\n" % (kmer, coverage))
    return uniqKmer

# hamming distance = 1
def snp_edges(uniqKmers, k, left_index, right_index, extendKmers, slabel):   

    edges, m = [], {}
    mid = int(k/2)
    print ("unique kmer number", len(uniqKmers))        
    for kmer in uniqKmers:
        key= kmer[:mid] + kmer[mid+1:]
        if key not in m:
            m[key] = []
        m[key].append( kmer )
    print ("total number possible pair kmer", len(m))    
    count1 = 0
    fout = open('snp_edges', "w")
    for key in m:
        mKeyLen = len(m[key])
        if mKeyLen == 1:
            count1 += 1
        for i in range(mKeyLen-1):
            for j in range(i+1, mKeyLen):
                k1, k2 = m[key][i], m[key][j]
                '''
                cov1, cov2 = uniqKmers[k1], uniqKmers[k2]
                if k1[1:] == k2[:-1] or k1[:-1]==k2[1:]:
                    print (k1, k2)
                    continue
                '''    
                ek1, ek2, supportPair, flag = extend_one_pair(k1, k2, left_index, right_index, k)
                if flag:
                    edges.append( (k1, k2, supportPair) )
                    extendKmers[ (k1,k2) ] = (ek1, ek2)              
                    
                    if k1 < k2:
                        fout.write("%s %s %s\n" % (k1, k2, supportPair) )
                    else:
                        fout.write("%s %s %s\n" % (k2, k1, supportPair) )
                       
                else:
                    if k1 < k2:
                        fout.write("%s %s %s\n" % (k1, k2, 0) )
                    else:
                        fout.write("%s %s %s\n" % (k2, k1, 0) )
                        
    fout.close()                
    print ("kmer cannot find pair number", count1)      
    print ("snp edges number", len(edges))       
    return edges 
    

'''
print ("before filter mapK")
uniqMapK = {}
for key in mapK:
    for (kmer, cov) in mapK[key]:
        leftHalf = kmer[:mid]
        rightHalf = kmer[mid+1:]
        if (tools.hamming_distance(rightHalf, kmer[mid : -1 ] ) <= 1 or # mutation => delete 
                tools.hamming_distance(leftHalf, kmer[1:mid+1]) <= 1 ):
            mapK[key].remove((kmer, cov))
    if len(mapK[key]) == 1:
        uniqMapK[key] = mapK[key]

print ("after filter mapK")
'''
def indel_edges(kmerCov, uniqK_1mers, k, left_index, right_index):  
    mid = int(k/2)
    mapK = {}
    edges = []
    print ( "uniq K-1 mer size", len(uniqK_1mers) )
    for temp in kmerCov:
        leftHalf = temp[:mid]
        rightHalf = temp[mid+1:]
        if leftHalf.count(leftHalf[0]) >= mid-1 or rightHalf.count(rightHalf[0]) >= mid-1:
            print (temp)
            continue
        key = leftHalf + rightHalf
        if key not in mapK:
            mapK[key] = []
        mapK[key].append( temp )
    for key in uniqK_1mers:
        if key in mapK and len( mapK[key] ) >= 1:
            for i in range( len(mapK[key]) ):
                kmer = mapK[key][i]
                ek1, ek2, supportPair, flag = extend_one_pair(kmer, key, left_index, right_index, k)
                if flag:
                    edges.append( (kmer, key, supportPair*1) )
     
    print ("indel edges number", len(edges))       
    return edges

def update_map_merge(left_i, left_j, key1, key2, mapMerge):
    
    k1 = left_i
    k2 = left_j
    minkey1, minkey2 = tools.get_smaller_pair_kmer(key1, key2)
    if (minkey1, minkey2) not in mapMerge:
        mapMerge[ (minkey1, minkey2) ] = []   
    Rk1, Rk2 = tools.reverse(k1), tools.reverse(k2)
    if k1.count(minkey1) == 1 and k2.count(minkey2) == 1:
        mapMerge[ (minkey1, minkey2) ].append( (k1, k2) )
    elif k1.count(minkey2) == 1 and k2.count(minkey1) == 1: 
        mapMerge[ (minkey1, minkey2) ].append( (k2, k1) )
    elif Rk1.count(minkey1) == 1 and Rk2.count(minkey2) == 1:   
        mapMerge[ (minkey1, minkey2) ].append( (Rk1, Rk2) )
    elif Rk1.count(minkey2) == 1 and Rk2.count(minkey1) == 1:
        mapMerge[ (minkey1, minkey2) ].append( (Rk2, Rk1) )
    else:
        print ("something wrong 1")
        sys.exit() 
    return     

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
        if hisMap[key] > 5: # rule 
            highRepeat.add(key)
    
    print ("high Repeat kmer number", len(highRepeat) )
    return mapMerge, highRepeat

def merge_pair(mapMerge, highRepeat, k, left_index, right_index, kmerCov, extendKmers):
    edges = []
    for (key1, key2) in mapMerge:
        mlen = len( mapMerge [ (key1, key2) ] )
        #if mlen > 2 or mlen==1: # one key only allow a pair
        #    continue
        if mlen == 1:
            continue
        for i in range(mlen-1):
            for j in range(i+1, mlen):
                Left1, Left2 = mapMerge[ (key1, key2) ][i]
                Right1, Right2 = mapMerge[ (key1, key2) ][j]
                l = len(key1)
                Lmin1 = min(Left1, tools.reverse(Left1))
                Lmin2 = min(Left2, tools.reverse(Left2))
                Rmin1 = min(Right1, tools.reverse(Right1))
                Rmin2 = min(Right2, tools.reverse(Right2))
                covL1, covL2 = kmerCov[Lmin1], kmerCov[Lmin2]
                covR1, covR2 = kmerCov[Rmin1], kmerCov[Rmin2]
                 
                # rule  
                if (Lmin1 in highRepeat or Lmin2 in highRepeat or
                        Rmin1 in highRepeat or Rmin2 in highRepeat):
                    continue
                
                # rule
                if covL1 > 1.5*covR1 or 1.5*covL1 < covR1:
                    continue
                if covL2 > 1.5*covR2 or 1.5*covL2 < covR2:
                    continue
         
                if Right1[-l:] == Left1[:l] and Right2[-l:] == Left2[:l]:
                    merge1 = Right1 + Left1[l:]
                    merge2 = Right2 + Left2[l:]
                    # rule
                    if merge1[1:] == merge2[:-1] or merge1[:-1]==merge2[1:]:
                        print (merge1, merge2)
                        continue
                    ek1L, ek2L, supportPairL, flagL = extend_one_pair(Lmin1, Lmin2, left_index, right_index, k)
                    ek1R, ek2R, supportPairR, flagR = extend_one_pair(Rmin1, Rmin2, left_index, right_index, k)
                    if flagL and flagR:
                        p1 = tools.hamming_distance2(ek1L, ek2L)
                        p2 = tools.hamming_distance2(ek1R, ek2R)
                        #print (Lmin1, Lmin2, ek1L, ek2L, supportPairL, p1)
                        #print (Rmin1, Rmin2, ek1R, ek2R, supportPairR, p2)
                        #if p1[0] > k/2 and p2[0] > k/2 and len(ek1L)-p1[1] > k/2 and len(ek1R)-p2[1] > k/2:
                        # make sure support pairs between non pair are exists
                        edges.append( (Lmin1, Lmin2, supportPairL) )
                        edges.append( (Rmin1, Rmin2, supportPairR) )
                    
                elif Left1[-l:] == Right1[:l] and Left2[-l:] == Right2[:l]:
                    merge1 = Left1 + Right1[l:]
                    merge2 = Left2 + Right2[l:]
                    # rule   
                    if merge1[1:] == merge2[:-1] or merge1[:-1]==merge2[1:]:
                        print (merge1, merge2)
                        continue
                    ek1L, ek2L, supportPairL, flagL = extend_one_pair(Lmin1, Lmin2, left_index, right_index, k)
                    ek1R, ek2R, supportPairR, flagR = extend_one_pair(Rmin1, Rmin2, left_index, right_index, k)
                    if flagL and flagR:
                        p1 = tools.hamming_distance2(ek1L, ek2L)
                        p2 = tools.hamming_distance2(ek1R, ek2R)
                        #print (Lmin1, Lmin2, ek1L, ek2L, supportPairL, p1)
                        #print (Rmin1, Rmin2, ek1R, ek2R, supportPairR, p2)
                        #if p1[0] > k/2 and p2[0] > k/2 and len(ek1L)-p1[1] > k/2 and len(ek1R)-p2[1] > k/2:
                        edges.append( (Lmin1, Lmin2, supportPairL) )
                        edges.append( (Rmin1, Rmin2, supportPairR) ) 
                else:
                    continue #print ("one side")
                 
                #TTTTTTTTTTTTTTTCAAAAAAAAAAAAAAAA # also useful SNP and indel 
                #TTTTTTTTTTTTTTTTCAAAAAAAAAAAAAAA 
                small1, small2 = tools.get_smaller_pair_kmer(merge1, merge2)
                if Lmin1 < Lmin2:
                    extendKmers[(Lmin1, Lmin2)] = (small1, small2)
                else:
                    extendKmers[(Lmin2, Lmin1)] = (small1, small2)

                if Rmin1 < Rmin2:
                    extendKmers[(Rmin1, Rmin2)] = (small1, small2)
                else:
                    extendKmers[(Rmin2, Rmin1)] = (small1, small2)
    return edges 


def non_snp_edges(kmerCov, k, left_index, right_index, extendKmers, slabel):

    mid = int(k/2)
    left = {}
    for kmer in kmerCov:
        leftKey = kmer[:mid] 
        if leftKey not in left:
            left[ leftKey ] = []  
        left[leftKey].append( (kmer) )
        Rkmer = tools.reverse(kmer)
        leftKey = Rkmer[:mid] 
        if leftKey not in left:
            left[ leftKey ] = []  
        left[leftKey].append( (Rkmer) )
    
    print ("left size", len(left) )
    mapMerge, highRepeat = build_map_merge(left, k)
    print ("map Merge size", len(mapMerge) )
    edges = merge_pair(mapMerge, highRepeat, k, left_index, right_index, kmerCov, extendKmers) 
    '''
    print ("non pair size", len(nonPair) )
    fout = open("non_pair", "w")
    sortedNon = sorted(nonPair)
    print ("non pair size", len(sortedNon) )
    for (k1, k2, ek1, ek2, c1, c2, c3, c4) in sortedNon:     
        fout.write("%s %s %s %s %s %s %s %s\n" % (k1,k2,ek1,ek2,c1,c2,c3,c4) )
    fout.close()        
    '''
    print ("non snp edges number", len(edges))       
    return edges




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

def init_groups(h1,h2,k):

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

    return group1, group2
   
def extend_one_pair(h1, h2, left_index, right_index, k):

    flag = True
    supportPairL, supportPairR = 0, 0
    group1, group2 = init_groups(h1, h2, k)
    temp1, add1 = extend_to_left(h1, left_index, right_index, k, group1)
    temp2, add2 = extend_to_left(h2, left_index, right_index, k, group2)
    minL = min ( len(add1), len(add2) )
    for i in range(1, minL+1):
        if add1[0-i] == add2[0-i]:
            supportPairL = supportPairL + 1
        else:
            break

    ''' 
    if minL!= 0 and add1[-minL:] != add2[-minL:]: #rule exist but conflict
        flag = False    
        if minL == int(k/2) and tools.hamming_distance(add1, add2) == 1:
            flag = True    
    '''
    minTemp = min( len(temp1), len(temp2)  )
    temp1, temp2 = temp1[-minTemp:], temp2[-minTemp:]
    
    ekmer1, add1R = extend_to_right(temp1, left_index, right_index, k, group1)
    ekmer2, add2R = extend_to_right(temp2, left_index, right_index, k, group2) 
    ''' 
    if flag == False:
        return ekmer1, ekmer2, 0, flag
    '''
    minR = min ( len(add1R), len(add2R) )
    for i in range(minR):
        if add1R[i] == add2R[i]:
            supportPairR = supportPairR + 1
        else:
            break
    '''    
    if minR!=0 and add1R[0:minR] != add2R[0:minR]: # rule exist but conflict
        flag = False
        # the distance of two snps larger than k/2 smaller than k    
        if minR == int(k/2) and tools.hamming_distance(add1R, add2R) == 1:
            flag = True
    '''
    minTemp = min( len(ekmer1), len(ekmer2)  )
    ekmer1, ekmer2 = ekmer1[:minTemp], ekmer2[:minTemp]   
    '''
    if supportPairL < 1 or supportPairR < 1: # rule, at least one side one support
        flag = False
    '''    
    #if max(minL, minR) <= 2: 
    #    flag = False
    #if min(minL, minR) <= threshold:
        #flag = False
     
    return ekmer1, ekmer2, min(supportPairL, supportPairR), flag

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
