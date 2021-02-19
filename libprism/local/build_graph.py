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


def snp_edges(m_index, k, left_index, right_index, extendKmers):   

    edges = []
    print ("total number possible pair kmer", len(m_index))   
    count1 = 0
    mid = int(k/2)
    #fout=open("snp_edges", "w")
    for key in iter( m_index ):
        mKeyLen = len(m_index[key])
        if mKeyLen == 1:
            count1 += 1
            continue
        for i in range(mKeyLen-1):
            for j in range(i+1, mKeyLen):
                k1, k2 = m_index[key][i], m_index[key][j]
                ek1, ek2, supportPair, flag = extend_one_pair(k1, k2, left_index, right_index, k)
                if flag: 
                    if k1 < k2:
                        edges.append( (k1, k2, supportPair) )
                        extendKmers[ (k1,k2) ] = (ek1, ek2)              
                        #fout.write("%s %s %s\n" % (tools.transfer_int_kmer(k1, k), tools.transfer_int_kmer(k2, k), supportPair) )
                    else:
                        edges.append( (k2, k1, supportPair) )
                        extendKmers[ (k2, k1) ] = (ek2, ek1)              
                        #fout.write("%s %s %s\n" % (tools.transfer_int_kmer(k2, k), tools.transfer_int_kmer(k1, k), supportPair) )
    print ("kmer cannot find pair number", count1)      
    print ("snp edges number", len(edges))      
    #fout.close()
    return edges 
    

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

def update_map_merge(k1, k2, key1, key2, mapMerge):
    
    minkey1, minkey2 = tools.get_smaller_pair_binary(key1, key2)
    if (minkey1, minkey2) not in mapMerge:
        mapMerge[ (minkey1, minkey2) ] = []   
    Rk1 = tools.reverse_binary_string(k1)
    Rk2 = tools.reverse_binary_string(k2)
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
    for key in iter(left):
        groupSize = len(left[key])
        for i in range(0, groupSize-1):
            for j in range(i+1, groupSize):
                k1, k2 = left[key][i], left[key][j]
                binaryk1 = bin(k1)[2:].zfill(2*k)
                binaryk2 = bin(k2)[2:].zfill(2*k) 
                if binaryk1[k-1:k+1] == binaryk2[k-1:k+1]: #or (k1 in mappedKmer) or (k2 in mappedKmer):
                    continue
                dis = 0
                cnt = k + 1
                diffPos = 0
                while cnt < 2*k:
                    if binaryk1[cnt:cnt+2] != binaryk2[cnt:cnt+2]:
                        dis += 1
                        diffPos = cnt
                    cnt += 2    
                    if dis >= 2:
                        break
                if cnt == 2*k and dis == 1:
                    key1, key2 = binaryk1[ diffPos - (k-1) : ] , binaryk2[ diffPos - (k-1) : ] 
                    mink1, mink2 = tools.get_smaller_pair_int(k1, k2, k)
                    #assert (mink1 == k1 and mink2==k2) or (mink1 == k2 and mink2 == k1) --> False
                    if mink1 not in hisMap:
                        hisMap[mink1] = 0
                    if mink2 not in hisMap:
                        hisMap[mink2] = 0
                    hisMap[mink1] += 1
                    hisMap[mink2] += 1
                    update_map_merge(binaryk1, binaryk2, key1, key2, mapMerge) # ????
                    #mappedKmer.add(k1)
                    #mappedKmer.add(k2)
                    #break #3 lines add 22 Aug. a kmer only allow one kmer hamming distance equal to 2
    
    left.clear()
    highRepeat = set()
    for key in hisMap:
        if hisMap[key] > 5: # rule 
            highRepeat.add(key)
    print ("high Repeat kmer number", len(highRepeat) )
    return mapMerge, highRepeat

def merge_pair(mapMerge, highRepeat, k, left_index, right_index, kmerCov, extendKmers):
    edges = []

    #fout = open("non_snp_edges" ,"w")
    for (key1, key2) in iter(mapMerge):
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
                Lmin1 = min(Left1, tools.reverse_binary_string(Left1))
                Lmin2 = min(Left2, tools.reverse_binary_string(Left2))
                Rmin1 = min(Right1, tools.reverse_binary_string(Right1))
                Rmin2 = min(Right2, tools.reverse_binary_string(Right2))
                Lint1, Lint2 = int(Lmin1, 2), int(Lmin2,2)
                Rint1, Rint2 = int(Rmin1, 2), int(Rmin2,2 )
                
                covL1, covL2 = kmerCov[Lint1], kmerCov[Lint2]
                covR1, covR2 = kmerCov[Rint1], kmerCov[Rint2]
                 
                # rule  
                if (Lint1 in highRepeat or Lint2 in highRepeat or
                        Rint1 in highRepeat or Rint2 in highRepeat):
                    continue
                
                # rule 
                if covL1 > 1.5*covR1 or 1.5*covL1 < covR1:
                    continue
                if covL2 > 1.5*covR2 or 1.5*covL2 < covR2:
                    continue

                if Right1[-l:] == Left1[:l] and Right2[-l:] == Left2[:l]:
                    merge1 = Right1 + Left1[l:] # this is binary string
                    merge2 = Right2 + Left2[l:]
                    # rule binary string should be 2, kmer string should be 1
                    if merge1[2:] == merge2[:-2] or merge1[:-2]==merge2[2:]: # debug 1->2 
                        print (merge1, merge2)
                        continue
                    ek1L, ek2L, supportPairL, flagL = extend_one_pair(Lint1, Lint2, left_index, right_index, k)
                    ek1R, ek2R, supportPairR, flagR = extend_one_pair(Rint1, Rint2, left_index, right_index, k)
                    if flagL and flagR:
                        #p1 = tools.hamming_distance2(ek1L, ek2L)
                        #p2 = tools.hamming_distance2(ek1R, ek2R)
                        #if p1[0] > k/2 and p2[0] > k/2 and len(ek1L)-p1[1] > k/2 and len(ek1R)-p2[1] > k/2:
                        # make sure support pairs between non pair are exists
                        edges.append( (Lint1, Lint2, supportPairL) )
                        edges.append( (Rint1, Rint2, supportPairR) )

                        #fout.write("%s %s %s\n" % (tools.transfer_int_kmer(Lint1, k), tools.transfer_int_kmer(Lint2, k), supportPairL) )
                        #fout.write("%s %s %s\n" % (tools.transfer_int_kmer(Rint1, k), tools.transfer_int_kmer(Rint2, k), supportPairR) )
                    
                elif Left1[-l:] == Right1[:l] and Left2[-l:] == Right2[:l]:
                    merge1 = Left1 + Right1[l:]
                    merge2 = Left2 + Right2[l:]
                    # rule   
                    if merge1[2:] == merge2[:-2] or merge1[:-2]==merge2[2:]: #debug 1->2
                        print (merge1, merge2)
                        continue
                    ek1L, ek2L, supportPairL, flagL = extend_one_pair(Lint1, Lint2, left_index, right_index, k)
                    ek1R, ek2R, supportPairR, flagR = extend_one_pair(Rint1, Rint2, left_index, right_index, k)
                    if flagL and flagR:
                        #p1 = tools.hamming_distance2(ek1L, ek2L)
                        #p2 = tools.hamming_distance2(ek1R, ek2R)
                        #if p1[0] > k/2 and p2[0] > k/2 and len(ek1L)-p1[1] > k/2 and len(ek1R)-p2[1] > k/2:
                        edges.append( (Lint1, Lint2, supportPairL) )
                        edges.append( (Rint1, Rint2, supportPairR) )
                        #fout.write("%s %s %s\n" % (tools.transfer_int_kmer(Lint1, k), tools.transfer_int_kmer(Lint2, k), supportPairL) )
                        #fout.write("%s %s %s\n" % (tools.transfer_int_kmer(Rint1, k), tools.transfer_int_kmer(Rint2, k), supportPairR) )
                else:
                    continue #print ("one side")
                 
                #TTTTTTTTTTTTTTTCAAAAAAAAAAAAAAAA # also useful SNP and indel 
                #TTTTTTTTTTTTTTTTCAAAAAAAAAAAAAAA 
                small1, small2 = tools.get_smaller_pair_binary(merge1, merge2)
                s1 = tools.transfer_binary_string_2_kmer(small1)
                s2 = tools.transfer_binary_string_2_kmer(small2)
                if Lmin1 < Lmin2:
                    extendKmers[(Lint1, Lint2)] = (s1, s2)
                else:
                    extendKmers[(Lint2, Lint1)] = (s1, s2)

                if Rmin1 < Rmin2:
                    extendKmers[(Rint1, Rint2)] = (s1, s2)
                else:
                    extendKmers[(Rint2, Rint1)] = (s1, s2)
    #fout.close()                
    return edges 


def non_snp_edges(kmerCov, k, left_index, right_index, extendKmers):

    left = {}
    for intkmer in iter(kmerCov):
        binarykmer = bin(intkmer)[2:].zfill(k*2) 
        leftKey = int(binarykmer[:k-1], 2)
        if leftKey not in left:
            left[ leftKey ] = []  
        left[leftKey].append( (intkmer) )
        Rkmer = tools.reverse_binary_string(binarykmer)
        leftKey = int(Rkmer[:k-1], 2)
        if leftKey not in left:
            left[ leftKey ] = []  
        left[leftKey].append( (int(Rkmer,2) ) )
         
    print ("left size", len(left) )
    logging.info( "left Memory %s bytes " % (sys.getsizeof(left) ) )
    mapMerge, highRepeat = build_map_merge(left, k)
    logging.info( "mapMerge Memory %s bytes " % (sys.getsizeof(mapMerge) ))
    logging.info( "highRepeat Memory %s bytes " % (sys.getsizeof(highRepeat) ))
    print ("map Merge size", len(mapMerge) )
    edges = merge_pair(mapMerge, highRepeat, k, left_index, right_index, kmerCov, extendKmers)
    print ("non snp edges number", len(edges))       
    return edges


def extend_to_left(h1, left_index, right_index, k, group): 
    
    mid = int((k-1)/2)
    h1_binary = bin(h1)[2:].zfill(2*k) 

    #logging.debug( "in extend to left, len h1_binary %s " % (len(h1_binary) ) )
    temp, Rtemp = h1_binary, tools.reverse_binary_string(h1_binary)

    #logging.debug( "first reverse binary string, len h1_binary %s " % (len(h1_binary) ) )
    key = int(temp[ : (2*k-2)], 2)
    Rkey = int(Rtemp[2:], 2)
    add, Radd = "", ""
    for i in range(0, mid):
        flag, flagR = False, False
        if key in right_index :
            temp = bin(right_index[key])[2:].zfill(2) + temp
            Rtemp = tools.reverse_binary_string(temp)
            flag = True

        if Rkey in left_index:
            Rtemp = Rtemp + bin(left_index[Rkey])[2:].zfill(2)
            temp = tools.reverse_binary_string(Rtemp)
            flagR = True
        
        if  flag == flagR:
            if flag == True:
                temp = temp[2:]
                Rtemp = Rtemp[:-2]
            break

        if flag == True and flagR == False:
            add = bin(right_index[key])[2:].zfill(2) + add
            Radd = tools.reverse_binary_string(add)
            group.append( int(temp[:2*k],2) )
        elif flag == False and flagR == True: 
            Radd = Radd + bin(left_index[Rkey])[2:].zfill(2)
            add = tools.reverse_binary_string(Radd)
            group.append( int(Rtemp[-2*k:], 2) )

        key = int(temp[: (2*k-2) ] , 2)
        Rkey = int(Rtemp[-(2*k-2):], 2) 
    return temp, add

def extend_to_right(h1_binary, left_index, right_index, k, group):

    #logging.info( "in extend to right ")
    mid = int(k/2)
    #key = h1[-(k-1):]
    #Rkey = tools.reverse(key)
    #temp, Rtemp = h1, tools.reverse(h1)
    #logging.debug( "in extend to right, len h1_binary %s " % (len(h1_binary) ) )
    temp, Rtemp = h1_binary, tools.reverse_binary_string(h1_binary) 
    key = int(temp[-(2*k-2):], 2)
    Rkey = int(Rtemp[ : (2*k-2)], 2)
    add, Radd = "", ""
    for i in range(0, mid):
        flag, flagR = False, False
        if key in left_index:
            temp = temp + bin(left_index[key])[2:].zfill(2)
            Rtemp = tools.reverse_binary_string(temp)
            flag = True

        if Rkey in right_index:
            Rtemp = bin(right_index[Rkey])[2:].zfill(2) + Rtemp
            temp = tools.reverse_binary_string(Rtemp)
            flagR = True

        elif flag == flagR:
            if flag == True:
                temp = temp[:-2]
                Rtemp = Rtemp[2:]
            break

        if flag == True and flagR == False:
            add = add + bin(left_index[key])[2:].zfill(2)
            Radd = tools.reverse_binary_string(add)
            group.append( int(temp[-2*k:], 2) )
        elif flag == False and flagR == True: 
            Radd = bin(right_index[Rkey])[2:].zfill(2) + Radd
            add = tools.reverse_binary_string(Radd)
            group.append( int(Rtemp[:2*k], 2) )
    
        key = int(temp[-(2*k-2) : ] , 2)
        Rkey = int(Rtemp[ : (2*k-2)], 2)  
    return temp, add

def init_groups(h1,h2,k):
   
    group1, group2 = [], []
    group1.append(h1)
    group2.append(h2)
    '''
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
    '''
    return group1, group2
   
def extend_one_pair(h1, h2, left_index, right_index, k):

    flag = True
    supportPairL, supportPairR = 0, 0
    group1, group2 = init_groups(h1, h2, k)
    temp1, add1 = extend_to_left(h1, left_index, right_index, k, group1)
    temp2, add2 = extend_to_left(h2, left_index, right_index, k, group2)
    minL = min ( len(add1), len(add2 ) )
    
    #assert (minL <= 30)
    
    add1, add2 = add1[-minL:], add2[-minL:] 
    for i in range(minL, 0, -2):
        if add1[i-2 : i] == add2[i-2 : i]:
            supportPairL = supportPairL + 1
        else:
            break
    #assert supportPairL <= 15        
    
    temp1, temp2 = temp1[-(supportPairL+k)*2:], temp2[-(supportPairL+k)*2:]

    #print ("debug supportPairL, temp1, temp2", supportPairL, temp1, temp2)
    
    ekmer1, add1R = extend_to_right(temp1, left_index, right_index, k, group1)
    ekmer2, add2R = extend_to_right(temp2, left_index, right_index, k, group2) 
    minR = min ( len(add1R), len(add2R) )

    #assert minR <= 30 
    for i in range(0, minR, 2):
        if add1R[i:i+2] == add2R[i:i+2]:
            supportPairR = supportPairR + 1
        else:
            break
    #assert supportPairR <= 15
    
    ekmer1, ekmer2 = ekmer1[:(supportPairL+k+supportPairR)*2], ekmer2[:(supportPairL+k+supportPairR)*2]

    #print ("debug binary string lens", len(ekmer1), len(ekmer2))
    ekmer1 = tools.transfer_binary_string_2_kmer(ekmer1)
    ekmer2 = tools.transfer_binary_string_2_kmer(ekmer2)

    #if max(minL, minR) <= 2: 
    #    flag = False
    #if min(minL, minR) <= threshold:
        #flag = False

    if supportPairL > int(k/2):
        print ("supportPairL", supportPairL, supportPairR)
        #sys.exit()
     
    if supportPairR > int(k/2):
        print ("supportPairR", supportPairL, supportPairR)
        #sys.exit()
    
    #print ("debug", len(ekmer1), len(ekmer2), supportPairL+supportPairR, flag)
    return ekmer1, ekmer2, supportPairL+supportPairR, flag

'''
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
            uniqKmer[kmer] = coverage
    return uniqKmer
'''
