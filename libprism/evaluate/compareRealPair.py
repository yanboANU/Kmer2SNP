#########################################################################
# File Name: compareRealRair.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 11:51:38 AEST
#########################################################################
#!/bin/bash

import sys
#from libprism.local 
import tools

def read_kmer(filename):
    pair_kmer = []
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            pair_kmer.append(words)
    return sorted(pair_kmer)


def compare_2_list(list2, list4):

    intersection=[]
    index2, index4=0,0
    len2, len4 = len(list2), len(list4)
    foutFP=open("Kmer2SNP-specific", "w") #not real SNP, but found
    foutTN=open("ref-specific", "w") #real SNP, not found
    foutTP=open("shared-snp-pair", "w") #real SNP, not found

    while index2 < len2 and index4 < len4:
        if (list2[index2][0:2] < list4[index4][0:2]):
            for ele in list2[index2]:
                foutTN.write( "%s " % ( ele ) )    
            foutTN.write( "\n" )
            index2 +=1
        elif list2[index2][0:2] > list4[index4][0:2]:   
            for ele in list4[index4]:
                foutFP.write( "%s " % ( ele ) )    
            foutFP.write( "\n" )
            index4 +=1
        elif list2[index2][0:2] == list4[index4][0:2]:
            for ele in list4[index4]:
                foutTP.write( "%s " % ( ele ) )    
            foutTP.write( "\n" )
            intersection.append( list4[index4] )
            index2 +=1
            index4 +=1
        else: 
            print ("unknow error")
            sys.exit()

    foutTP.close()        
    foutFP.close()
    foutTN.close()
    return intersection, len2, len4


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print ("real.31mer find.31mer")
        sys.exit()

    realFile=sys.argv[1]
    pairFile=sys.argv[2]

    #print ("input %s %s" % (realFile, pairFile) )
    list2 = read_kmer(realFile)
    list4 = read_kmer(pairFile)

    intersection, len2, len4 = compare_2_list(list2, list4) 

    len_inter = len(intersection)   
    if len4 != 0:

        print ("#(real kmer pair), #(find kmer pair), #(share kmer pair), recall, precision")
        print ( "%s\t%s\t%s\t%.2f\t%.2f  " % (len2, len4, len_inter, len_inter*100.0/len2,  len_inter*100.0/len4) )
    else:
        print ("#(real kmer pair), #(find kmer pair), #(share kmer pair)")
        print ( "%s\t%s\t%s  " % (len2, len4, len_inter) )

