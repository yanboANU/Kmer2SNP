#########################################################################
# File Name: compareRealRair.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 11:51:38 AEST
#########################################################################
#!/bin/bash
import sys

def read_kmer(filename):
    pair_kmer = []
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            pair_kmer.append(words)
    return pair_kmer      


if len(sys.argv) < 2:
    print ("hete/homo/repeat.31mer uniq.31mer")
    sys.exit()

realFile=sys.argv[1]
#uniqFile="/home/yulin/bio/VariationCalling/experiment/Simulate/30x_30x/kmercalling/K31/sorted_uniq_kmer"
uniqFile=sys.argv[2]

#print ("input %s %s" % (realFile, pairFile) )
list2 = read_kmer(realFile)
list4 = read_kmer(uniqFile)

intersection=[]
index2, index4=0,0
len2, len4 = len(list2), len(list4)
foutFP=open("FPPair", "w") #not real SNP, but found
foutTN=open("TNPair", "w") #real SNP, not found
foutTP=open("TPPair", "w") #real SNP, not found


while index2 < len2 and index4 < len4:
    if (list2[index2][0:1] < list4[index4][0:1]):
        foutTN.write( "%s\n" % (list2[index2][0]) )
        index2 +=1
    elif list2[index2][0:1] > list4[index4][0:1]:   
        #foutFP.write( "%s %s %s %s\n" % ( list4[index4][0], list4[index4][1], list4[index4][2], list4[index4][3] ) )
        #foutFP.write( "%s %s " % ( list4[index4][0], list4[index4][1] ) )
        #if len(list4[index4]) >= 2:
            #foutFP.write( "%s %s \n" % ( list4[index4][2], list4[index4][3] ) )
        #else:
            #foutFP.write( "\n" )
        for ele in list4[index4]:
            foutFP.write( "%s " % ( ele ) )    
        foutFP.write( "\n" )
        index4 +=1
    elif list2[index2][0:1] == list4[index4][0:1]:

        #foutTP.write( "%s %s " % ( list4[index4][0], list4[index4][1] ) )
        #if len(list4[index4]) >= 4:
            #foutTP.write( "%s %s \n" % ( list4[index4][2], list4[index4][3] ) )
        #else:
            #foutTP.write( "\n" )
        for ele in list4[index4]:
            foutTP.write( "%s " % ( ele ) )    
        foutTP.write( "\n" )

        intersection.append( list4[index4] )
        index2 +=1
        index4 +=1
    else: 
        print "unknow"
        print list2[index2]
        print list4[index4]
        index2 +=1
        index4 +=1

foutTP.close()        
foutFP.close()
foutTN.close()
#print "real pair kmer number, find pair kmer number, intersection pair kmer number"
len_inter = len(intersection)   
#print len2, len4, len_inter
#print "TP FP SEN(%) PREC(%)"
#print ( " %s & %s & %.2f & %.2f  " % (len_inter, len4 - len_inter, len_inter*100.0/len2,  len_inter*100.0/len4) )

print ( "%s\t%s\t%s\t%s\t%.2f\t%.2f  " % (len2, len4, len_inter, len4 - len_inter, len_inter*100.0/len2,  len_inter*100.0/len4) )
        

