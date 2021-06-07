#########################################################################
# File Name: kmer2snp.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Fri 04 Jun 2021 09:28:10 PM AEST
#########################################################################
#!/bin/bash
import os
import sys
import argparse
import logging
import time
import networkx as nx

from libprism.local import build_graph
from libprism.local import kmercalling
from libprism.local import tools


parser = argparse.ArgumentParser(description="Variation calling for Single-individual based on NGS reads")
parser.add_argument('--k', help='kmer size', required=True)
parser.add_argument('--c', help='homozygous coverage', required=True)
parser.add_argument('--fastaq', help='NGS reads in fasta or fastq format, either gziped or not. To pass several files as input, separate file names by a comma (",") ', required=True)
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
logging.basicConfig(stream=sys.stderr, level=logging.INFO)
#############################################################################

parser.add_argument('--b', help='b=0 only call isolated SNV (only one snv in a k-mer), b=1 call at most 2 SNVs in a k-mer  (Default b = 1)', required=False)
parser.add_argument('--t1', help='kmer freq txt file (Default is DSK result)', required=False)
parser.add_argument('--c1', help='low coverage (Default is FindGSE result)', required=False)
parser.add_argument('--c2', help='high coverage (Default is FindGSE result)', required=False)
parser.add_argument('--r', help='heterozygous rate (Default is FindGSE result)', required=False)
args = parser.parse_args()


if args.t1 == None:
    command="sh /path2Kmer2SNP/script/run_dsk.sh " + args.k + " " + args.c + " " + args.fastaq
    print (command)
    os.system( command )
    args.t1 = "chr_k" + args.k  + ".txt"

if args.c1 == None or args.c2 == None or args.r == None:
    command="sh /path2Kmer2SNP/script/run_findgse.sh " + args.k + " " + args.c + " " + args.fastaq
    print (command)
    os.system( command )
    c1, c2, r = tools.read_findGSE_result("hete.para") 

if args.c1 == None:
    args.c1 = c1

if args.c2 == None:
    args.c2 = c2

if args.r == None:
    args.r = r

if args.b != None:
    b = int(args.b)
else:
    b = 1


k = int(args.k)
c = float(args.c)
lowCov, highCov = int(args.c1), int(args.c2)
heteRate=float(args.r)

print ("Estimated heteroyzgous rate is", heteRate)
print ("heterozygous kmer coverage range", lowCov, highCov)
logging.info("heterozygous kmer coverage range %s %s" % (lowCov, highCov))

time1 = time.clock()
kmerCov = kmercalling.pick_smaller_unique_kmer( args.t1, 
           lowCov, highCov)
time2 = time.clock()

logging.info( "picked heterozygous Kmer number: %s" % len(kmerCov) )
logging.info( "finish reading (kmer cov) file, cost %.2f seconds" % (time2-time1) )

#---------------------------------------------------------------
left_index, right_index, m_index = kmercalling.build_index(kmerCov.keys(), k)
edges = []
extendKmers = {}
edges.extend( build_graph.snp_edges(m_index, k, left_index, right_index, extendKmers) )
m_index.clear()
time3 = time.clock()
logging.info( "add snp edges, cost %.2f second" % (time3 - time2) )
logging.info( "edge number %s" % ( len(edges) )  )

if b == 1:
    edges.extend( build_graph.non_snp_edges(kmerCov, k, left_index, right_index, extendKmers) )
    time4 = time.clock()
    logging.info( "add non snp edges, cost %.2f second" % (time4 - time3) )
    logging.info( "edge number %s" % ( len(edges) )  )
else:
    time4 = time.clock()

kmerCov.clear()
left_index.clear()
right_index.clear()
#----------------------------------------------------------------------
G=nx.Graph()
G.add_weighted_edges_from(edges)

time6 = time.clock()
logging.info( "build graph, cost %.2f seconds" % (time6 - time4) )
print ("number of components: ", nx.number_connected_components(G) )
time7 = time.clock()
logging.info("number of components: %s" % nx.number_connected_components(G) )


graphs = list(nx.connected_component_subgraphs(G)) # <=2.3 networkx
#graphs = list(G.subgraph(c) for c in nx.connected_components(G))
#---------------------------------------------------------------

componentSize = {}
selectedMates = {}
weightDis = {}
for g in graphs:
    mate = nx.max_weight_matching(g) 

    gSize = g.size()
    if (gSize >= 1):
        #print ("connected component size: ", g.size())
        #print ("heterozygous kmer pair number: ", len(mate))
        if gSize not in componentSize:
            componentSize[gSize] = 1
        else:
            componentSize[gSize] += 1

    for (ele1, ele2) in iter(mate):
        w = g.get_edge_data(ele1, ele2)['weight']
        selectedMates[(ele1, ele2)] = w
        if w not in weightDis:
            weightDis[w] = 0            
        weightDis[w] += 1


time8 = time.clock()
logging.info( "compute max weight matching, cost %.2f seconds" % (time8 - time7) )

#sortedWeightDis = sorted(weightDis.items())
#print("debug weight distribution:", sortedWeightDis)
#print ( "component size", sorted(componentSize.items()) )

weightThreshold = tools.calc_weightThreshold(heteRate,k)
#print ("weight threshold", weightThreshold)  


#---------------------------------------------------------
marknon = {}
indelPair, snpPair, nonPair = [], [], []
mateKeys = selectedMates.keys()

snpExtendFile = "k_" + str(k) + "_pair.snp.extend"   # pick back extended kmer and those kmer pair
snpExtendOut = open(snpExtendFile, "w")

logging.info( "selected %s kmer pairs" % len(mateKeys) )
for (e1, e2) in iter(mateKeys):
    w= selectedMates[(e1,e2)] 
    ele1 = tools.transfer_int_kmer(e1, k)
    ele2 = tools.transfer_int_kmer(e2, k)
    if tools.hamming_distance(ele1, ele2)==1:
        if w < weightThreshold:
            continue
        if ele1 <= ele2:
            snpPair.append( ( ele1, ele2, w ) )
            ek1, ek2 = extendKmers[ (e1, e2) ]
        else:
            ek1, ek2 = extendKmers[ (e2, e1) ]
            snpPair.append( (ele2, ele1, w ) )
        snpExtendOut.write("%s %s %s\n" % (ek1, ek2, w) )
    if tools.hamming_distance(ele1, ele2)==2:
        if w < 2*weightThreshold:  
            continue
        if ele1 <= ele2:
            temp = extendKmers[(e1, e2)]
        else:
            temp = extendKmers[(e2, e1)]
        if temp not in marknon:
            nonPair.append( (temp[0], temp[1], w ) )
            marknon[temp] = 1
        else:
            marknon[temp] += 1

snpExtendOut.close()
#------------------------------------------------------------------------


sSnpPair, sNonPair = sorted(snpPair), sorted(nonPair)

snpFile = "k_" + str(k) + "_pair.snp"
snpOut = open(snpFile, "w")
#for (ele1, ele2, w) in iter(sSnpPair):
for (ele1, ele2, w) in iter(snpPair): # _pair.snp and _pair.snp.extend one correspond one
    snpOut.write("%s %s %s\n" % (ele1, ele2, w) )
snpOut.close()

if b == 1:
    nonFile = "k_" + str(k) + "_pair.non"
    nonOut = open(nonFile, "w")
    nonSepPairKmer = set()
    logging.debug("len nonPair %s" % (len(sNonPair)))
    for (ele1, ele2, w) in iter(sNonPair):
        if marknon[(ele1, ele2)] >= 2: # both selected, then think it's true
            nonOut.write("%s %s %s\n" % (ele1, ele2, w)) 
            h1Pre, h2Pre, h1Suf, h2Suf = ele1[:k], ele2[:k], ele1[-k:], ele2[-k:]
            smallerH1P, smallerH2P = tools.get_smaller_pair_kmer(h1Pre, h2Pre)
            smallerH1S, smallerH2S = tools.get_smaller_pair_kmer(h1Suf, h2Suf)
            c1 = tools.calc_ATCG_zero_num( smallerH1P )
            c2 = tools.calc_ATCG_zero_num( smallerH2P )
            c3 = tools.calc_ATCG_zero_num( smallerH1S )
            c4 = tools.calc_ATCG_zero_num( smallerH2S )
            if c1 == 2 or c2 == 2 or c3 == 2 or c4 == 2:
                continue
            nonSepPairKmer.add( (smallerH1P, smallerH2P, w) )
            nonSepPairKmer.add( (smallerH1S, smallerH2S, w) )
    nonOut.close()


    nonSepFile = "k_" + str(k) + "_pair.non.sep"
    nonSepOut = open(nonSepFile, "w")
    sNon = sorted(list(nonSepPairKmer))
    for (ele1, ele2, w) in iter(sNon):
        nonSepOut.write("%s %s %s\n" % (ele1, ele2, w) )
    nonSepOut.close()

#--------------
