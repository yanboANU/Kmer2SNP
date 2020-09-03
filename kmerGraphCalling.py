#########################################################################
# File Name: variationCalling.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Mon 18 Mar 2019 13:59:49 AEDT
#########################################################################
#!/bin/bash
import argparse
import os
import sys
from libprism.local import build_graph
from libprism.local import kmercalling
from libprism.local import tools
import logging
import time
import networkx as nx
#############################################################################

parser = argparse.ArgumentParser(description="Variation calling for Single-individual based on NGS reads")
parser.add_argument('--t1', help='kmer freq txt file (dsk result)', required=True)
parser.add_argument('--c1', help='low coverage', required=True)
parser.add_argument('--c2', help='high coverage', required=True)
parser.add_argument('--k', help='kmer size', required=True)
parser.add_argument('--r', help='heterozygous rate', required=True)
args = parser.parse_args()
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
logging.basicConfig(stream=sys.stderr, level=logging.INFO)
#############################################################################
time1 = time.clock()
lowCov, highCov = int(args.c1), int(args.c2)
k = int(args.k)
heteRate=float(args.r)
#extendLen = 1 
extendLen = tools.calc_extendLen((lowCov+highCov)/2, heteRate) 
lowCov = lowCov - extendLen
highCov = highCov + extendLen

print ("Estimated heteroyzgous rate is", heteRate)
print ("extendLen", extendLen)
print ("heterozygous kmer coverage range", lowCov, highCov)
logging.info("heterozygous kmer coverage range %s %s" % (lowCov, highCov))
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

edges.extend( build_graph.non_snp_edges(kmerCov, k, left_index, right_index, extendKmers) )
time4 = time.clock()
logging.info( "add non snp edges, cost %.2f second" % (time4 - time3) )
logging.info( "edge number %s" % ( len(edges) )  )
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


graphs = list(nx.connected_component_subgraphs(G))
#---------------------------------------------------------------

# for writing paper
#kmerInGraph =set()
#kmerInSelected = set()
selectedMates = {}
weightDis = {}
for g in graphs:
    mate = nx.max_weight_matching(g) 
    for (ele1, ele2) in iter(mate):
        w = g.get_edge_data(ele1, ele2)['weight']
        selectedMates[(ele1, ele2)] = w
        if w not in weightDis:
            weightDis[w] = 0            
        weightDis[w] += 1


time8 = time.clock()
logging.info( "compute max weight matching, cost %.2f seconds" % (time8 - time7) )

sortedWeightDis = sorted(weightDis.items())
#print(sortedWeightDis)

weightThreshold = tools.calc_weightThreshold(heteRate,k)

#weightThreshold = 4 # for simulated dataset
    

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

snpFile = "k_" + str(k) + "_pair.snp"
nonFile = "k_" + str(k) + "_pair.non"
nonSepFile = "k_" + str(k) + "_pair.non.sep"
snpOut = open(snpFile, "w")
nonOut = open(nonFile, "w")
nonSepOut = open(nonSepFile, "w")


sSnpPair, sNonPair = sorted(snpPair), sorted(nonPair)
for (ele1, ele2, w) in iter(sSnpPair):
    snpOut.write("%s %s %s\n" % (ele1, ele2, w) )
snpOut.close()


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


sNon = sorted(list(nonSepPairKmer))
for (ele1, ele2, w) in iter(sNon):
    nonSepOut.write("%s %s %s\n" % (ele1, ele2, w) )
nonSepOut.close()
