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
args = parser.parse_args()
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

#############################################################################
lowCov, highCov = int(args.c1), int(args.c2)
k = int(args.k)
extendLen = 1 
lowCov = lowCov - extendLen
highCov = highCov + extendLen

logging.info("heterozygous kmer coverage range %s %s" % (lowCov, highCov))

kmerCov = kmercalling.pick_smaller_unique_kmer( args.t1, 
           lowCov, highCov)


logging.info( "picked heterozygous Kmer number: %s" % len(kmerCov) )
left_index, right_index = kmercalling.build_left_right_kmer_index(kmerCov)
edges = []
extendKmers = {}
edges.extend( build_graph.snp_edges(kmerCov, k, left_index, right_index, extendKmers) )

edges.extend( build_graph.non_snp_edges(kmerCov, k, left_index, right_index, extendKmers) )

G=nx.Graph()
G.add_weighted_edges_from(edges)
logging.info("number of components: %s" % nx.number_connected_components(G) )


graphs = list(nx.connected_component_subgraphs(G))


# for writing paper
#kmerInGraph =set()
#kmerInSelected = set()
#componentSize = {}


indelPair, snpPair, nonPair = [], [], []
marknon = {}
count = 0
selectedMates = {}
weightDis = {}
for g in graphs:
    mate = nx.max_weight_matching(g)
    '''
    gSize = g.size()
    if (gSize >= 1):
        if gSize not in componentSize:
            componentSize[gSize] = 1
        else:
            componentSize[gSize] += 1
        for node in g.nodes():
            kmerInGraph.add( min(node, tools.reverse(node)) )
    '''        
    for (ele1, ele2) in mate:
        w = g.get_edge_data(ele1, ele2)['weight']
        selectedMates[(ele1, ele2)] = w
        if w not in weightDis:
            weightDis[w] = 0            
        weightDis[w] += 1


sortedWeightDis = sorted(weightDis.items())
l=len(sortedWeightDis)
#print(sortedWeightDis)
temp = sortedWeightDis[0][1]
deleta = []
#weightThreshold = 0
for i in range(1,l):
    deleta.append(sortedWeightDis[i-1][1] - sortedWeightDis[i][1])

#print ("deleta", deleta)

weightThreshold = 2
for i in range(1,l-2):
    if deleta[i] < deleta[i-1] and deleta[i] < deleta[i+1]:
        weightThreshold = sortedWeightDis[i+1][0] 
        break
    
    if deleta[i] < 0:
        weightThreshold = sortedWeightDis[i][0] 
        break
    

logging.info("weight threshold %s" % weightThreshold)   
for (ele1, ele2) in selectedMates:
    w= selectedMates[(ele1,ele2)] 
    if w < weightThreshold:  
        continue
    if len(ele1) > len(ele2):
        indelPair.append( ( ele1, ele2, w ) )
    elif len(ele1) < len(ele2):
        indelPair.append( ( ele2, ele1, w ) )
    else:    
        if tools.hamming_distance(ele1, ele2)==1:
            if ele1 <= ele2:
                snpPair.append( ( ele1, ele2, w ) )
            else:
                snpPair.append( (ele2, ele1, w ) )

        if tools.hamming_distance(ele1, ele2)==2:
            if ele1 <= ele2:
                temp = extendKmers[(ele1, ele2)]
            else:
                temp = extendKmers[(ele2, ele1)]
            if temp not in marknon:
                nonPair.append( (temp[0], temp[1], w ) )
                marknon[temp] = 1
            else:
                marknon[temp] += 1

# for write paper
#fout = open("kmerInGraph.txt", "w")
#sortedKmerInGraph = sorted( list(kmerInGraph) )
#for kmer in sortedKmerInGraph:
    #fout.write("%s\n" % kmer)
#fout.close()
#
#
#fout = open("kmerInSelected.txt", "w")
#sortedKmerInSelected = sorted( list(kmerInSelected) )
#for kmer in sortedKmerInSelected:
#    fout.write("%s\n" % kmer)
#fout.close()



snpFile = "k_" + str(k) + "_pair.snp"
nonFile = "k_" + str(k) + "_pair.non"
nonSepFile = "k_" + str(k) + "_pair.non.sep"
#indelFile = "k_" + str(k) + "_pair.indel"
snpOut = open(snpFile, "w")
nonOut = open(nonFile, "w")
nonSepOut = open(nonSepFile, "w")
#indelOut = open(indelFile, "w")


sSnpPair, sNonPair = sorted(snpPair), sorted(nonPair)
for (ele1, ele2, w) in sSnpPair:
    snpOut.write("%s %s %s\n" % (ele1, ele2, w) )
snpOut.close()


nonSepPairKmer = set()
for (ele1, ele2, w) in sNonPair:
    if marknon[(ele1, ele2)] >= 2: # both selected, then think it's true
        nonOut.write("%s %s %s\n" % (ele1, ele2, w)) 
        h1Pre, h2Pre, h1Suf, h2Suf = ele1[:k], ele2[:k], ele1[-k:], ele2[-k:]
        smallerH1P, smallerH2P = tools.get_smaller_pair_kmer(h1Pre, h2Pre)
        smallerH1S, smallerH2S = tools.get_smaller_pair_kmer(h1Suf, h2Suf)
        nonSepPairKmer.add( (smallerH1P, smallerH2P) )
        nonSepPairKmer.add( (smallerH1S, smallerH2S) )
nonOut.close()

sNon = sorted(list(nonSepPairKmer))
for (ele1, ele2) in sNon:
    nonSepOut.write("%s %s\n" % (ele1, ele2) )
nonSepOut.close()

