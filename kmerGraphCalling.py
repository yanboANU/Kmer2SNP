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
#from libprism.local.prepare import clouds_from_refhap, merge_clouds, print_clouds
from libprism.evaluate import read
#from math import log, exp
#from pysam
import logging
import time
import networkx as nx
#from networkx.generators.atlas import *
#from networkx.algorithms.isomorphism.isomorph import graph_could_be_isomorphic as isomorphic
#import random
import draw
#############################################################################

parser = argparse.ArgumentParser(description="Variation calling for Single-individual based on NGS reads")

#parser.add_argument('--bam', help='path to alignment bam file', required=True)
#parser.add_argument('--ref', help='reference or scaffolds', required=True)
parser.add_argument('--t1', help='kmer freq txt file (dsk result)', required=True)
parser.add_argument('--t2', help='k-1mer freq txt file (dsk result)', required=False)
parser.add_argument('--c1', help='low coverage', required=True)
parser.add_argument('--c2', help='high coverage', required=True)
parser.add_argument('--k', help='kmer size', required=True)
#parser.add_argument('--b', help='for simulate data 0, real data 1, can dicide indel threshold', required=True)
parser.add_argument('--b', help='for NGS 0, TGS 1, can dicide indel threshold', required=True)
#parser.add_argument('--o', help='output file prefix', required=True)

args = parser.parse_args()
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

#############################################################################
time1 = time.clock()
lowCov, highCov = int(args.c1), int(args.c2)
k = int(args.k)
extendLen = 1 
#if k>45:
#    extendLen = 4
#elif k>41:
#    extendLen = 3
#elif k>31:
#    extendLen = 2
    
lowCov = lowCov - extendLen
highCov = highCov + extendLen


print ("heterozygous kmer coverage range", lowCov, highCov)

#kmerCov = kmercalling.pick_smaller_unique_kmer( args.t1, 
#           lowCov, highCov)
#k_1merCov = kmercalling.pick_smaller_unique_kmer( args.t2, 
#             lowCov, highCov)

time2 = time.clock()
#uniq.kmer
kmerCov = read.read_2_columns(args.t1) 
#k_1merCov = read.read_2_columns(args.t2)

logging.info( "uniq Kmer size: %s" % len(kmerCov) )
#logging.info( "uniq (K-1)mer size: %s" % len(k_1merCov) )
logging.info( "finish reading (kmer cov) and (k-1mer cov) file, cost %.2f seconds" % (time2-time1))
left_index, right_index = kmercalling.build_left_right_kmer_index(kmerCov)
edges = []
extendKmers = {}
edges.extend( build_graph.snp_edges(kmerCov, k, left_index, right_index, extendKmers, args.b) )
time3 = time.clock()
logging.info( "add snp edges, cost %.2f second" % (time3 - time2) )
print ("extendKmers size", len(extendKmers))
#if args.b == '0':
edges.extend( build_graph.non_snp_edges(kmerCov, k, left_index, right_index, extendKmers, args.b) )
time4 = time.clock()
logging.info( "add non snp edges, cost %.2f second" % (time4 - time3) )
print ("extendKmers size", len(extendKmers))

'''
edges.extend( build_graph.indel_edges(kmerCov, k_1merCov, k, left_index, right_index, extendKmer) )
time5 = time.clock()
logging.info( "add indel edges, cost %.2f seconds" % (time5 - time4) )
'''

G=nx.Graph()
G.add_weighted_edges_from(edges)
#draw.first_try(G)
#sys.exit()
time6 = time.clock()
logging.info( "build graph, cost %.2f seconds" % (time6 - time3) )
print ("number of components: ", nx.number_connected_components(G) )
time7 = time.clock()
logging.info( "compute components, cost %.2f seconds" % (time7 - time6) )

graphs = list(nx.connected_component_subgraphs(G))

indelPair, snpPair, nonPair = [], [], []
marknon = {}
count = 0

# for writing paper
kmerInGraph =set()
#kmerInSelected = set()
componentSize = {}
for g in graphs:
    mate = nx.max_weight_matching(g)
    # for writing paper
    gSize = g.size()
    if (gSize >= 1):
        #print ("connected component size: ", g.size())
        #print ("heterozygous kmer pair number: ", len(mate))
        if gSize not in componentSize:
            componentSize[gSize] = 1
        else:
            componentSize[gSize] += 1

        for node in g.nodes():
            kmerInGraph.add( min(node, tools.reverse(node)) )
    '''  
    if (g.size() > 1):
        count += 1
        draw.first_try(g, "fig"+str(count))
    '''
    for (ele1, ele2) in mate:
        w = g.get_edge_data(ele1, ele2)['weight']
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

                #    nonPair.append( (ele1, ele2, w ) )
                #else:
                #    nonPair.append( (ele2, ele1, w ) )

#print ( "marnon", marknon )
print ( "component size", sorted(componentSize.items()) )
# for write paper
fout = open("kmerInGraph.txt", "w")
sortedKmerInGraph = sorted( list(kmerInGraph) )
for kmer in sortedKmerInGraph:
    fout.write("%s\n" % kmer)
fout.close()
#
#
#fout = open("kmerInSelected.txt", "w")
#sortedKmerInSelected = sorted( list(kmerInSelected) )
#for kmer in sortedKmerInSelected:
#    fout.write("%s\n" % kmer)
#fout.close()


time8 = time.clock()
logging.info( "compute max weight matching, cost %.2f seconds" % (time8 - time7) )

#fout.write(">kmer_snp%s_1_cov_%s\n" % (ID, c1))
#fout.write("%s\n" % ( k1 ) )
#fout.write(">kmer_snp%s_2_cov_%s\n" % (ID, c2))
#fout.write("%s\n" % ( k2 ) )

snpFile = "k_" + str(k) + "_pair.snp"
nonFile = "k_" + str(k) + "_pair.non"
nonSepFile = "k_" + str(k) + "_pair.non.sep"
indelFile = "k_" + str(k) + "_pair.indel"
snpOut = open(snpFile, "w")
nonOut = open(nonFile, "w")
nonSepOut = open(nonSepFile, "w")
indelOut = open(indelFile, "w")


sIndelPair, sSnpPair, sNonPair = sorted(indelPair), sorted(snpPair), sorted(nonPair)
for (ele1, ele2, w) in sSnpPair:
    #kmerInSelected.add( min(ele1, tools.reverse(ele1)) )
    #kmerInSelected.add( min(ele2, tools.reverse(ele2)) )
    snpOut.write("%s %s %s\n" % (ele1, ele2, w) )
snpOut.close()

for (ele1, ele2, w) in sIndelPair:
    indelOut.write("%s %s %s\n" % (ele1, ele2, w))
indelOut.close()

nonSepPairKmer = set()
for (ele1, ele2, w) in sNonPair:
    if marknon[(ele1, ele2)] >= 2: # both selected, then think it's true
        nonOut.write("%s %s %s\n" % (ele1, ele2, w)) 
        h1Pre, h2Pre, h1Suf, h2Suf = ele1[:k], ele2[:k], ele1[-k:], ele2[-k:]
        smallerH1P, smallerH2P = tools.get_smaller_pair_kmer(h1Pre, h2Pre)
        smallerH1S, smallerH2S = tools.get_smaller_pair_kmer(h1Suf, h2Suf)
        nonSepPairKmer.add( (smallerH1P, smallerH2P) )
        nonSepPairKmer.add( (smallerH1S, smallerH2S) )
        #kmerInSelected.add( min(h1Pre, tools.reverse(h1Pre)) )
        #kmerInSelected.add( min(h2Pre, tools.reverse(h2Pre)) )
        #kmerInSelected.add( min(h1Suf, tools.reverse(h1Suf)) )
        #kmerInSelected.add( min(h2Suf, tools.reverse(h2Suf)) )
nonOut.close()

sNon = sorted(list(nonSepPairKmer))
for (ele1, ele2) in sNon:
    nonSepOut.write("%s %s\n" % (ele1, ele2) )
nonSepOut.close()

