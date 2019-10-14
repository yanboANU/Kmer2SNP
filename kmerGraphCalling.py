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
#############################################################################

parser = argparse.ArgumentParser(description="Variation calling for Single-individual based on NGS reads")

#parser.add_argument('--bam', help='path to alignment bam file', required=True)
#parser.add_argument('--ref', help='reference or scaffolds', required=True)
parser.add_argument('--t1', help='kmer freq txt file (dsk result)', required=True)
parser.add_argument('--t2', help='k-1mer freq txt file (dsk result)', required=True)
parser.add_argument('--c1', help='low coverage', required=True)
parser.add_argument('--c2', help='high coverage', required=True)
parser.add_argument('--k', help='kmer size', required=True)
parser.add_argument('--b', help='for simulate data 0, real data 1, can dicide indel threshold', required=True)
#parser.add_argument('--o', help='output file prefix', required=True)

args = parser.parse_args()
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

#############################################################################
time1 = time.clock()
lowCov, highCov = int(args.c1), int(args.c2)
lowCov = lowCov - 1
highCov = highCov + 1

k = int(args.k)
'''
kmerCov = kmercalling.pick_smaller_unique_kmer( args.t1, 
             lowCov, highCov)
k_1merCov = kmercalling.pick_smaller_unique_kmer( args.t2, 
             lowCov, highCov)
'''
time2 = time.clock()
#uniq.kmer
kmerCov = read.read_2_columns(args.t1) 
k_1merCov = read.read_2_columns(args.t2)

logging.info( "uniq Kmer size: %s" % len(kmerCov) )
logging.info( "uniq (K-1)mer size: %s" % len(k_1merCov) )
logging.info( "finish reading (kmer cov) and (k-1mer cov) file, cost %.2f seconds" % (time2-time1))
left_index, right_index = kmercalling.build_left_right_kmer_index(kmerCov)
edges = []
extendKmers = {}
edges.extend( build_graph.snp_edges(kmerCov, k, left_index, right_index, extendKmers) )
time3 = time.clock()
logging.info( "add snp edges, cost %.2f second" % (time3 - time2) )

edges.extend( build_graph.non_snp_edges(kmerCov, k, left_index, right_index, extendKmers) )
time4 = time.clock()
logging.info( "add non snp edges, cost %.2f second" % (time4 - time3) )

'''
edges.extend( build_graph.indel_edges(kmerCov, k_1merCov, k, left_index, right_index, extendKmer) )
time5 = time.clock()
logging.info( "add indel edges, cost %.2f seconds" % (time5 - time4) )
'''

G=nx.Graph()
G.add_weighted_edges_from(edges)
time6 = time.clock()
logging.info( "build graph, cost %.2f seconds" % (time6 - time4) )
print ("number of components: ", nx.number_connected_components(G) )
time7 = time.clock()
logging.info( "compute components, cost %.2f seconds" % (time7 - time6) )

graphs = list(nx.connected_component_subgraphs(G))

indelPair, snpPair, nonPair = [], [], []
for g in graphs:
    print ("connected component size: ", g.size())
    mate = nx.max_weight_matching(g)
    print ("heterozygous kmer pair number: ", len(mate))

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
                temp = extendKmers[(ele1, ele2)]
                nonPair.append( (temp[0], temp[1], w ) )
                #if ele1 <= ele2:
                #    nonPair.append( (ele1, ele2, w ) )
                #else:
                #    nonPair.append( (ele2, ele1, w ) )

time8 = time.clock()
logging.info( "compute max weight matching, cost %.2f seconds" % (time8 - time7) )

#fout.write(">kmer_snp%s_1_cov_%s\n" % (ID, c1))
#fout.write("%s\n" % ( k1 ) )
#fout.write(">kmer_snp%s_2_cov_%s\n" % (ID, c2))
#fout.write("%s\n" % ( k2 ) )

snpFile = "k_" + str(k) + "_pair.snp"
nonFile = "k_" + str(k) + "_pair.non"
indelFile = "k_" + str(k) + "_pair.indel"
snpOut = open(snpFile, "w")
nonOut = open(nonFile, "w")
indelOut = open(indelFile, "w")

sIndelPair, sSnpPair, sNonPair = sorted(indelPair), sorted(snpPair), sorted(nonPair)
for (ele1, ele2, w) in sSnpPair:
    snpOut.write("%s %s %s\n" % (ele1, ele2, w) )
snpOut.close()

for (ele1, ele2, w) in sIndelPair:
    indelOut.write("%s %s %s\n" % (ele1, ele2, w))
indelOut.close()

for (ele1, ele2, w) in sNonPair:
    nonOut.write("%s %s %s\n" % (ele1, ele2, w))
nonOut.close()

