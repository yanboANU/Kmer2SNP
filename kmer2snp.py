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
import time
import networkx as nx

from libprism.local import build_graph
from libprism.local import kmercalling
from libprism.local import tools
from libprism.local import process
from libprism.local import generate_VCF_without_reference

def read_coverage(filename):
    covMap = {}
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split(' ')
            sampleID = words[0]
            cov=words[1]
            covMap[sampleID] = cov         
    return covMap

def run_population(args):
    
    covMap = read_coverage(args.cfile)
    prePath = os.path.abspath(args.output_dir)    
    with open(args.faqfile, "r") as f:
        for line in f:
            words = line.strip().split(' ')
            sampleID = words[0]
            faqFile = words[1]
            if args.b == None:
                args.b = 0
            print ("debug", args.output_dir)
            path = prePath +"/"+sampleID
            os.system('mkdir ' + path)
            #os.system("cd "+ path)
            os.chdir(path)
            os.system('pwd')
            #sys.exit()
            run_single(args.k, covMap[sampleID], faqFile, args.b)

    os.chdir(prePath)
    generate_VCF_without_reference.run(args.cfile)
    return
    
    



def run_single(argsk, argsc, argsfastaq, argsb=1, argst1=None, argsc1=None, argsc2=None, argsr=None):

    if argst1 == None:
        command="sh " + os.path.dirname(os.path.realpath(__file__)) + "/script/run_dsk.sh " + argsk + " " + argsc + " " + argsfastaq
        print (command)
        os.system( command )
        argst1 = "chr_k" + argsk  + ".txt"

    if argsc1 == None or argsc2 == None or argsr == None:
        command="sh " + os.path.dirname(os.path.realpath(__file__)) + "/script/run_findgse.sh " + argsk + " " + argsc + " " + argsfastaq
        print (command)
        os.system( command )
        c1, c2, r = tools.read_findGSE_result("hete.para") 
        if argsc1 == None:
            argsc1 = c1
        if argsc2 == None:
            argsc2 = c2
        if argsr == None:
            argsr = r
    if argsb == None:
        argsb = 1
    process.run(argst1, argsc1, argsc2, argsr, argsk, int(argsb))




parser = argparse.ArgumentParser(description="Variation calling for Single-individual based on NGS reads")

# Modes:
sample_parsers = parser.add_subparsers(title="sample_mode")

#

parser_single = sample_parsers.add_parser('single')
parser_single.add_argument('--k', help='kmer size', required=True)
parser_single.add_argument('--c', help='homozygous coverage', required=True)
parser_single.add_argument('--fastaq', help='NGS reads in fasta or fastq format, gzipped acceptable. Separate multiple file paths by a comma (",") ', required=True)
parser_single.add_argument('--output_dir', help='output directory (Default is ./)', required=False)
parser_single.add_argument('--b', help='b=0 call only isolated SNP (one SNP in a k-mer), b=1 call at most 2 SNPs in a k-mer (default b = 1)', required=False)
parser_single.add_argument('--t1', help='kmer frequency txt file (Default is DSK result)', required=False)
parser_single.add_argument('--c1', help='low coverage threshold (default: FindGSE result)', required=False)
parser_single.add_argument('--c2', help='high coverage threshold (default: FindGSE result)', required=False)
parser_single.add_argument('--r', help='heterozygous rate (default: FindGSE result)', required=False)
parser_single.set_defaults(which='single')


# population sample mode

parser_population = sample_parsers.add_parser('population')
parser_population.add_argument('--k', help='kmer size', required=True)
parser_population.add_argument('--cfile', help='a file contains homozygous coverage of NGS reads for each sample', required=True)
parser_population.add_argument('--faqfile', help='a file contains the Next Generation Sequencing (NGS) reads files for each sample', required=True)
parser_population.add_argument('--output_dir', help='output directory (Default is ./)', required=False)
parser_population.add_argument('--b', help='b=0 only call isolated SNV (only one snv in a k-mer), b=1 call at most 2 SNVs in a k-mer (Default b = 0)', required=False)
#parser_population.add_argument('--t1file', help=' a file contains kmer freq txt file (Default is DSK result)', required=False)
#parser_population.add_argument('--c1file', help='a file contains the estimated minimum low coverage (Default is FindGSE result)', required=False)
#parser_population.add_argument('--c2file', help='a file contains the estimated maximum  coverage (Default is FindGSE result)', required=False)
#parser_population.add_argument('--rfile', help='a file contains the estimated heterozygous rate (Default is FindGSE result)', required=False)
parser_population.add_argument('--ref', help='a reference genome file', required=False)

parser_population.set_defaults(which='population')




args = parser.parse_args()

if args.which == 'population':
    print (args)
    if args.output_dir == None:
        args.output_dir = './'
    run_population(args)
    #return
else:
    print (args)
    if args.output_dir == None:
        args.output_dir = './'
    #os.system('cd ' + args.output_dir)
    os.chdir( args.output_dir )
    os.system('pwd')
    run_single(args.k, args.c, args.fastaq, 
            args.b, args.t1, args.c1, args.c2, args.r)
    
