#########################################################################
# File Name: transfer_bowtie2sam_into_vcf.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Mon 07 Jun 2021 10:12:21 PM AEST
#########################################################################
#!/bin/bash
import os
import sys

def write_VCF_header(fout):
    fout.write("##fileformat=VCFv4.1\n")
    ##filedate=20191017
    fout.write("##REF=<ID=REF,Number=1,Type=String,Description=\"Allele of reference genome\">\n")
    fout.write("##REF=<ID=ALT,Number=1,Type=String,Description=\"Allele of non-reference genome\">\n")
    fout.write("##INFO=<ID=Ty,Number=1,Type=String,Description=\"SNP, INS, DEL or .\">\n")
    fout.write("##INFO=<ID=Rk,Number=1,Type=Float,Description=\"SNP rank\">\n")
    ##INFO=<ID=UL,Number=1,Type=Integer,Description="length of the unitig left">
    ##INFO=<ID=UR,Number=1,Type=Integer,Description="length of the unitig right">
    ##INFO=<ID=CL,Number=1,Type=Integer,Description="length of the contig left">
    ##INFO=<ID=CR,Number=1,Type=Integer,Description="length of the contig right">
    ##INFO=<ID=Genome,Number=1,Type=String,Description="Allele of the reference;for indel reference is . ">
    fout.write("##INFO=<ID=Sd,Number=1,Type=Integer,Description=\"Reverse (-1) or Forward (1) Alignement\">\n")
    fout.write("##INFO=<ID=XA,Number=.,Type=String,Description=\"Other mapping positions (chromosome_position). Position is negative in case of Reverse alignment. The position designs the starting position of the alignment, not the position of the variant itself.\">\n")
    fout.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Cumulated depth accross samples (sum)">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled Genotype Likelihoods">
    ##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Depth of each allele by sample">
    fout.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	G1\n")
    return

def get_snp_pos(s1, s2):
    pos = []
    lenS = len(s1)
    assert lenS == len(s2)
    for i in range(lenS):
        if s1[i] == -1 or s2[i] == -1:
            continue
        if s1[i] != s2[i]:
            pos.append(i)
    return pos
if __name__ == '__main__':   

    if len(sys.argv) < 3:
        print ("kmerPair.sam kmerPair.vcf")
        sys.exit()
 
    fout = open(sys.argv[2], "w")
    write_VCF_header(fout)

    cnt = 0
    with open(sys.argv[1], "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            words = line.strip().split()
            cnt += 1
            if cnt % 2 == 1:
                chrom1, pos1 = words[2], int(words[3])
                seq1 = words[9]
                k = str(len(seq1))
                ID1 = words[0].split('_')[0]
                for ele in words[11:]:
                    if ele.startswith("MD"):
                        MD1 = ele.split(':')[-1]
            else:
                chrom2, pos2 = words[2], int( words[3])
                seq2 = words[9]
                for ele in words[11:]:
                    if ele.startswith("MD"):
                        MD2 = ele.split(':')[-1]
                ID2 = words[0].split('_')[0]
                assert ID1 == ID2
                pos = get_snp_pos(seq1, seq2)
                assert len(pos) == 1
                p = pos[0]
                if (MD1 != k and MD2 != k):
                    print (ID1, chrom1, pos1, MD1, chrom2, pos2, MD2)
                    #sys.exit()
                    if MD1 == MD2: 
                        fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1, 
                            pos1, ID1, MD1[int(len(MD1)/2)] , seq1[p]+"," +seq2[p], 
                            '.', '.', "Ty=SNP,Rk=0", "GT", "1/2") )
                if MD1 == k:
                    ref = seq1
                    alt = seq2
                    fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1, 
                        pos1, ID1, ref[p], alt[p], '.', '.', "Ty=SNP,Rk=0", "GT", "0/1") )
                elif MD2 == k:
                    ref = seq2
                    alt = seq1
                    fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1, 
                        pos1, ID2, ref[p], alt[p], '.', '.', "Ty=SNP,Rk=0", "GT", "0/1") )
    fout.close()                

        
