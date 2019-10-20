#########################################################################
# File Name: eva.sh
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Mon 14 Oct 2019 11:42:27 AEDT
#########################################################################
#!/bin/bash
python ~/software/Kmer2SNP/libprism/evaluate/compareRealPair.py /media/yanbo/Data/VariationCalling/experiment/realPair/chr22/chr22.snp.real.31mer k_31_pair.snp

python ~/software/Kmer2SNP/libprism/evaluate/compareRealPair.py /media/yanbo/Data/VariationCalling/experiment/realPair/chr22/chr22.nonIsolated.snp.real.31mer k_31_pair.non

#python ~/software/Kmer2SNP/libprism/evaluate/compareRealPair.py /media/yanbo/Data/VariationCalling/experiment/realPair/chr22/chr22.indel.real.31mer k_31_pair.indel

