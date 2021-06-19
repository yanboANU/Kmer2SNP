#########################################################################
# File Name: evaluate.sh
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Thu 17 Jun 2021 04:14:01 PM AEST
#########################################################################
#!/bin/bash


if [ $# != 5 ]; then
	echo "\$1:kmer-pair-file \$2: ref.fasta (or h1.fasta) \$3: group-truth vcf (or h2.fasta) \4: k-mer k size \5: chrID"
	exit 
fi	


findPair=$1
h1Fa=$2
refFa=$2
h2Fa=$3  
vcf=$3
path=$(readlink -f "$(dirname "$3")")
chrID=$5
k=$4


if [ ! -d "$path/${chrID}_kmerPair" ]; then
    mkdir $path/${chrID}_kmerPair
fi

path=$path/${chrID}_kmerPair

if echo "$3" | grep -q -E '\.vcf'; then
	
	if [ ! -f "$path/k_${k}_pair.snp" ]; then
		python3 /path2Kmer2SNP/libprism/evaluate/getKmerPair_FromVCF_REF.py $refFa $k $vcf $path
	fi
else
	if [ ! -f "$path/ref1_ref2.blasr" ]; then
		/usr/bin/time -f "%U %M %e" blasr $h1Fa $h2Fa -m 5 --nproc 10 --minMatch 15 --maxMatch 20 --advanceExactMatches 10 --fastSDP > $path/ref1_ref2.blasr
	fi

	if [ ! -f "$path/k_${k}_pair.snp" ]; then 
		python3 /path2Kmer2SNP/libprism/supplementary_tools/pick_KmerPair_sep_kmer_from_blasrM5.py $k $path/ref1_ref2.blasr $path >pick_${k}mer.log
	fi
fi	
	
python3 /path2Kmer2SNP/libprism/evaluate/compareRealPair.py $path/k_${k}_pair.snp $1
