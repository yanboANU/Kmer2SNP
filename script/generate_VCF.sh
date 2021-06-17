#########################################################################
# File Name: generate_VCF.sh
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Mon 07 Jun 2021 09:59:54 PM AEST
#########################################################################
#!/bin/bash


if [ $# != 3 ]; then
	echo "\$1:kmer-pair-file \$2: ref.fasta \$3: threads" # \4: k-mer k"
	exit 
fi	
echo "\$1:kmer-pair-file \$2: ref.fasta \$3: threads" # \4: k-mer k"
echo $1 $2 
#bowtie2
kmerPair=$1
ref=$2
index=`ls $ref | cut -f 1 -d '.'`
threads=$3
#k=$4
echo $ref $index
exit
bowtie2-build $2 $index
python3 /path2Kmer2SNP/libprism/supplementary_tools/write_kmerPair_2fastq.py $kmerPair kmerPair.fastq
bowtie2 --threads $threads  --ignore-quals -i L,1,0 -N 0 -L 13 --mp 1,1 --ma 0 --np 0 --score-min L,0,-0.15 --end-to-end -R 19 -I 0 -X 1 --ff --reorder --mm --seed 12321231 -x $index --interleaved kmerPair.fastq -S output.sam

samtools sort -O SAM output.sam >sorted_output.sam
python3 /path2Kmer2SNP/libprism/supplementary_tools/transfer_bowtie2sam_into_vcf.py output.sam kmerpair.vcf 




#python3 /path2Kmer2SNP/libprism/supplementary_tools/write_kmerPair_2fastq.py k_31_pair.snp.extend extend.kmerPair.fastq
#bowtie2 --threads 12  --ignore-quals -i L,1,0 -N 0 -L 13 --mp 1,1 --ma 0 --np 0 --score-min L,0,-0.15 --end-to-end -R 19 -I 0 -X 1 --ff --reorder --mm --seed 12321231 -x $index --interleaved extend.kmerPair.fastq -S extend.output.sam
