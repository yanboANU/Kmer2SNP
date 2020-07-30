#########################################################################
# File Name: runkmercalling.sh
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Thu 04 Jul 2019 16:52:18 AEST
#########################################################################
#!/bin/bash
set -e

if [ $# != 3 ]; then
	echo "\$1:kmer-size \$2:homo coverage \$3:input.fq[,input2.fq][,input3.fq]"
	exit 
fi	

echo "\$1:kmer-size \$2:homo coverage \$3:input.fq"
echo $1 $2 $3

if [ ! -f "chr_k$1.txt" ]; then 
	/path2dsk/dsk/build/bin/dsk -nb-cores 10 -file $3 -histo 1 -out chr_k$1 -kmer-size $1
	/path2dsk/dsk/build/bin/dsk2ascii -nb-cores 10 -file chr_k$1 -out chr_k$1.txt
fi

if [ ! -f "hete.peak.k$1" ]; then 
	Rscript /path2Kmer2SNP/Kmer2SNP/libprism/local/runfindGSE.r chr_k$1.histo $1 ./ $2 >hete.peak.k$1
fi

if [ -f "hete.peak.k$1" ]; then
	left=`cat hete.peak.k$1 | grep "het_xfit_left" | awk '{print $6}'`
	right=`cat hete.peak.k$1 | grep "het_xfit_right" | awk '{print $6}'`
fi

if [ -f "chr_k$1.h5" ]; then 
    rm chr_k$1.h5
fi

echo $left
echo $right

command="python3 /path2Kmer2SNP/Kmer2SNP/kmerGraphCalling.py --t1 chr_k$1.txt --c1 $left --c2 $right --k $1 >vc_k$1.log"
echo $command
python3 /path2Kmer2SNP/Kmer2SNP/kmerGraphCalling.py --t1 chr_k$1.txt --c1 $left --c2 $right --k $1 >vc_k$1.log


