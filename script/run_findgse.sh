#########################################################################
# File Name: script/run_dsk.sh
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Mon 07 Jun 2021 12:46:48 PM AEST
#########################################################################
#!/bin/bash


#findGSE:
if [ ! -f "hete.peak.k$1" ]; then
    out=$(Rscript /path2Kmer2SNP/libprism/local/runfindGSE.r chr_k$1.histo $1 ./ $2) && echo "$out" > hete.peak.k$1
fi

heteRate=`cat v1*txt | grep "Het_rate" | awk '{print $2}' | head -1`
echo "heterozygous rate is" $heteRate >hete.para


if [ -f "hete.peak.k$1" ]; then
	left=`cat hete.peak.k$1 | grep "het_xfit_left" | awk '{print $6}'`
	right=`cat hete.peak.k$1 | grep "het_xfit_right" | awk '{print $6}'`
fi

echo "heterozygous lowCoverage is" $left >>hete.para
echo "heterozygous highCoverage is" $right  >> hete.para
