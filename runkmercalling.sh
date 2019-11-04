#########################################################################
# File Name: runkmercalling.sh
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Thu 04 Jul 2019 16:52:18 AEST
#########################################################################
#!/bin/bash
set -e
#cd $$3x/chr$1/
#mkdir kmercalling
#cd kmercalling

if [ $# != 6 ]; then
	echo "\$1:chrID \$2:k \$3:k-1 \$4:hom cov \$5:*fq \$6: NGS/TGS "
	exit 
fi	

echo "\$1:chrID \$2:k \$3:k-1 \$4:hom cov \$5:*fq \$6: NGS/TGS"
echo $1 $2 $3 $4 $5 $6
start=`date +%s`

if [ ! -f "chr$1_k$2.txt" ]; then 
	/home/yulin/software/dsk/build/bin/dsk -nb-cores 10 -file $5 -histo 1 -out chr$1_k$2 -kmer-size $2
	/home/yulin/software/dsk/build/bin/dsk2ascii -nb-cores 10 -file chr$1_k$2 -out chr$1_k$2.txt

	#/home/yulin/software/dsk/build/bin/dsk -nb-cores 2 -file $5 -histo 1 -out chr$1_k$3 -kmer-size $3
	#/home/yulin/software/dsk/build/bin/dsk2ascii -nb-cores 2 -file chr$1_k$3 -out chr$1_k$3.txt
fi
#exec echo hello
if [ "$6" = "1" ]; then
	left=8
    right=24
elif [ ! -f "hete.peak.k$2" ]; then 
	#findGSE:
	Rscript /home/yulin/bio/VariationCalling/libprism/local/runfindGSE.r chr$1_k$2.histo $2 ./ $4 >hete.peak.k$2
fi

if [ -f "hete.peak.k$2" ]; then
	left=`cat hete.peak.k$2 | grep "het_xfit_left" | awk '{print $6}'`
	right=`cat hete.peak.k$2 | grep "het_xfit_right" | awk '{print $6}'`
fi

if [ -f "chr$1_k$2.h5" ]; then 
    rm chr$1_k$2.h5
fi

echo $left
echo $right

#python /home/yanbo/software/Kmer2SNP/variationCalling.py --t1 chr$1_k$2.txt --t2 chr$1_k$3.txt --c1 $left --c2 $right --k $2 --b $6 >vc_k$2.log

command="/home/yulin/py36/bin/python3 /home/yanbo/software/Kmer2SNP/kmerGraphCalling.py --t1 chr$1_k$2.txt --c1 $left --c2 $right --k $2 --b $6 >vc_k$2.log"
echo $command
#/home/yulin/py36/bin/python3 /home/yanbo/software/Kmer2SNP/kmerGraphCalling.py --t1 chr$1_k$2.txt --c1 $left --c2 $right --k $2 --b $6 >vc_k$2.log
/home/yulin/py36/bin/python3 /home/yanbo/software/Kmer2SNP/kmerGraphCalling.py --t1 chr$1_k$2.uniq.kmer --c1 $left --c2 $right --k $2 --b $6 >vc_k$2.log


end=`date +%s`
runtime=$((end-start))
echo run $runtime seconds
##hour=$(echo "$runtime/3600" | bc)
##echo run $hour hours
#ret=$?; times; exit "$ret"
##trap times EXIT
#bash -c 'trap times EXIT; : {1..1000000}'
#zsh -c 'trap times EXIT; : {1..1000000}'

