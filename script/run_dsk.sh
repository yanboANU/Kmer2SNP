#########################################################################
# File Name: script/run_dsk.sh
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Mon 07 Jun 2021 12:46:48 PM AEST
#########################################################################
#!/bin/bash


if [ ! -f "chr_k$1.txt" ]; then 
	/path2dsk/dsk -nb-cores 10 -file $3 -histo 1 -out chr_k$1 -kmer-size $1
	/path2dsk/dsk2ascii -nb-cores 10 -file chr_k$1 -out chr_k$1.txt
fi

if [ -f "chr_k$1.h5" ]; then 
    rm chr_k$1.h5
fi
