#!/home/yulin/software/R-3.6.0/bin/Rscript

args <- commandArgs(TRUE)
library("findGSE")
findGSE(histo=args[1], sizek=args[2], outdir=args[3], exp_hom=strtoi(args[4]))
