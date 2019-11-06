#########################################################################
# File Name: getData.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Fri 16 Aug 2019 15:22:29 AEST
#########################################################################
#!/bin/bash

import os
import sys

indel_sen = []
indel_pre = []
snp_sen = []
snp_pre = []
with open(sys.argv[1]) as f:
    count = 1
    for line in f:
        words = line.split()
        if count%3 == 1:
            iso_true = int(words[0])
            iso_find = int(words[1])
            iso_TP = int(words[2])
        if count%3 == 2:
            non_true = int(words[0])*2
            non_find = int(words[1])*2
            non_TP = int(words[2])*2
            snp_TP = non_TP + iso_TP
            snp_true = non_true + iso_true
            snp_find = non_find + iso_find
            snp_sen.append(round(float(snp_TP)/snp_true*100.0, 2))

            snp_pre.append(round(float(snp_TP)/snp_find*100.0, 2))
        if count%3 == 0:
            indel_sen.append(float(words[4]))
            indel_pre.append(float(words[5]))
        count += 1 

print snp_sen[0], "&", snp_pre[0], "&", snp_TP, "&", snp_find-snp_TP, "&", snp_true-snp_TP
print snp_sen[0], "&", snp_pre[0], "&", round(2*snp_sen[0]*snp_pre[0]/(snp_sen[0] + snp_pre[0]), 2)

f1 = []
l = len(snp_sen)
for i in range(l):
    temp = round(2*snp_sen[i]*snp_pre[i]/(snp_sen[i] + snp_pre[i]), 2)
    f1.append(temp)

print snp_sen
print snp_pre
print f1


'''
for ele in snp_sen:
    print ele,
print 
for ele in snp_pre:
    print ele,
'''    
#print indel_sen
#print indel_pre
