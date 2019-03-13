#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python SubFunction for assigin TF read data to genes using published Transcription Start Site data (see Boerlin et. al., https://doi.org/10.1093/femsyr/foy128).
Part of the ChIP-exo data analysis pipeline
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""
import numpy as np
import sys
from collections import defaultdict

pathToTSS=sys.argv[1]

#import TSS annotations and take only dominat cluster
TSSData=[[x for x in line.rstrip('\n\r').split('\t')] for line in open(pathToTSS)][1:]
TSSData_genes=[x[0] for x in TSSData]

#Load data
data=[[x for x in line.rstrip('\n\r').split(' ')] for line in open(sys.argv[2])]
data=data[1:]
dataList=[]
for i in range(0,len(data)):
	if data[i][0]=='variableStep':
		chromTMP=data[i][1][6:]
	else:
		dataList.append([chromTMP,data[i][0],data[i][1]])          

#Split Data into Chromosomes
dataList_split = defaultdict(list)
for item in dataList:
	dataList_split[item[0]].append([int(item[1]),round(float(item[2]),6)])

#assign data to genes and create wiglike files
f = open(sys.argv[3], 'w')

for i in range(0,len(TSSData)):
	#check strand of gene and set range accordingly
	if TSSData[i][3]=='+':
		selPos=range(int(TSSData[i][2])-1000,int(TSSData[i][2])+1000)
	elif TSSData[i][3]=='-':
		selPos=range(int(TSSData[i][2])+1000,int(TSSData[i][2])-1000,-1)

	tmpData=np.array([0.0]*2000)
	#Get data for the selected chromosome
	chrData=dataList_split[TSSData[i][1]]
	k=-1
	while k+1 < len(chrData) and chrData[k+1][0] not in selPos:
		k+=1
	k+=1
	while k < len(chrData) and chrData[k][0] in selPos:
		f.write(TSSData[i][0]+'\t'+str(selPos.index(chrData[k][0]))+'\t'+str(chrData[k][1])+'\n')
		k+=1

f.close()