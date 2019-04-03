#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python SubFunction for generating the sample correlation figures.
Part of the ChIP-exo data analysis pipeline
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""


import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import sys

selectedTF=sys.argv[1]
pathToTF=sys.argv[2]
outputPath=sys.argv[3]
pathToChromSizes=sys.argv[4]
date=sys.argv[5]
samples=sys.argv[6:]

#get number of samples, calculate number of pairwise comparisions
numSamples=len(samples)
numPlots=int(numSamples*(numSamples-1)/2)

#get plot and label positions
plotPositions=[]
plotPartners=[]
xlabels={}
ylabels={}
for rowCount in range(1,numSamples):
	for columnCount in range(0,rowCount):
		plotPositions.append((rowCount-1)*(numSamples-1)+columnCount+1)
		plotPartners.append([columnCount,rowCount])
		if rowCount==numSamples-1:
			xlabels[(rowCount-1)*(numSamples-1)+columnCount+1]=samples[columnCount]
		if columnCount==0:
			ylabels[(rowCount-1)*(numSamples-1)+columnCount+1]=samples[rowCount]

#load genome sizes
genomeSizeTMP=[[x for x in line.rstrip('\n\r').split('\t')] for line in open(pathToChromSizes)]
genomeSize={}
for chrTMP, n in genomeSizeTMP:
	if chrTMP[0]=='c':
		genomeSize[chrTMP]=int(n)
del genomeSizeTMP

#load TF data and convert it to plotting data
sampleData={}
sampleDataBinned={}
for sample in samples:
	TFDataRaw=[[x for x in line.rstrip('\n\r').split(' ')] for line in open(pathToTF.replace('SAMPLE',sample))]
	TFDataRaw=TFDataRaw[1:]
	TFDataTMP=[]
	for i in range(0,len(TFDataRaw)):
		if TFDataRaw[i][0]=='variableStep':
			chromTMP=TFDataRaw[i][1][6:]
		else:
			TFDataTMP.append([chromTMP,TFDataRaw[i][0],TFDataRaw[i][1]])          

	sampleData[sample]={}
	for chrI in sorted(list(set([x[0] for x in TFDataTMP]))):
		sampleData[sample][chrI]=[]
	for item in TFDataTMP:
		sampleData[sample][item[0]].append([int(item[1]),float(item[2])])
	del TFDataRaw, TFDataTMP

	#Sum up every 1kb interval for correlation
	intervalSize=1000
	sampleDataBinned[sample]=[]
	for chrTMP in genomeSize.keys():
		for i in range(0,genomeSize[chrTMP],intervalSize):
			sampleDataBinned[sample].append(sum([x[1] for x in sampleData[sample][chrTMP] if x[0] >=i and x[0] < min(i+intervalSize,genomeSize[chrTMP])]))

#go through all positions and only take them if there are reads in at least one sample
sampleDataBinnedFiltered={}
for sample in sampleDataBinned.keys():
	sampleDataBinnedFiltered[sample]=[]
for i in range(0,len(sampleDataBinned[samples[0]])):
	if sum([sampleDataBinned[x][i] for x in sampleDataBinned.keys()])>0:
		for sample in sampleDataBinned.keys():
			sampleDataBinnedFiltered[sample].append(sampleDataBinned[sample][i])

#make plot
fig=plt.figure(figsize=(4*(numSamples-1),4*(numSamples-1)))
ax={}

for i in range(0,len(plotPositions)):
	ax[i]=fig.add_subplot(numSamples-1,numSamples-1,plotPositions[i])
	with np.errstate(divide='ignore'):
		sample1=np.log2(sampleDataBinnedFiltered[samples[plotPartners[i][0]]])
		sample1[np.isinf(sample1)]=0
		sample2=np.log2(sampleDataBinnedFiltered[samples[plotPartners[i][1]]])
		sample2[np.isinf(sample2)]=0
	ax[i].plot(sample1,sample2,'bo',ms=5,alpha=0.3,mew=0)
	PCC=round(pearsonr(sampleDataBinnedFiltered[samples[plotPartners[i][0]]],sampleDataBinnedFiltered[samples[plotPartners[i][1]]])[0],2)
	plt.title(samples[plotPartners[i][0]]+' vs ' +samples[plotPartners[i][1]]+' - PCC: '+str(PCC),fontSize=16)
	if xlabels.get(plotPositions[i],0)!=0:
		plt.xlabel(xlabels[plotPositions[i]]+' read count [log2]',fontSize=16)
	if ylabels.get(plotPositions[i],0)!=0:
		plt.ylabel(ylabels[plotPositions[i]]+' read count [log2]',fontSize=16)     

fig.savefig(outputPath+'/'+selectedTF+'_PairwiseSampleComparision_'+date+'.png',dpi=300,bbox_inches="tight")