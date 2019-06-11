#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python SubFunction for generating the peak centered read profile figures.
Part of the ChIP-exo data analysis pipeline
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

selectedTF=sys.argv[1]
cond=sys.argv[2]
pathToGEM=sys.argv[3]
date=sys.argv[6]
pathToTF={'plus':sys.argv[4].replace('STRAND','plus'),'minus':sys.argv[4].replace('STRAND','minus')}
outputPath=sys.argv[5].replace('COND',cond)

strandProfile={}
TFData={}
for strand in ['plus','minus']:
	strandProfile[strand]=pd.DataFrame(0,index=[],columns=range(-50,51))
	#import TF data and split into Chromosomes
	TFDataRaw=[[x for x in line.rstrip('\n\r').split(' ')] for line in open(pathToTF[strand])]
	TFDataRaw=TFDataRaw[1:]
	TFDataTMP=[]
	for i in range(0,len(TFDataRaw)):
		if TFDataRaw[i][0]=='variableStep':
			chromTMP=TFDataRaw[i][1][6:]
		else:
			TFDataTMP.append([chromTMP,TFDataRaw[i][0],TFDataRaw[i][1]])

	TFData[strand]={}
	for chrI in sorted(list(set([x[0] for x in TFDataTMP]))):
		TFData[strand][chrI]={}
	for item in TFDataTMP:
		TFData[strand][item[0]][int(item[1])]=float(item[2])
	del TFDataRaw, TFDataTMP

#Load GEM Data
gemData=pd.read_csv(pathToGEM,sep='\t')
#select only peaks present int that condition
gemData=gemData.loc[gemData.loc[:,'Cond'+cond+'_IP']/gemData.loc[:,'Cond'+cond+'_Expectd']>2,'Position']
if len(gemData)>0:
    gemData=[['chr'+x[0],int(x[1])] for x in [[x for x in line.split(':')] for line in gemData] if 'chr'+x[0] in TFData['plus'].keys()]
    
    #iterate through GEM peaks
    peakCounter=1
    for selChr,selPos in gemData:
    	for strand in ['plus','minus']:
    		#find position and save the values
    		posList=[x for x in TFData[strand][selChr].keys() if x>=selPos-50 and x<=selPos+50]
    		strandProfileTMP=pd.DataFrame(0,index=[peakCounter],columns=range(-50,51))
    		for pos in posList:
    			strandProfileTMP.loc[peakCounter,pos-selPos]=TFData[strand][selChr][pos]
    		strandProfile[strand]=strandProfile[strand].append(strandProfileTMP)
    	peakCounter+=1
    
    #sort peaks after read count
    sortIndex=(strandProfile['plus']+strandProfile['minus']).sum(axis=1).sort_values(0,ascending=False).index
    strandProfile['plus']=strandProfile['plus'].reindex(sortIndex,axis='index')
    strandProfile['minus']=strandProfile['minus'].reindex(sortIndex,axis='index')
    
    #create plot data for plus strand / blue color
    blueData=np.zeros([strandProfile['plus'].shape[0],strandProfile['plus'].shape[1],4])
    blueData[:,:,2]=np.ones(blueData[:,:,2].shape)
    blueData[:,:,3]=strandProfile['plus']/(np.max(strandProfile['plus'].values)*0.10)
    #everything above 0.75 is clipped to 0.75
    blueData[blueData[:,:,3]>0.75,3]=0.75
    
    #create plot data for minus strand / red color
    redData=np.zeros([strandProfile['minus'].shape[0],strandProfile['minus'].shape[1],4])
    redData[:,:,0]=np.ones(redData[:,:,2].shape)
    redData[:,:,3]=strandProfile['minus']/(np.max(strandProfile['minus'].values)*0.10)
    #everything above 0.75 is clipped to 0.75
    redData[redData[:,:,3]>0.75,3]=0.75
    
    #plot histogram
    fig, ax = plt.subplots(figsize=(10, 12))
    plt.xticks([0,24,50,75,100],[-50,-25,0,25,50])
    plt.xlabel('Distance from Peak [bp]',fontsize=16)
    plt.yticks([],[])
    plt.title('Read distribution around '+str(len(gemData))+' peaks for '+selectedTF+' in '+cond,fontsize=16)
    ax.imshow(blueData,interpolation=None,aspect='auto')
    ax.imshow(redData,interpolation=None,aspect='auto')
    fig.savefig(outputPath+'_PeakHistogram_'+date+'.png',dpi=300,bbox_inches="tight")
    
    #Plot overview
    maxValue=round(max([max(strandProfile['plus'].sum(axis=0)),max(strandProfile['minus'].sum(axis=0))])*1.05)
    fig = plt.figure(figsize=(10,6),dpi=80)
    plt.plot(range(-50,51),strandProfile['plus'].sum(axis=0),color='blue',label='Forward strand')
    plt.fill_between(range(-50,51), strandProfile['plus'].sum(axis=0), y2=0,color='blue',alpha=0.25)
    plt.plot(range(-50,51),-strandProfile['minus'].sum(axis=0),color='red',label='Reverse strand')
    plt.fill_between(range(-50,51), -strandProfile['minus'].sum(axis=0), y2=0,color='red',alpha=0.25)
    plt.ylim(-maxValue,maxValue)
    plt.ylabel('Read count',fontsize=16)
    plt.xlabel('Distance from Peak [bp]',fontsize=16)
    plt.title('Read distribution around '+str(len(gemData))+' peaks for '+selectedTF+' in '+cond,fontsize=16)
    plt.xticks([-50,-25,0,25,50],[-50,-25,0,25,50])
    plt.legend(fontsize=16,loc=1)
    fig.savefig(outputPath+'_PeakProfile_'+date+'.png',dpi=300,bbox_inches="tight")
else:
    print('No detected Peaks for '+cond)
