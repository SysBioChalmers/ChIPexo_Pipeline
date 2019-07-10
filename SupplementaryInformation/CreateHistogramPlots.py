#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python SubFunction for creating the HeatMap Data
Part of the Supplementary Information
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

selectedTFandConds={'Ino2':'Glu','Stb5':'Nit','Gcn4':'Glu','Cbf1':'Ana'}

path='C:/Users/borlinc/Documents/Projects/190219_ChipExoPipeline/SupplementaryInformation/'
pathToGEM=path+'GEM_Files/'
pathToTF=path+'ReadData/'

for TF,cond in selectedTFandConds.items():   
    
    #load histogram data
    histogramData=pd.read_csv(path+TF+'_'+cond+'_PeakHistogramData.csv',sep='\t',index_col=0)
    
    #sort peaks after read count
    sortIndex=histogramData.loc[:,[str(x)+'_plus' for x in range(-50,51)]+[str(x)+'_minus' for x in range(-50,51)]].sum(axis=1).sort_values(0,ascending=False).index 
    histogramData=histogramData.reindex(sortIndex,axis='index')
    
    #extract strandProfiles   
    strandProfile={}
    strandProfile['plus']=histogramData.loc[:,[str(x)+'_plus' for x in range(-50,51)]]
    strandProfile['minus']=histogramData.loc[:,[str(x)+'_minus' for x in range(-50,51)]]
     
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
    plt.title('Read distribution around '+str(len(histogramData))+' peaks for '+TF+' in '+cond,fontsize=16)
    ax.imshow(blueData,interpolation=None,aspect='auto')
    ax.imshow(redData,interpolation=None,aspect='auto')
    
    fig.savefig(path+TF+'_'+cond+'_PeakHistogram_ReadCountSorted.png',dpi=300,bbox_inches="tight")
    
    
    #sort peaks after peakStrength
    sortIndex=histogramData.loc[:,'Peak Strength'].sort_values(0,ascending=False).index 
    histogramData=histogramData.reindex(sortIndex,axis='index')
    
    #extract strandProfiles   
    strandProfile={}
    strandProfile['plus']=histogramData.loc[:,[str(x)+'_plus' for x in range(-50,51)]]
    strandProfile['minus']=histogramData.loc[:,[str(x)+'_minus' for x in range(-50,51)]]
     
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
    plt.title('Read distribution around '+str(len(histogramData))+' peaks for '+TF+' in '+cond,fontsize=16)
    ax.imshow(blueData,interpolation=None,aspect='auto')
    ax.imshow(redData,interpolation=None,aspect='auto')
    
    fig.savefig(path+TF+'_'+cond+'_PeakHistogram_PeakStrengthSorted.png',dpi=300,bbox_inches="tight")
