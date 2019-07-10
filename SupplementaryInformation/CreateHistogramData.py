#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python SubFunction for exporting the HeatMap Data
Part of the Supplementary Information
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""
import pandas as pd

selectedTFandConds={'Ino2':'Glu','Stb5':'Nit','Gcn4':'Glu','Cbf1':'Ana'}

path='C:/Users/borlinc/Documents/Projects/190219_ChipExoPipeline/SupplementaryInformation/'
pathToGEM=path+'GEM_Files/'
pathToTF=path+'ReadData/'

for TF,cond in selectedTFandConds.items():
    TFData={}
    for strand in ['plus','minus']:
    	#import TF data and split into Chromosomes
    	TFDataRaw=[[x for x in line.rstrip('\n\r').split(' ')] for line in open(pathToTF+TF+'_'+cond+'_'+strand+'_singlePos_combRep.wig')]
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
    gemData=pd.read_csv(pathToGEM+TF+'_GEM.GEM_events.txt',sep='\t')
    #Filter out peaks not on Chr1 to Chr16
    filterIndex=[]
    for i in range(0,len(gemData)):
        chrTMP,posTMP=gemData.loc[:,'Position'][i].split(':')
        if chrTMP in [str(x) for x in range(0,17)]:
            filterIndex.append(i)
    gemData=gemData.iloc[filterIndex]
    #select only peaks present int that condition
    gemData=gemData.loc[gemData.loc[:,'Cond'+cond+'_IP']/gemData.loc[:,'Cond'+cond+'_Expectd']>2,['Position','Cond'+cond+'_IP','Cond'+cond+'_Expectd']]
    gemData['Strength']=gemData.loc[:,'Cond'+cond+'_IP']/gemData.loc[:,'Cond'+cond+'_Expectd']
      
    #creat output table
    outputTable=pd.DataFrame(0,index=gemData.index,columns=['Peak Strength']+[str(x)+'_plus' for x in range(-50,51)]+[str(x)+'_minus' for x in range(-50,51)])
         
    for peakID in gemData.index:
        outputTable.loc[peakID,'Peak Strength']=gemData.loc[peakID,'Strength']
        selChr,selPos=gemData.loc[peakID,'Position'].split(':')
        selChr='chr'+selChr
        selPos=int(selPos)
        for strand in ['plus','minus']:
    		#find position and save the values
            posList=[x for x in TFData[strand][selChr].keys() if x>=selPos-50 and x<=selPos+50]
            for pos in posList:
                outputTable.loc[peakID,str(pos-selPos)+'_'+strand]=TFData[strand][selChr][pos]
   
    #sort peaks after read count
    sortIndex=outputTable.loc[:,[str(x)+'_plus' for x in range(-50,51)]+[str(x)+'_minus' for x in range(-50,51)]].sum(axis=1).sort_values(0,ascending=False).index
    
    outputTable=outputTable.reindex(sortIndex,axis='index')
 
    #save table
    outputTable.to_csv(path+TF+'_'+cond+'_PeakHistogramData.csv',sep='\t')