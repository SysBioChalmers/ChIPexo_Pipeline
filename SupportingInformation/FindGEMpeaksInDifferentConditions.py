# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:09:59 2019

@author: borlinc
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

TFList=['Cat8','Cbf1','Ert1','Gcn4','Gcr1','Gcr2','Hap1','Ino2','Ino4','Leu3','Oaf1','Pip2','Rds2','Rgt1','Rtg1','Rtg3','Sip4','Stb5','Sut1','Tye7']
allowedChromosomes=[str(x) for x in range(1,17)]


path='C:/Users/borlinc/Documents/Projects/190219_ChipExoPipeline/SupportingInformation/'
pathToGEM=path+'GEM_Files/'

conditionBasis='Glu'
conditionOther=['Ana','Eth','Nit']

selectedTFs=['Ino2','Cbf1','Gcn4','Stb5']

distancesClosestPeakSelectedTFOtherConditions={}
distancesClosestPeakOtherTFs={}

filteredPeaks={x:0 for x in selectedTFs}
numPeaks={x:0 for x in selectedTFs}

for selectedTF in selectedTFs:
    otherTFs=[x for x in TFList if x != selectedTF]

    #STEP 1: Get peaks for selected TF
    gemDataSelectedTF=pd.read_csv(pathToGEM+selectedTF+'_GEM.GEM_events.txt',sep='\t')
    
    #Split into basis and other condition
    #take only peaks with SNR >2 in basis Condition
    gemDataSelectedBasis=gemDataSelectedTF.loc[gemDataSelectedTF.loc[:,'Cond'+conditionBasis+'_IP']/gemDataSelectedTF.loc[:,'Cond'+conditionBasis+'_Expectd']>2,'Position']
    
    #count peaks
    numPeaks[selectedTF]=len(gemDataSelectedTF)
    
    gemPeaksSelectedTFBasis={}
    for i in range(1,17):
        gemPeaksSelectedTFBasis[i]=[]
    
    for i in range(0,len(gemDataSelectedBasis)):
        chrTMP,posTMP=gemDataSelectedBasis.iloc[i].split(':')
        if chrTMP in allowedChromosomes:
            gemPeaksSelectedTFBasis[int(chrTMP)].append(int(posTMP))

    #take only peaks with SNR >2 in the other conditions
    
    selectionIndex=[]
    for cond in conditionOther:
        selectionIndex+=gemDataSelectedTF.index[gemDataSelectedTF.loc[:,'Cond'+cond+'_IP']/gemDataSelectedTF.loc[:,'Cond'+cond+'_Expectd']>2].tolist()
        
    selectionIndex=list(set(selectionIndex))
    gemDataSelectedOther=gemDataSelectedTF.loc[selectionIndex,'Position']

    gemPeaksSelectedTFOther={}
    for i in range(1,17):
        gemPeaksSelectedTFOther[i]=[]

    for i in range(0,len(gemDataSelectedOther)):
        chrTMP,posTMP=gemDataSelectedOther.iloc[i].split(':')
        if chrTMP in allowedChromosomes:
            gemPeaksSelectedTFOther[int(chrTMP)].append(int(posTMP))
    

    #STEP 2: Load GEM data for the other TFs in basis condition 
    gemPeaksOtherTFs={}
    for i in range(1,17):
        gemPeaksOtherTFs[i]=[] 
            
    for otherTF in otherTFs:
        #load GEM data and get only selected condition
        gemDataOtherTF=pd.read_csv(pathToGEM+otherTF+'_GEM.GEM_events.txt',sep='\t')
        #take only peaks with SNR >2
        gemDataOtherTF=gemDataOtherTF.loc[gemDataOtherTF.loc[:,'Cond'+conditionBasis+'_IP']/gemDataOtherTF.loc[:,'Cond'+conditionBasis+'_Expectd']>2,'Position']
            
        for i in range(0,len(gemDataOtherTF)):
            chrTMP,posTMP=gemDataOtherTF.iloc[i].split(':')
            if chrTMP in allowedChromosomes:
                gemPeaksOtherTFs[int(chrTMP)].append(int(posTMP))
    
    
    #sort peaks and make unique ones    
    for chrTMP in gemPeaksOtherTFs.keys():
        gemPeaksOtherTFs[chrTMP]=sorted(list(set(gemPeaksOtherTFs[chrTMP])),key=int)
        
        
    #STEP 3: Go through each detected peak in the basis condition of the selected TF and find the closest Peak in the selected TF in the other conditions as well as in the other TFs
    distancesClosestPeakSelectedTFOtherConditions[selectedTF]=[]
    distancesClosestPeakOtherTFs[selectedTF]=[]
    for chrTMP in gemPeaksSelectedTFBasis.keys():
        for posTMP in gemPeaksSelectedTFBasis[chrTMP]:
            distanceTMP=min([abs(x-posTMP) for x in gemPeaksSelectedTFOther[chrTMP]])
            #Check if there is a peak within 1kb. Otherwise the target gene is just not targeted in any of the other conditions
            if distanceTMP<=1000:
                distancesClosestPeakSelectedTFOtherConditions[selectedTF].append(distanceTMP)
                distancesClosestPeakOtherTFs[selectedTF].append(min([abs(x-posTMP) for x in gemPeaksOtherTFs[chrTMP]]))
            else:
                filteredPeaks[selectedTF]+=1
                
            
    

#STEP 4: Plot
boxPlotData=[]
boxPlotLabels=[]
boxPlotPValues=[]
for selectedTF in selectedTFs:
    boxPlotData.append(distancesClosestPeakSelectedTFOtherConditions[selectedTF])
    boxPlotLabels.append(selectedTF+'-otherCond')
    boxPlotData.append(distancesClosestPeakOtherTFs[selectedTF])
    boxPlotLabels.append(selectedTF+'-otherTF')
    
    #Perform two-sided TTest
    print(str(filteredPeaks[selectedTF])+' ('+str(round(filteredPeaks[selectedTF]/numPeaks[selectedTF]*100,2))+'%) Peaks have been filtered out')
    ttest=ttest_ind(distancesClosestPeakSelectedTFOtherConditions[selectedTF],distancesClosestPeakOtherTFs[selectedTF])
    if ttest.statistic<0:
        print('One sided TTest for '+selectedTF+' gives pValue of '+str(ttest.pvalue/2))
        boxPlotPValues.append("{0:.1e}".format(ttest.pvalue/2))
    else:
        print('One sided TTest for '+selectedTF+' shows that Other is smaller')
        boxPlotPValues.append('NA')
    
    
    
plt.figure(figsize=(16,6))
axes = plt.gca()
plt.boxplot(boxPlotData,labels=boxPlotLabels)
plt.ylabel('Peak Distance',fontSize=16)
plt.xticks(fontsize=14)
ylim=axes.get_ylim()
plt.ylim([ylim[0],ylim[1]*1.2])
for i in range(0,len(boxPlotPValues)):
    plt.plot([2*(i+1)-1,2*(i+1)-1],[ylim[1]*1.05,ylim[1]*1.1],color='black',linewidth=1.5)
    plt.plot([2*(i+1),2*(i+1)],[ylim[1]*1.05,ylim[1]*1.1],color='black',linewidth=1.5)
    plt.plot([2*(i+1)-1,2*(i+1)],[ylim[1]*1.1,ylim[1]*1.1],color='black',linewidth=1.5)
    plt.text(2*(i+1)-0.7,ylim[1]*1.12,boxPlotPValues[i],fontsize=15)

plt.savefig(path+'GEMpeaksInDifferentConditions.png',dpi=300)
