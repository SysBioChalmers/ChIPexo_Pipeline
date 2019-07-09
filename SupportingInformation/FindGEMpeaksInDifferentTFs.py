# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:09:59 2019

@author: borlinc
"""

import pandas as pd
import matplotlib.pyplot as plt

TFList=['Cat8','Cbf1','Ert1','Gcn4','Gcr1','Gcr2','Hap1','Ino2','Ino4','Leu3','Oaf1','Pip2','Rds2','Rgt1','Rtg1','Rtg3','Sip4','Stb5','Sut1','Tye7']
cond='Glu'

path='C:/Users/borlinc/Documents/Projects/190219_ChipExoPipeline/SupportingInformation/'
pathToGEM=path+'GEM_Files/'


gemPeaks={}
for i in range(1,17):
    gemPeaks[i]=[]
    
allowedChromosomes=[str(x) for x in range(1,17)]

numTotalPeaks=0
peakCounts=[]

for tf in TFList:
    #load GEM data and get only selected condition
    gemData=pd.read_csv(pathToGEM+tf+'_GEM.GEM_events.txt',sep='\t')
    #take only peaks with SNR >2
    gemData=gemData.loc[gemData.loc[:,'CondGlu_IP']/gemData.loc[:,'CondGlu_Expectd']>2,'Position']
    peakCounts.append(len(gemData))
    
    
    for i in range(0,len(gemData)):
        chrTMP,posTMP=gemData.iloc[i].split(':')
        if chrTMP in allowedChromosomes:
            gemPeaks[int(chrTMP)].append(int(posTMP))
            numTotalPeaks+=1


#sort peaks and make unique ones    
for chrTMP in gemPeaks.keys():
    gemPeaks[chrTMP]=sorted(list(set(gemPeaks[chrTMP])),key=int)
    
    
#do clustering using greedy algorithm
#create intital cluters
gemClusterTMP={}
for chrTMP in gemPeaks.keys():
    gemClusterTMP[chrTMP]=[]
    for chrPos in gemPeaks[chrTMP]:
        gemClusterTMP[chrTMP].append(chrPos)    

gemCluster={}
for chrTMP in gemPeaks.keys():
    gemCluster[chrTMP]={}
   
minDistance=50
    
for chrTMP in gemClusterTMP.keys():
    #do first cluster
    tempI=0
    gemCluster[chrTMP][gemClusterTMP[chrTMP][0]]=[gemClusterTMP[chrTMP][0]]
        
    #check for cluster 1 to n-1    
    for i in range(1,len(gemClusterTMP[chrTMP])-1):
        if gemClusterTMP[chrTMP][i]-gemClusterTMP[chrTMP][i-1]>minDistance:
            tempI=i
            gemCluster[chrTMP][gemClusterTMP[chrTMP][i]]=[gemClusterTMP[chrTMP][i]]
        else:
            gemCluster[chrTMP][gemClusterTMP[chrTMP][tempI]].append(gemClusterTMP[chrTMP][i])
            
    #check for last cluster
    if gemClusterTMP[chrTMP][-1]-gemClusterTMP[chrTMP][-2]>minDistance:
         gemCluster[chrTMP][gemClusterTMP[chrTMP][-1]]=[gemClusterTMP[chrTMP][-1]]
    else:
        gemCluster[chrTMP][gemClusterTMP[chrTMP][tempI]].append(gemClusterTMP[chrTMP][-1])
        
         
#STEP 2: GO through all 
gemClusterPresent={}
for chrTMP in gemCluster.keys():
    for posTMP in gemCluster[chrTMP].keys():
        gemClusterPresent[str(chrTMP)+':'+str(posTMP)]=0


for tf in TFList:
    #load GEM data and get only selected condition
    gemData=pd.read_csv(pathToGEM+tf+'_GEM.GEM_events.txt',sep='\t')
    #take only peaks with SNR >2
    gemData=gemData.loc[gemData.loc[:,'CondGlu_IP']/gemData.loc[:,'CondGlu_Expectd']>2,'Position']
    
    for chrTMP in gemCluster.keys():
        for posTMP in gemCluster[chrTMP].keys():
            if len([x for x in gemCluster[chrTMP][posTMP] if str(chrTMP)+':'+str(x) in gemData.tolist()])>0:
                gemClusterPresent[str(chrTMP)+':'+str(posTMP)]+=1
    

#STEP 3: Plotting
x=[]
y=[]
for i in range(1,len(TFList)+1):
    x.append(i)
    y.append(len([x for x in gemClusterPresent.values() if x ==i]))

plt.figure(figsize=(8,4))
#plt.plot(x,[0]*len(x),'--r')
plt.plot(x,y,'x')
plt.xticks(x)
plt.xlabel('Found in x independent ChIP experiments',fontSize=18)
plt.ylabel('Number of Peaks',fontSize=18)
plt.savefig(path+'GEMpeaksInDifferentTFs.png',dpi=300)


#STEP 4: Analyse
clusterSpan=[]
for chrTMP in gemCluster.keys():
    for c in gemCluster[chrTMP].values():
        if len(c)>1:
            clusterSpan.append(c[-1]-c[0])
            
print('Longest cluster is '+str(max(clusterSpan))+' bp long')
