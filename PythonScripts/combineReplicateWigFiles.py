# -*- coding: utf-8 -*-
"""
Python SubFunction to normalize and then combine replicates.
Part of the ChIP-exo data analysis pipeline
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""

import numpy as np
import sys

dataList_split={}
for repI in [1,2]:
    #import TF data
    data=[[x for x in line.rstrip('\n\r').split(' ')] for line in open(sys.argv[repI])]
    data=data[1:]
    dataList=[]
    for i in range(0,len(data)):
        if data[i][0]=='variableStep':
            chromTMP=data[i][1][6:]
        else:
            dataList.append([chromTMP,data[i][0],data[i][1]])          
    #print('Done Import TF Data - '+datetime.now().strftime('%H:%M:%S'))
               
    #remove part on chrom 12, plasmid-2-micron, mit or unplaced telemore
    filterList=[]
    for i in range(0,len(dataList)):
        if dataList[i][0]=='chr12':
            if int(dataList[i][1])>=429000 and int(dataList[i][1])<=438500:
                filterList.append(i)
        elif dataList[i][0]=='plasmid-2-micron' or dataList[i][0]=='mit' or dataList[i][0]=='unplaced_telomere_1' or dataList[i][0]=='unplaced_telomere_2':
            filterList.append(i)
    for index in sorted(filterList, reverse=True):
        del dataList[index]
    del chromTMP, data, index, filterList
    #print('Done Chr12 Filtering - '+datetime.now().strftime('%H:%M:%S'))

    #Calculate Background Value
    TFsumTMP=0
    for i in range(0,len(dataList)):
        TFsumTMP+=float(dataList[i][2])
    
    #Divide by genome length and scale by 10
    TFbackground=TFsumTMP/12147764*10
    del TFsumTMP
    print('Background level for '+sys.argv[repI]+' : '+str(round(TFbackground,2)))
                     
    #Split Data into Chromosomes
    dataList_split[repI] = {}
    for chrI in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16']:
        dataList_split[repI][chrI] = {}    
    for item in dataList:
        dataList_split[repI][item[0]][int(item[1])]=round(float(item[2])/TFbackground,6)
    #print('Done Data splitting - '+datetime.now().strftime('%H:%M:%S'))   
            
#find all positions
combinedList={}
for chrI in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16']:
    combinedList[chrI]=[]
    allKeys=list(dataList_split[1][chrI].keys())+list(dataList_split[2][chrI].keys())
    for k in allKeys:
        combinedList[chrI].append([k,round(np.mean([dataList_split[1][chrI].get(k,0),dataList_split[2][chrI].get(k,0)]),6)])
        
            
#write processed wigFile
f = open(sys.argv[3], 'w')
f.write('track type=track1\n')
for chrI in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16']:
    if len(combinedList[chrI])>0:
        f.write('variableStep chrom='+chrI+'\n')
        f.write(str(combinedList[chrI][0][0])+' '+str(combinedList[chrI][0][1])+'\n')
        for dataI in range(1,len(combinedList[chrI])):
            if combinedList[chrI][dataI][0] == combinedList[chrI][dataI-1][0] + 1:
                f.write(str(combinedList[chrI][dataI][0])+' '+str(combinedList[chrI][dataI][1])+'\n')
            else:
                f.write('variableStep chrom='+chrI+'\n')
                f.write(str(combinedList[chrI][dataI][0])+' '+str(combinedList[chrI][dataI][1])+'\n')
            
f.close()  