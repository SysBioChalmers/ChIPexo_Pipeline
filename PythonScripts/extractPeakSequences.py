# -*- coding: utf-8 -*-
"""
Python SubFunction for extracting the peak sequences for motif discovery
Part of the ChIP-exo data analysis pipeline
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""

import pandas as pd
import os.path
import sys
import csv

TF=sys.argv[1]
cond=sys.argv[2]
pathToGEM=sys.argv[3]
pathForOutput=sys.argv[4]

if os.path.isfile(pathToGEM):

    #Load GEM Data
    gemData=pd.read_csv(pathToGEM,sep='\t')
    
    #filter from chromosomes
    allowedChromosomes=[str(x) for x in range(1,17)]
    filterList=[]
    for i in range(0,len(gemData)):
        chrTMP,posTMP=gemData.loc[gemData.index[i],'Position'].split(':')
        if chrTMP in allowedChromosomes:
            filterList.append(i)
     
    gemData=gemData.iloc[filterList,:]  
    
    gemDataTMP=gemData.copy()
    gemDataTMP['Strength']=gemDataTMP.loc[:,'Cond'+cond+'_IP']/gemDataTMP.loc[:,'Cond'+cond+'_Expectd']
    
    #filter for strength
    gemDataTMP=gemDataTMP.loc[gemDataTMP['Strength']>2,['Position','Strength']]
     
    #write out peak positions as bed file +-50 bp around peak
    with open(pathForOutput, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for i in range(0,len(gemDataTMP)):
            chrTMP,posTMP=gemDataTMP.loc[gemDataTMP.index[i],'Position'].split(':')
            writer.writerow(['chr'+chrTMP,int(posTMP)-30,int(posTMP)+30])

else:
    print('GEM file does not exist')