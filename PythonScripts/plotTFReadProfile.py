#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python SubFunction for generating the TF read  profile figures.
Part of the ChIP-exo data analysis pipeline
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

selectedTF=sys.argv[1]
selectedCond=sys.argv[2]
pathToTF=sys.argv[3]
outputPath=sys.argv[4]
date=sys.argv[5]

#Load raw data and get targeted genes
TFDataRaw=[[x for x in line.rstrip('\n\r').split('\t')] for line in open(pathToTF)]
TFDataGenes=list(set([x[0] for x in TFDataRaw]))

#assign reads into dataframe for plotting
TFData=pd.DataFrame(0.0,index=TFDataGenes,columns=list(range(0,2000)))
for t in TFDataRaw:
    TFData.loc[t[0],int(t[1])]=float(t[2])
      
#plot
fig = plt.figure(figsize=(10,6),dpi=80)
plt.plot(TFData.columns,[1]*len(TFData.columns),'--',color='grey',lineWidth=2,label='Average')
plt.plot(TFData.mean(axis=0)/TFData.mean(axis=0).mean(),'g',lineWidth=2,label=selectedTF+'_'+selectedCond)
plt.xticks([0,500,1000,1500,2000],['-1000bp','-500bp','TSS','+500bp','+1000bp'])
plt.ylabel('Normalized Read Count',fontsize=16)
plt.legend(fontsize=16)
plt.title('Read distribution around TSS for '+selectedTF+' in '+selectedCond,fontsize=16)
fig.savefig(outputPath+'/'+selectedTF+'_'+selectedCond+'_ReadProfile_'+date+'.png',dpi=300,bbox_inches="tight")