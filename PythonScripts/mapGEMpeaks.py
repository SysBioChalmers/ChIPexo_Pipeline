#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python SubFunction for mapping GEM peaks to genes using published Transcription Start Site data (see Boerlin et. al., https://doi.org/10.1093/femsyr/foy128).
Part of the ChIP-exo data analysis pipeline
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""

import pandas as pd
import os.path
import sys
import re

tf=sys.argv[1]
pathToGEM=sys.argv[2]
pathToTSS=sys.argv[3]
pathForOutput=sys.argv[4]
date=sys.argv[5]

#import TSS annotations and take only dominat cluster
TSSData=[[x for x in line.rstrip('\n\r').split('\t')] for line in open(pathToTSS)][1:]
TSSData_genes=[x[0] for x in TSSData]

if os.path.isfile(pathToGEM):
    #load GEM data
    gemData=pd.read_csv(pathToGEM,sep='\t')
    #get conditions
    condListGEM=[x for x in gemData.columns.tolist() if bool(re.match('Cond[a-zA-Z0-9]*_Present',x))]
    condList=[x.replace('Cond','').replace('_Present','') for x in condListGEM]
    #gemData=gemData.loc[:,['Position']+condListGEM]
    #create outputTable
    outputData=pd.DataFrame(0,index=TSSData_genes,columns=condList)
    outputDataSNR=pd.DataFrame(0,index=TSSData_genes,columns=condList) 
    outputDataDetailed=pd.DataFrame(0,index=[],columns=['Gene','Chr','Pos','Strand','DistanceTSS','Strength','Condition'])
    #iterate through GEM peaks, assign to genes
    for i in range(0,len(gemData)):
        selChr='chr'+str(gemData.iloc[i]['Position'].split(':')[0])
        selPos=int(gemData.iloc[i]['Position'].split(':')[1])
        #find gene
        for gene in TSSData:
            if selChr==gene[1] and abs(selPos-int(gene[2]))<=1000:
                for cond in condList:
                    if gemData.iloc[i]['Cond'+cond+'_IP']/gemData.iloc[i]['Cond'+cond+'_Expectd']>2:
                        outputData.loc[gene[0],cond]+=gemData.iloc[i]['Cond'+cond+'_Present']
                        outputDataSNR.loc[gene[0],cond]+=gemData.iloc[i]['Cond'+cond+'_IP']/gemData.iloc[i]['Cond'+cond+'_Expectd']
                        if gene[3]=='+':
                            outputDataDetailed.loc[len(outputDataDetailed)+1]=[gene[0],selChr,selPos,gene[3],selPos-int(gene[2]),gemData.iloc[i]['Cond'+cond+'_IP']/gemData.iloc[i]['Cond'+cond+'_Expectd'],cond]
                        else:
                            outputDataDetailed.loc[len(outputDataDetailed)+1]=[gene[0],selChr,selPos,gene[3],int(gene[2])-selPos,gemData.iloc[i]['Cond'+cond+'_IP']/gemData.iloc[i]['Cond'+cond+'_Expectd'],cond]

    outputData=outputData.loc[outputData.sum(axis=1)>0]
    outputDataSNR=outputDataSNR.loc[outputDataSNR.sum(axis=1)>0]
    outputDataSNR.to_csv(pathForOutput+'/'+tf+'_geneTargetListSNR_'+date+'.csv')
    outputData.to_csv(pathForOutput+'/'+tf+'_geneTargetList_'+date+'.csv')
    outputDataDetailed.to_csv(pathForOutput+'/'+tf+'_GEManalysis_'+date+'.csv')
else:
    print('GEM file does not exist')
