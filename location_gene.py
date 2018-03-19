#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 16:33:11 2018

@author: kkwang
"""

import numpy as np
import pandas as pd 

gene_loc_info=pd.read_table('mart_export.txt',header=0)

gene_loc_info.columns=['id','chr','start','end','gene']
#gene_info.index=gene_info['gene']
#del gene_info['gene']

row,col=np.where(VFC_bool==True)
gene1=[differ_gene[x] for x in row]
gene2=[differ_gene[x] for x in col]
d1={'gene':gene1}
d2={'gene':gene2}
df1=pd.DataFrame(d1)
#DF1 = df1.join(gene_info, on='gene')
df1=pd.merge(df1, gene_info, how='left', on=['gene'])

df1['chr']='chr'+df1['chr']
df2=pd.DataFrame(d2)
#DF2=df2.join(gene_info, on='gene')
df2=pd.merge(df2, gene_info, how='left', on=['gene'])
df2['chr']='chr'+df2['chr']


circos_data=pd.concat([df1,df2],axis=1)

circos=circos_data.iloc[:,[2,3,4,7,8,9]]

#circos.columns=['chromosomenameA', 'chromStartA', 'chromEndA', 'chromosomenameB', 'chromStartB', 'chromEndB']
hhhhh=circos.dropna()
hhhhh.to_csv('circos.tsv',sep='\t')        
