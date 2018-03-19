#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 19:37:40 2018

find TF gene in tissue 

and find correlation of TF and differ_genes(differ in tissue but not common differexpresssion)


 @author: kkwang
"""


import pandas as pd 
TF_gene=pd.read_csv('TFCheckpoint_download_180515.txt',sep='\t')
TF_gene=TF_gene[TF_gene['DbTF']=='yes']
VFC_data=pd.read_csv('VFC_log2_fiter_matrix.tsv',sep='\t',header=0)
VFC_TF=VFC_data[VFC_data['gene_symbol'].isin(TF_gene['gene_symbol'])]
differ_data_VFC=pd.read_csv('VFC_limma.tsv',header=0,sep='\t')
VFC_diff=VFC_data[VFC_data['gene_symbol'].isin(differ_data_VFC.index)]





from scipy.stats.stats import pearsonr
VFC_corr=pd.DataFrame(index=VFC_TF['gene_symbol'],columns=VFC_diff['gene_symbol'])

for ind,gtf in enumerate(VFC_TF['gene_symbol'].values):
    for index,gd in enumerate(VFC_diff['gene_symbol'].values):
        pcc,p_val=pearsonr(VFC_TF.iloc[ind,:-1].values, VFC_diff.iloc[index,:-1].values)
        VFC_corr.loc[gtf,gd]=pcc
        
                
import numpy as np 

VFC_corr[VFC_corr==1]=0.99
VFC_corr[VFC_corr>0.9]=1
VFC_corr[VFC_corr<-0.9]=1



with open('VFC_DFC.txt','r') as f:
    VFC_DFC=f.readlines()
    VFC_DFC=[gene[:-1] for gene in VFC_DFC]

VFC_corr_T=VFC_corr.T
VFC_DFC_corr=VFC_corr_T[VFC_corr_T.index.isin(VFC_DFC)]
   
row,col=np.where(VFC_DFC_corr==1)


gene_loc_info=pd.read_table('mart_export.txt',header=0)


gene1=[VFC_diff['gene_symbol'].values[x] for x in row]
gene2=[VFC_TF['gene_symbol'].values[x] for x in col]
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

