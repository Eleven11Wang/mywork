#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
""" 
date: 2018/03/17 01:04:06 下午 CST
@author :kkwang

#what make tissue differernt:

for each tissue :
    relation of high corr gene to chromosol 
    high corr gene common TF 


#what make some tissue same or share some common function or how they interact 

for each clustering of the sample 
find the commen gene of each tissue (remove the the shared by all tissue)

dict:'A+B':[gene_list]
for gene in gene_list:
    


""" 

import pandas as pd

VFC_data=pd.read_csv('VFC_log2_fiter_matrix.tsv',sep='\t',header=0)

VFC_gene_symbol=VFC_data['gene_symbol']
VFC_data_T=VFC_data.iloc[:,:-1].T
VFC_corr=VFC_data_T.corr()
del VFC_data_T

VFC_corr=VFC_corr.abs()
VFC_corr[VFC_corr<0.8]=0
VFC_bool=VFC_corr.astype(bool)
del VFC_corr


f=open('VFC_highcorr_gene.tsv','w')
f.write('gene1\tgene2')

for row,gene1 in enumerate(VFC_gene_symbol):
    for col,gene2 in enumerate(VFC_gene_symbol):
        if VFC_bool.iloc[row,col]:
            f.write('{0}\t{1}'.format(gene1,gene2))

f.close()

