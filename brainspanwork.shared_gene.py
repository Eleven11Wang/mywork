#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 21:40:58 2018

@author: kkwang
"""

import pandas as pd
TF_gene=pd.read_csv('trrust_rawdata.human.tsv',sep='\t',names=['TF','gene','xx','xxx'])

differ_data_VFC=pd.read_csv('VFC_limma.tsv',header=0,sep='\t')
differ_gene_VFC=differ_data_VFC.index.tolist()

differ_data_DFC=pd.read_csv('DFC_limma.tsv',header=0,sep='\t')
differ_gene_DFC=differ_data_DFC.index.tolist()

del differ_data_VFC
del differ_data_DFC


with open('common_gene_file.txt', 'r') as f:
    myNames = f.readlines()
common_gene=[x[:-1] for x in myNames]

VFC_DFC=set(differ_gene_DFC).intersection(set(differ_gene_VFC))-set(common_gene)


with open('VFC_DFC.txt','w') as f:
    f.write('\n'.join(VFC_DFC))
    
    
    
TF_VFC_DFC=TF_gene[TF_gene['gene'].isin(VFC_DFC)]
