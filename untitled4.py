#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 23:09:04 2018

@author: kkwang
"""
import pandas as pd
mt =['MD', 'CBC','AMY', 'HIP','STR','VFC', 'A1C','S1C','IPC','M1C', 'DFC','MFC', 'OFC' ,'ITC','V1C','STC']
tissue_genes={}     

for mm in mt:  # put all genes into gene_array
    tissue_genes[mm]=pd.read_csv('/Users/kkwang/mywork/gene_array_matrix_csv/{}_limma.tsv'.format(mm),index_col=0,header=0,sep='\t').index.values    

from itertools import combinations
choice=combinations(tissue_genes.keys(),2)

with open('/Users/kkwang/mywork/gene_array_matrix_csv/common_gene_file.txt' ,'r')as f:
    f.readlines()
    common_genes=[x[:-1] for x in f]
tissue_shared_genes={}# a dict of gene shared by tissue 
for com in choice:
    tissue_shared_genes[com]=set(tissue_genes[com[0]]).intersection(set(tissue_genes[com[1]]))-set(common_genes)


with open('VFC_DFC.txt','r') as f:
    VFC_DFC=f.readlines()
    VFC_DFC=[gene[:-1] for gene in VFC_DFC]
    
gene_loc_info=pd.read_table('mart_export.txt',header=0)
gene_loc_info.columns=['id','chr','start','end','gene']

VFC_DFC_loc=gene_loc_info[gene_loc_info['gene'].isin(VFC_DFC)]


from collections import Counter

hhhh=Counter(VFC_DFC_loc['chr'])
X=hhhh['X']
del hhhh['X']
sk=sorted(map(int,hhhh.keys()))
val=[hhhh[str(k)] for k in sk]

# library
import matplotlib.pyplot as plt
 
# create data
label=['chr'+str(k) for k in sk]
label.append('chrX')
val.append(X)

explode = [0.1]* len(label)
 
fig, ax = plt.subplots(figsize=(10,8))
ax.axis('equal')
mypie, _ = ax.pie(val, radius=1.3 ,labels=label,explode=explode,labeldistance=0.8)
plt.setp( mypie, width=0.3, edgecolor='white')
plt.show()
 

