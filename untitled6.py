#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 20:13:52 2018

@author: kkwang
"""

import pandas as pd
gene_12_chart=pd.read_table('common_12_chart.txt',sep='\t',header=0)



choice_chart=gene_12_chart[gene_12_chart['Fold Enrichment']>2]
choice_chart=choice_chart[choice_chart['Benjamini']<0.05]
choice_chart=choice_chart[choice_chart['Count']>20]

choice_chart['Term']=['Calcium', 'Ion transport', 'Cell junction',
       'calcium ion binding', 'Cell cerface junction',
       'cell surface', 'Ion channel', 'Synapse',
       'nervous system development',
       'chemical synaptic transmission',
       'Chromosomal rearrangement', 'neuronal cell body']


bar_plot_char=choice_chart.loc[:,['Term','Count','Fold Enrichment']]
bar_plot_char=bar_plot_char.sort_values('Count',ascending=False)
import numpy as np
import matplotlib.pyplot as plt
Term=bar_plot_char['Term'].values
er=bar_plot_char['Fold Enrichment'].values
er=['{0:.2f}'.format(x) for x in er]
Term=zip(Term,er)
x=np.arange(12)
plt.figure(figsize=(9,12))
Y = bar_plot_char['Count'].values
K=zip(x,Y)
plt.barh(x,Y,facecolor='#ff9999', edgecolor='white')

er=bar_plot_char['Fold Enrichment'].values
for K,term in zip(K,Term):
    plt.text(1,K[0]-0.2, term, ha='left', va= 'bottom')

plt.xlabel('Numbers of genes')
plt.yticks([])
plt.savefig('hhhh.png')

tissue='AMY'
VFC_data=pd.read_csv('gene_array_matrix_csv/{}_log2_fiter_matrix.tsv'.format(tissue),sep='\t',header=0)
choice_gene=choice_chart['Genes']

choice_term=choice_chart['Term'].values



fig, axes = plt.subplots(nrows=4, ncols=3,figsize=(18,12))

col_info=pd.read_csv('gene_array_matrix_csv/columns_metadata.csv',header=0)

age=col_info[col_info['column_num'].isin(map(int,VFC_data.columns.values[:-1]))]['age']
for inx,term in enumerate(choice_term):
    row=int(inx/3)
    col=inx%3
    choice_genes=choice_gene.iloc[inx].lstrip().split(',')
    genes=[x[1:] for x in choice_genes[1:]]
    genes[0:0]=choice_genes[0]


    VFC_C=VFC_data[VFC_data['gene_symbol'].isin(genes)].iloc[:,:-1]
    VFC_C=VFC_C.diff(periods=1, axis=1)
    plt.tight_layout() 
    
    
    ax=VFC_C.T.plot(legend=False,title=term,ax=axes[row,col])
    ax.xaxis.set_major_locator(plt.MaxNLocator(len(VFC_C.columns)))
    #ax.set_xticks(VFC_C.index)
    ax.set_xticklabels(list(age), rotation=45)
    
plt.savefig('coomon_gene_expression_pattern_{}.pdf'.format(tissue),transparent=True)
   
fig=plt.figure(figsize=(12,9)  )
inx=11
choice_genes=choice_gene.iloc[inx].lstrip().split(',')
genes=[x[1:] for x in choice_genes[1:]]
genes[0:0]=choice_genes[0]
VFC_C=VFC_data[VFC_data['gene_symbol'].isin(genes)].iloc[:,:-1]
VFC_C=VFC_C.diff(periods=1, axis=1)
_x = VFC_C['gene_symbol'].values

VFC_C.columns=age

plt.tight_layout()   
ax=VFC_C.T.plot(legend=False,title=choice_term[inx])
ax.xaxis.set_major_locator(plt.MaxNLocator(len(VFC_C.columns)))
    #ax.set_xticks(VFC_C.index)
ax.set_xticklabels(list(age), rotation=45)
plt.ylabel('Difference of two stage')
plt.savefig('AMY_neuron.pdf')
