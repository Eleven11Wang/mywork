#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 21:54:38 2018

@author: kkwang
"""

import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
from mpl_toolkits.mplot3d import Axes3D




## filter the enrichment gene
gene_12_chart=pd.read_table('common_12_chart.txt',sep='\t',header=0)
choice_chart=gene_12_chart[gene_12_chart['Fold Enrichment']>2]
choice_chart=choice_chart[choice_chart['Benjamini']<0.05]
choice_chart=choice_chart[choice_chart['Count']>20]

## get of specific pathway 
gene_want=choice_chart.loc[37,'Genes']
gene_want=gene_want.replace(' ','').split(',')
gene_want_10=gene_want[:10]


## imput specific tissue_data
tissue='AMY'
tissue_data=pd.read_csv('gene_array_matrix_csv/{}_log2_fiter_matrix.tsv'.format(tissue),sep='\t',header=0)

tissue_genes=tissue_data[tissue_data['gene_symbol'].isin(gene_want)].iloc[:,:-1]
tissue_diff=tissue_genes.diff(periods=1, axis=1).abs()
tissue_diff[tissue_diff<2]=0.01
tissue_diff.fillna(0.2,inplace=True)
tissue_diff_10=tissue_diff.iloc[:10,:]
##impot age_data
col_info=pd.read_csv('gene_array_matrix_csv/columns_metadata.csv',header=0)
age=col_info[col_info['column_num'].isin(map(int,tissue_data.columns.values[:-1]))]['age']

fig=plt.figure(figsize=(24,10),dpi=300,facecolor=None)
fig.patch.set_alpha(0.5)
ax1 = fig.add_subplot(121, projection='3d')

xlabels = age
ylabels = gene_want_10
_x = np.arange(len(age))
_y = np.arange(len(gene_want_10))
_xx, _yy = np.meshgrid(_x, _y)
x, y = _xx.ravel(), _yy.ravel()



top=tissue_diff_10.values.flatten()
bottom = np.zeros_like(top)
width = depth = 0.5

#ax1.view_init(30, 120)
ax1.w_xaxis.set_ticks(_x + width/2.)
ax1.w_xaxis.set_ticklabels(xlabels,rotation=90)

ax1.w_yaxis.set_ticks(_y + depth/2.)
ax1.w_yaxis.set_ticklabels(ylabels,rotation=-45)
ax1.set_zlim3d(0, top.max())
ax1.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax1.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
values = np.linspace(0.2, 1., x.ravel().shape[0])
colors = plt.cm.rainbow(values)
ax1.grid(False)
ax1.set_zticks([])
ax1.bar3d(x, y, bottom, width, depth, top, shade=False,color=colors)
#ax1.set_title('Ion_channel_enrichmen_gene')
plt.savefig('AMY_neuron_cell_body_10_r.png',transparent=True)