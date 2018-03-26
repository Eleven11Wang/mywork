#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 13:14:45 2018

@author: kkwang
"""

import pandas as pd
import numpy as np 


#common_gene_enrich_in_network
tissue='AMY'
tissue_data=pd.read_csv('gene_array_matrix_csv/{}_log2_fiter_matrix.tsv'.format(tissue),sep='\t',header=0)
tissue_gene=tissue_data['gene_symbol']

tissue_T=tissue_data.iloc[:,:-1].T
tissue_corr=tissue_T.corr()
tissue_corr[tissue_corr==1]=0.1
tissue_corr[tissue_corr.abs()>0.8]=2
rows,cols=np.where(tissue_corr==2)

import networkx as nx 
edges = zip(tissue_gene[rows.tolist()],tissue_gene[cols.tolist()])

gr = nx.Graph()
gr.add_edges_from(edges)

with open('common_gene_file_12.txt','r') as f:
    gene_12=f.readlines()
    gene_12=[g[:-1] for g in gene_12]
    
node,degree=zip(*gr.degree())
#network_info=pd.DataFrame(index=node,columns=['degree','betwenness','centrality'])
network_info=pd.DataFrame(index=node,columns=['degree'])

#node,betwenness=zip(*nx.betweenness_centrality(gr))
#node,centrality=zip(*nx.closeness_centrality(gr))

network_info['degree']=degree
#network_info['betw']
network_info.sort_values(by='degree',inplace=True)
degree_s=network_info['degree'].values

degree_per_n=[]
per=10
while per>0  :
    cut_num=degree_s[(per-1)*int(len(network_info)/10):per*int(len(network_info)/10)]
    degree_per_n.append((cut_num.max(),cut_num.min()))
    per-=1

per_list=[] 
x_formats=[]
median_info=[]
def get_median(data):
    data.sort()
    half = len(data) // 2
    return (data[half] + data[~half]) / 2

for x ,deg in enumerate(degree_per_n):
    on_cut_info=network_info[ (network_info['degree']<deg[0]) & (network_info['degree']>=deg[1])]
    common_in_info=on_cut_info[on_cut_info.index.isin(gene_12)]
    common_in_gene=len(common_in_info)    
    total=len(on_cut_info)
    median_t=get_median(on_cut_info['degree'].values)
    median_c=get_median(common_in_info['degree'].values)
    median_p=median_c/median_t
    median_info.append(median_p)
    x_format=total/len(network_info)
    x_formats.append(x_format)
    per=common_in_gene/total
    per_list.append(per)
#between=nx.betweenness_centrality(gr)

x_formats=[sum(x_formats[:n]) for n in np.arange(len(x_formats))]
x_formats=[round(x,2) for x in x_formats]
x_formats.append(1)

import numpy as np
import matplotlib.pyplot as plt
plt.figure(figsize=(12,9),dpi=200)
barWidth = 0.5
height = per_list
bars = ['{}%<=X<={}%'.format(x_formats[i]*100,x_formats[i+1]*100) for i in np.arange(len(per_list))]

y_pos = np.arange(len(bars))

color=['grey']*6
color[0:0]=['r']*2
color[-1:]=['g']*3

plt.bar(y_pos, height,width=barWidth,color=color)
plt.xticks(y_pos, bars,rotation=90)
plt.ylabel('propertion of tissue shared gene')
plt.xlabel('Degree catalogy(top%)')
plt.show()
plt.savefig('common_gene_Degree_propertion_AMY.png')


tissue_corr=tissue_T.corr()
tissue_corr.index=tissue_gene
tissue_corr.columns=tissue_gene
corr_12=tissue_corr[tissue_corr.index.isin(gene_12)]
corr_12_T=corr_12.T
corr_12_12=corr_12_T[corr_12_T.index.isin(gene_12)]

hi
# set width of bar

#r2 = [x + barWidth for x in r1]
#r3 = [x + barWidth for x in r2]

#plt.bar(r1, bars1, color='#7f6d5f', width=barWidth, edgecolor='white', label='degree')
#plt.bar(r2, bars2, color='#557f2d', width=barWidth, edgecolor='white', label='betweenness')
#plt.bar(r3, bars3, color='#2d7f5e', width=barWidth, edgecolor='white', label='centrality')

#plt.xlabel('group', fontweight='bold')
#plt.xticks([r + barWidth for r in range(len(bars1))], ['top{}_{}'.format(x_format[i])])
 
# Create legend & Show graphic
