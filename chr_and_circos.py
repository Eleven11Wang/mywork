#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 18:19:12 2018

@author: kkwang
"""

# match chromosol and gene and circular plot 
#
#
import pandas as pd 
import numpy as np
gene_loc_info=pd.read_table('mart_export.txt',header=0)
gene_loc_info.columns=['id','chr','start','end','gene']
#2241 genes

tissue='STR'
tissue_DE_gene=pd.read_csv('/Users/kkwang/mywork/gene_array_matrix_csv/{}_limma.tsv'.format(tissue),
            index_col=0,header=0,sep='\t').index.values
tissue_data=pd.read_csv('/Users/kkwang/mywork/gene_array_matrix_csv/{}_log2_fiter_matrix.tsv'.format(tissue),sep='\t')                    
DE_data=tissue_data[tissue_data['gene_symbol'].isin(tissue_DE_gene)]
DE_data_gene=DE_data['gene_symbol'].values
##get DE gene exp data


DE_data_T=DE_data.iloc[:,:-1].T
DE_corr=DE_data_T.corr()
DE_corr.index=DE_data_gene
DE_corr.columns=DE_data_gene
##Corr_regulatory network


DE_corr[DE_corr==1]=0.999
DE_corr[DE_corr==-1]=-0.999
# fisher z transform z’ = 0.5[ln(1+r) – ln(1-r)]
DE_corr_z=0.5 *(np.log(DE_corr+1)-np.log(1-DE_corr))
DC_corr_h=DE_corr_z[DE_corr_z.abs()<2.7]=0
rows,cols=np.where(DE_corr_z>2)
gene_pairs =pd.DataFrame(columns=['gene1','gene2'])
gene_pairs['gene1']=DE_corr.index[rows]
gene_pairs['gene2']=DE_corr.columns[cols]




##filter high corr genes
gene_loc_info=pd.read_table('mart_export.txt',header=0)
gene_loc_info.columns=['id','chr','start','end','gene']

chr_df=pd.DataFrame(index=gene_loc_info['gene'],columns=['chr'])
chr_df['chr']=gene_loc_info['chr'].values
chr_df=chr_df[chr_df['chr'].str.len()<3]
chr_dict=chr_df['chr'].to_dict()

#map gene to chr
def chr_map(gene):
    try:
        return chr_dict[gene]
    except:
        return np.NaN
chr_c=gene_pairs.applymap(chr_map).dropna()
chr_i=zip(chr_c['gene1'],chr_c['gene2'])

# count chr
from collections import Counter
coun=Counter(chr_i)
genes=coun.keys()
gene1,gene2=zip(*list(genes))
num=[coun[g] for g in genes]


# output of data
hhhhh=pd.DataFrame(index=np.arange(len(coun)),columns=['gene1','gene2','num'])
hhhhh['gene1']=gene1
hhhhh['gene2']=gene2
hhhhh['num']=num
#hhhhh=pd.DataFrame(index=list(set([x[0] for x in coun.keys()])),columns=list(set([x[0] for x in coun.keys()])))
#for keys,it in coun.items():
#    hhhhh.loc[keys[0],keys[1]]=it
#hhhhh.fillna(0)

hhhhh.to_csv('chr_con_{}.tsv'.format(tissue),sep='\t')




