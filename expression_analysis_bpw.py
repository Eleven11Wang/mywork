#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 14:15:25 2018

@author: kkwang
"""


import pandas as pd 

kk=pd.read_csv('/Users/kkwang/mywork/all_log2_fiter_matrix.tsv',header=0,sep='\t')

period_info=pd.read_csv('/Users/kkwang/mywork/syllabus.2.csv',header=0)  # get the period information 
period_info_dict=period_info['col_numbers'].to_dict()  

disease_file=pd.read_csv('spci_df.tsv',sep='\t',index_col=0,header=0)


tissue_cla=disease_file['tissue'].unique()



common_gens=open('/Users/kkwang/mywork/specific_gene_file.txt','r')
common_gens=[x[:-1] for x in common_gens.readlines()]

kk.common=kk[kk['gene_symbol'].isin(common_gens)]



#tissue_cla=['MD', 'CBC','AMY', 'HIP','STR','VFC', 'A1C','S1C','IPC','M1C', 'DFC','MFC', 'OFC' ,'ITC','V1C','STC']

col_info=pd.read_csv('/Users/kkwang/mywork/gene_array_matrix_csv/columns_metadata.csv',index_col=0,header=0)
col_tissue_info= col_info['structure_acronym'].tolist()


mean_exp_df=pd.DataFrame(index=tissue_cla,columns=['period'+str(x) for x in range(12)])

tissue_dict={tissue:[] for tissue in tissue_cla}
for inx,tis in enumerate(col_tissue_info):
   try:
       tissue_dict[tis].append(inx)
   except:
       continue
     

network_gene=['KIF23',	'CASC5'	,'KIAA0101',	'RRM2',	'NCAPH',	'KIF18A','TOP2A',	'SGOL1',	'ASPM',	'KIF15',	
'TPX2',	'KIF20A'	,'MLF1IP',	
'BUB1B',	'NDC80', 	'BRIP1',	'NEIL3',	'CCNB2',	'DTL',	'TTK',	'KIF14',	'PBK',	
'KIF11',	'DLGAP5','BUB1',	'FANCI',	'CDK1','CKAP2L','ERBB3',	'C11orf9']	



for ab in period_info_dict.items():
    period='period'+str(ab[0])
    col_range=ab[1].split(',')
    min_col=int(col_range[0])
    max_col=int(col_range[1])
    columns_want=[int(col) for col in kk.common.iloc[:,:-1].columns.tolist() if min_col<=int(col)<=max_col]
    for tissue in tissue_cla:
        tissue_col=tissue_dict[tissue]
        both_col=list(set(columns_want).intersection(set(tissue_col)))     
        
        genes_want=disease_file[disease_file['tissue']==tissue]
    
        #cl_want=kk[kk['gene_symbol'].isin(network_gene)]
        cl_want=kk[kk['gene_symbol'].isin(genes_want['gene_symbol'])]
        #cl_want=kk[kk['gene_symbol'].isin(disease_file['gene_symbol'].unique())]
        mean_exp=cl_want.iloc[:,both_col].values.mean()
        mean_exp_df.loc[tissue,period]=mean_exp


mean_exp_df.fillna(mean_exp_df.mean().mean(),inplace=True)

import seaborn as sns; sns.set()
mean_exp_df[mean_exp_df<8]=0
sns.clustermap(mean_exp_df,col_cluster=False,row_cluster=False,cmap='plasma')




