#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
""" 
date: 2018/02/09 10:51:21 下午 CST
@author :kkwang
"""

import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import os
import numpy as np 
import seaborn as sns
from scipy.spatial.distance import cdist,pdist,squareform
filepath='/Users/kkwang/mywork/gene_array_matrix_csv'
class brainspanwork(object): 
    
    def __init__(self,filepath=filepath):
        self.filepath=filepath
    
    def import_colums_metadata(self):
        
        #print('{0}/columns_metadata.csv'.format(filepath))
        columns_metadata=pd.read_csv('{0}/columns_metadata.csv'.format(self.filepath),index_col=0)
        return columns_metadata
    
      
    def import_row_metadata(self):

        row_metadata=pd.read_csv('{0}/rows_metadata.csv'.format(filepath),index_col=0)
        return row_metadata
    
    def import_expression_matrix(self):
        expression_matrix_col=self.import_colums_metadata().index.tolist()
        expression_matrix=pd.read_csv('{0}/expression_matrix.csv'.format(filepath),index_col=0,names=expression_matrix_col)
        gene_symbol=self.import_row_metadata().loc[:,'gene_symbol'].tolist()
        expression_matrix['gene_symbol']=gene_symbol
        expression_matrix=expression_matrix.drop_duplicates()
        return expression_matrix


class sample_charactor(brainspanwork):
    def __init__(self):
        super().__init__()
    
    def get_tissue_names(self):
        """get sample tissue
           return a list with tissue name in str
           len=all columns
        """
        filee=brainspanwork.import_colums_metadata(self)
        tissue_names=filee.loc[:,'structure_acronym'].tolist()
        return tissue_names
    
    def find_sample_num(self):
        """get the sample number of each tissue return a dict"""
        tissue_names=self.get_tissue_names()
        
        brain_tissue=set(tissue_names)
        brain_tissue_number={}
        for tissue in brain_tissue:
            number=tissue_names.count(tissue)
            brain_tissue_number[tissue]=number
        sorted_brain_tissue_number=sorted(brain_tissue_number.items(),key=lambda x: x[1],reverse=True)
        sorted_brain_tissue_number=dict(sorted_brain_tissue_number)
        return sorted_brain_tissue_number 
    def filter_rare_tissue(self,sample_number,tissue_names,expression_matrix):
        """filter tissue that have less then 5 sample 
           return a matrix  
        """
        tissue_sample_drop=[k for k,v in sample_number.items() if v <=5]
        extract_list=[i+1 for i,t in enumerate(tissue_names) if t in tissue_sample_drop]
        expression_matrix.drop(columns=extract_list,inplace=True)
        return expression_matrix
    

    def find_columns_tissue(self,sample_number,tissue_names):
        """find which columns is which tissue 

        """
        fields={key:[] for key in sample_number.keys()}

        for keys,values in sample_number.items(): 
            for i,p in enumerate(tissue_names):
                if keys == p :
                    fields[keys].append(i+1) # i +1 enumerate start with 0 position start with 1 
            try :
                values == len(fields[keys])
            except:
                print('something is wrong')
        return fields

    def seperate_by_sex(self):
        column_data=brainspanwork.import_colums_metadata(self)
        sex_dic={'female':[],'male':[]}
        for i,p in enumerate(column_data.loc[:,'gender']):
            if p == 'M':
                sex_dic['male'].append(i+1)
            else:
                sex_dic['female'].append(i+1)
        return sex_dic
    def seperate_by_time(self):
        column_data=brainspanwork.import_colums_metadata(self)
        age_dic={'prenatal':['8 pcw','9 pcw','12 pcw','13 pcw','16 pcw','17 pcw','19 pcw','21 pcw','24 pcw','25 pcw','26 pcw'],'early child':['4 mos','10 mos','1 yrs','2 yrs','3 yrs','4 yrs','8 yrs'],'puberty':['13 yrs','15 yrs','18 yrs'],'adult':['21 yrs','23 yrs','30 yrs','36 yrs','37 yrs','40 yrs']}
        time_dic={'prenatal':[],'early child':[],'puberty':[],'adult':[]}

        column_data_age=column_data.loc[:,'age'].tolist()
        
        for i,age in enumerate(column_data_age):
            for time in time_dic.keys():
                if age in age_dic[time]:
                    time_dic[time].append(i+1)
        return time_dic
    def add_pearsonal_info(self,matrix):
        column_data=brainspanwork.import_colums_metadata(self)
        ps_info=list(zip(column_data.loc[matrix.iloc[:,:-1].columns.tolist(),'age'].tolist(),
                             column_data.loc[matrix.iloc[:,:-1].columns.tolist(),'gender'].tolist()))
        ps_info.append('p_info')
        
        return ps_info
        
class expression_deal_with(sample_charactor):
    def __init__(self):
        super().__init__()
    
    def tissue_matrix_workon(self,columns_tissue,tissue):
        """seperate the whole matrix into tissue matrix 
        """
        col_name=columns_tissue[tissue]
        expression_matrix=brainspanwork.import_expression_matrix(self)
        tissue_matrix=expression_matrix.loc[:,col_name]
        tissue_matrix['gene_symbol']=expression_matrix['gene_symbol'].tolist()
        return tissue_matrix
    def boxplot_of_matrix(self,tissue_matrix,name):
        """box plot of the data/tissue matrix 
    
        """
        fig=plt.figure()
        bp=tissue_matrix.iloc[:,:-1].boxplot(sym='r*',patch_artist=True,meanline=True,showfliers=True,return_type='dict')
        for box in bp['boxes']:
            box.set(color='#7570b3',linewidth=1)
            box.set(facecolor='#1b9e77')
        for whisker in bp['whiskers']:
            whisker.set(color='#7570b3',linestyle='--')
        for caps in bp['caps']:
            caps.set(color='g',linewidth=3)
        for median in bp['medians']:
            median.set(color="DarkBlue",linewidth=3)
        for flier in  bp['fliers']:
            flier.set(marker='o', color='#e7298a',alpha=0.5)
        plt.grid(False)
        plt.tight_layout(2,1)
        plt.xticks(rotation=45)
        plt.savefig("{0}.boxplot.pdf".format(name),figsize=(12,7))
    def filter_matrix_workon(self,tissue_matrix):
        """if a gene detected in 80% of the sample
            keep it 
        """
        df_copy=tissue_matrix.copy()
        df_copy_bool=df_copy.astype(bool)
        df_copy=df_copy_bool.astype(int)
        df_copy['sum']=df_copy.sum(axis=1)
        df_copy=df_copy[df_copy['sum']>0.8*len(df_copy.columns)]
        df_index=df_copy.index.tolist()
        df_want=tissue_matrix.loc[df_index,:]
        print(df_want.shape)
        return df_want
    def upperquantile_normalization(self,tissue_matrix):
        """ upperquantile normalization of the matrix 
        """
        tissue_matrix_c=tissue_matrix.copy().iloc[:,:-1]
        upperquantile=[]
        for col in tissue_matrix_c.columns:
            Q3=tissue_matrix_c.loc[:,col].quantile(0.75)
            #print('Q3:{0}'.format(Q3))
            upperquantile.append(Q3)
            tissue_matrix_c.loc[:,col]=tissue_matrix.loc[:,col]/Q3
        median_upperquantile=np.median(upperquantile)
        tissue_matrix_c=tissue_matrix_c * median_upperquantile
        return tissue_matrix_c
    def log2_transform(self,filted_tissue_matrix):
        filted_tissue_matrix_log2=np.log2(filted_tissue_matrix.iloc[:,:-1])
        filted_tissue_matrix_log2['gene_symbol']=filted_tissue_matrix['gene_symbol']
        return filted_tissue_matrix
    def matrix_corr(self,log2_tissue_matrix):
        df_T=log2_tissue_matrix.iloc[:,:-1].T
        df_corr=df_T.corr()
        df_corr['gene_symbol']=log2_tissue_matrix['gene_symbol']
        df_corr_T=df_corr.T
        df_corr_T['gene_symbol']=log2_tissue_matrix['gene_symbol']
        return df_corr_T
class basic_analysis_plot(sample_charactor):
    def __init__(self):
        super().__init__()
    
    def KDE_plot(self,log2_tissue_matrix,tissue,columns_tissue):
        """Kde plot of the matrix,expression distribution
        """
        fig=plt.figure()
        if tissue=='all':
            for keys in columns_tissue.keys():
                #print(columns_tissue[keys])
                x=log2_tissue_matrix.loc[:,columns_tissue[keys]].as_matrix().ravel()
                sns.kdeplot(x,label=keys)
        else:
            x=log2_tissue_matrix.iloc[:,:-1].as_matrix().ravel()
            sns.kdeplot(x)
        plt.legend(fontsize='xx-small')
        plt.xlabel('log2 expression')
        plt.ylabel('Distribution')
        plt.savefig('kde.distribution.pdf')

        
    def check_expression(self,log2_tissue_matrix):
        expression_filter=log2_tissue_matrix[log2_tissue_matrix>6].sum(axis=1).notna()
        exp_matrix=log2_tissue_matrix[expression_filter]
        print(exp_matrix.shape)
        return exp_matrix
    def euclidean_distance(self,log2_tissue_matrix):
        
        dist=pdist(log2_tissue_matrix.iloc[:,:-1],metric='euclidean')
        euclidean_distance_matrix=pd.DataFrame(squareform(dist),columns=log2_tissue_matrix.index,index=log2_tissue_matrix.index)
        print(euclidean_distance_matrix.shape)
        return euclidean_distance_matrix
    
    def sch_clustering(self,log2_tissue_matrix):
        plt.figure('sch_plot')
        dist=pdist(log2_tissue_matrix.iloc[:,:-1].T,metric='euclidean')
        import scipy.cluster.hierarchy as sch
        Z=sch.linkage(dist,method='average') 
        #将层级聚类结果以树状图表示出来并保存为plot_dendrogram.png
        from scipy.cluster.hierarchy import fcluster
        #clusters = fcluster(Z, t=120, criterion='distance')
        #print(clusters)

        P=sch.dendrogram(Z,leaf_rotation=90,leaf_font_size=1,truncate_mode='lastp')
        plt.savefig('plot_dendrogram.pdf')

def main():
    work=brainspanwork() # shili
    expression_matrix=work.import_expression_matrix()

    sample=sample_charactor() # shili
    tissue_names=sample.get_tissue_names() # return a list of full columns lens which tissue it is 
    sample_number=sample.find_sample_num() # return a sorted tissue and number of it type:dict
    #sexual_data=sample.seperate_by_sex() # a dict key: male female
    #time_data=sample.seperate_by_time()
    
    deal_with=expression_deal_with() # shili
    
    work_on_tissue=True #<-- True: work on tissue matrix /False: work on whole matrix 
    tissue='all'
    if work_on_tissue is True:
        filted_expression_matrix=sample.filter_rare_tissue(sample_number,tissue_names,expression_matrix) 
        #filter the tissue that less then five 
        sample_number={key:value for key,value in sample_number.items() if value > 5}
        columns_tissue=sample.find_columns_tissue(sample_number,tissue_names) #find which cloumns is which tissue (

        tissue='CBC'  #<--change you tissue there
        workon_matrix=deal_with.tissue_matrix_workon(columns_tissue,tissue)
        deal_with.boxplot_of_matrix(workon_matrix,tissue) #boxplot of workon_matrix
    else:
        workon_matrix=expression_matrix
    
    workon_matrix_nor=deal_with.upperquantile_normalization(workon_matrix)
    deal_with.boxplot_of_matrix(workon_matrix_nor,"{0}_upperquantile_normalize".format(tissue))
    filted_workon_matrix=deal_with.filter_matrix_workon(workon_matrix)
    log2_workon_matrix=deal_with.log2_transform(filted_workon_matrix)
    
    basic_analysis=basic_analysis_plot() #shili
    log2_workon_matrix=basic_analysis.check_expression(log2_workon_matrix)
    ps_info=sample.add_pearsonal_info(log2_workon_matrix)
    print(ps_info)
    log2_workon_matrix.to_csv('{0}_log2_fiter_matrix.tsv'.format(tissue),sep='\t',index=False)
    euclidean_distance_matrix=basic_analysis.euclidean_distance(log2_workon_matrix)
    basic_analysis.sch_clustering(log2_workon_matrix)

    #basic_analysis.KDE_plot(log2_workon_matrix,tissue,columns_tissue)

    #cor_workon_matrix=deal_with.matrix_corr(log2_workon_matrix)
    #print(cor_workon_matrix.shape)

    




if __name__ == "__main__":
    main()
