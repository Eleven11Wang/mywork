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

        

def main():
    work=brainspanwork()
    expression_matrix=work.import_expression_matrix()
    
    sample=sample_charactor()
    tissue_names=sample.get_tissue_names() # return a list of full columns lens which tissue it is 
    sample_number=sample.find_sample_num() # return a sorted tissue and number of it type:dict

    filter_by_five=False #<-- change whether you what to filter rare tissue 

    if filter_by_five is True:
        filted_expression_matrix=sample.filter_rare_tissue(sample_number,tissue_names,expression_matrix) 
        #filter the tissue that less then five 
        sample_number={key:value for key,value in sample_number.items() if value > 5}
    
    columns_tissue=sample.find_columns_tissue(sample_number,tissue_names) #find which cloumns is which tissue (
    sexual_data=sample.seperate_by_sex()
   
    print(sexual_data)
    time_data=sample.seperate_by_time()
    print(time_data)
if __name__ == "__main__":
    main()
