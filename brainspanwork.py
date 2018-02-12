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

filepath='/Users/kkwang/mywork/genes_matrix_csv'
class brainspanwork(object): 
    
    def __init__(self,filepath=filepath):
        self.filepath=filepath
    
    def import_colums_metadata(self):
        
        #print('{0}/columns_metadata.csv'.format(filepath))
        columns_metadata=pd.read_csv('{0}/columns_metadata.csv'.format(self.filepath),index_col=0)
        return columns_metadata
    
    def get_tissue_names(self):
        """get sample tissue
           return a list with tissue name in str
           len=all columns
        """
        filee=self.import_colums_metadata()
        tissue_names=filee.loc[:,'structure_acronym'].tolist()
        return tissue_names
    
    def import_row_metadata(self):

        row_metadata=pd.read_csv('{0}/rows_metadata.csv'.format(filepath),index_col=0)
        return row_metadata
    
    def import_expression_matrix(self):
        expression_matrix=pd.read_csv('{0}/expression_matrix.csv'.format(filepath),index_col=0,header=None)
        expression_matrix=pd.DataFrame(expression_matrix,index=self.import_row_metadata().index.tolist())
        expression_matrix=expression_matrix.drop_duplicates()
        return expression_matrix

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



def main():
    work=brainspanwork()
    expression_matrix=work.import_expression_matrix()
    tissue_names=work.get_tissue_names() # return a list of full columns lens which tissue it is 
    sample_number=work.find_sample_num() # return a sorted tissue and number of it type:dict

    filter_by_five=False #<-- change whether you what to filter rare tissue 

    if filter_by_five is True:
        
        filted_expression_matrix=work.filter_rare_tissue(sample_number,tissue_names,expression_matrix) 
        #filter the tissue that less then five 
        sample_number={key:value for key,value in sample_number.items() if value > 5}
    
    columns_tissue=work.find_columns_tissue(sample_number,tissue_names) #find which cloumns is which tissue (
   
    print(columns_tissue)
    print(expression_matrix.index)
    print(expression_matrix.columns)
if __name__ == "__main__":
    main()
