#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 23:01:47 2018

@author: kkwang
"""

import pandas as pd 

STR=pd.read_csv('/Users/kkwang/mywork/chr_con_STR.tsv',sep='\t',index_col=0)
HIP=pd.read_csv('/Users/kkwang/mywork/chr_con_HIP.tsv',sep='\t',index_col=0)
AMY=pd.read_csv('/Users/kkwang/mywork/chr_con_AMY.tsv',sep='\t',index_col=0)
tis=[STR,HIP,AMY]
chrs=list(set(list(STR['gene1'])))

colnames=['STR','HIP','AMY']

c_df_s=pd.DataFrame(index=chrs,columns=colnames)
c_df_d=pd.DataFrame(index=chrs,columns=colnames)


chr_gene_num={'chr1':3000,'chr2':2500,'chr3':1900,'chr4':1600,
              'chr5':1700,'chr6':1900,'chr7':1800,'chr8':1400,
              'chr9':1400,'chr10':1400,'chr11':2000,'chr12':1600,
              'chr13':800,'chr14':1200,'chr15':1200,'chr16':1300,
              'chr17':1600,'chr18':600,'chr19':1700,'chr20':900,
              'chr21':400,'chr22':800,'chrX':1400}

for i,tissue in enumerate(tis):
    
    
    c_dict_s = {ch: 0 for ch in chrs}
    c_dict_d=  {ch: 0 for ch in chrs}
    
    for tuo in tissue.itertuples():
        if tuo[1]==tuo[2]:
            
            c_dict_s[tuo[1]]=c_dict_s[tuo[1]]+float(tuo[3]/2)
            c_dict_s[tuo[2]]+=float(tuo[3]/2)
        else:
            c_dict_d[tuo[1]]=c_dict_d[tuo[1]]+float(tuo[3]/2)
            c_dict_d[tuo[2]]+=float(tuo[3]/2)
    
    for k,v in c_dict_s.items():
        c_df_s.loc[k,colnames[i]]=v
    for k,v in c_dict_d.items():
        c_df_d.loc[k,colnames[i]]=v


c_df=c_df_s+c_df_d    
c_df.index=['chr'+str(x) for x in c_df.index.values]

c_df['chr_gene_num']=pd.Series(c_df.index).map(chr_gene_num).values


gene_num=c_df.iloc[:,:3].sum()
gene_p=gene_num/c_df['chr_gene_num'].sum()



c_df['STR']=c_df['STR']/c_df['chr_gene_num']
c_df['HIP']=c_df['HIP']/c_df['chr_gene_num']
c_df['AMY']=c_df['AMY']/c_df['chr_gene_num']

from scipy import stats
for ch in c_df.index.values:
    
    st,p_vals=stats.ttest_ind(c_df.loc[ch,colnames],gene_p)
    c_df.loc[ch,'t-st']=st
    c_df.loc[ch,'p-val']=p_vals

sig_df=c_df[c_df['p-val']<0.05]
sign_up=sig_df[sig_df['t-st']>0].index.values
sign_down=sig_df[sig_df['t-st']<0].index.values


c_df_c=c_df.copy()
c_df=c_df.iloc[:-1,:]
c_df['chr']=[x[3:] for x in c_df.index.values]
c_df['chr']=list(map(int,c_df['chr']))
c_df.sort_values(by='chr',inplace=True)
c_df.append(c_df_c.loc['chrX',:])

import numpy as np
import matplotlib.pyplot as plt
fig=plt.figure(figsize=(12,9))
barWidth = 0.2
bars1 = c_df['STR'].tolist()
bars2 = c_df['HIP'].tolist() 
bars3=c_df['AMY'].tolist()
 
# The x position of bars
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
c1=[]
c2=[]
c3=[]
for ch in c_df.index:
    if ch in sign_up:
        c1.append('red')
        c2.append('orange')
        c3.append('peru')
    elif ch in sign_down:
        c1.append('green')
        c2.append('cyan')
        c3.append('deepskyblue')
    else:
        c1.append('grey')
        c2.append('grey')
        c3.append('grey')


bar1=plt.bar(r1, bars1, width = barWidth, color = c1,   capsize=7, label='STR')

bar2=plt.bar(r2, bars2, width = barWidth, color = c2,   label='HIP')
bar3=plt.bar(r3, bars3, width = barWidth, color = c3, label='AMY')


# general layout
plt.xticks([r + barWidth for r in range(len(bars1))], c_df.index.tolist(),rotation=45)
plt.ylabel('Percent of chromosomal genes')
plt.xlabel('Significant chromosomes p-val<0.05')
plt.title('Limbic Area')
import matplotlib.patches as mpatches
green_patch = mpatches.Patch(color='green', label='STR_down')
cyan_patch = mpatches.Patch(color='cyan', label='HIP_down')
sky_patch = mpatches.Patch(color='deepskyblue', label='AMY_down')
red_patch = mpatches.Patch(color='red', label='STR_up')
orange_patch = mpatches.Patch(color='orange', label='HIP_up')
peru_patch = mpatches.Patch(color='peru', label='AMY_up')


plt.legend(handles=[green_patch,red_patch,cyan_patch,orange_patch,sky_patch,peru_patch],fontsize='x-large')
plt.savefig('limbic_chr.pdf')
# Show graphic
plt.show()
"""
c_df_dc=c_df_d.copy()
c_df_d=c_df_d.iloc[:-1,:]
c_df_d['chr']=c_df_d.index.values

#c_df['chr']=list(map(int,c_df['chr']))
c_df.sort_values(by='chr',inplace=True)
c_df.append(c_df_dc.loc['chrX',:])

#c_df.iloc[:,:-1].plot.bar()
c_df_d['sum']=c_df_d.sum(axis=1)

labels =['chr'+str(x) for x in c_df_d.index]
sizes = c_df_d['sum'].values

#explode = (0.1, 0, 0, 0)  # explode 1st slice
#explode=explode,
# Plot

explode = [0.1]* len(c_df_d)
 
fig, ax = plt.subplots(figsize=(10,10))
ax.axis('equal')

plt.pie(sizes,explode=explode,  labels=labels, shadow=True)
ax.axis('equal') 
plt.axis('equal')
plt.show()
"""