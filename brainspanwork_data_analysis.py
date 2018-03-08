#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 17:09:54 2018



mt(list)
    tissue name
    
tissue_genes(dict) 
    tissue_genes[tissue]: differential expression gene of this tissue


tissue_data(dict)
    tissue_data[tissue]: specific differtial expression genes  daraframe 
    index: index number  of original df
    columns: sample 
    
period_data(dict)
    period_data[tissue]:period specific differential expression genes df
    index : index number  of original df
    columns : period+ genes +tissue 
    
        
specific_genes,common_genes(list):
    genes that are shared in <8 tissue/>8 tissue

res(df)
    concat all tissue period data together 
    index :number
    columns: period +gene+ tissue
    
    
@author: kkwang
"""


#import the package needed
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations
import networkx as nx
import seaborn as sns; sns.set()
from sklearn.cluster import KMeans
from PIL import Image




mt =['MD', 'CBC','AMY', 'HIP','STR','VFC', 'A1C','S1C','IPC','M1C', 'DFC','MFC', 'OFC' ,'ITC','V1C','STC']

tissue_genes={}     
tissue_data={}
period_data={}
gene_array=np.empty((0))


for mm in mt:  # put all genes into gene_array
    tissue_genes[mm]=pd.read_csv('/Users/kkwang/mywork/{}_limma.tsv'.format(mm),index_col=0,header=0,sep='\t').index.values    
    gene_array=np.append(gene_array,tissue_genes[mm])
 
    
unique, counts = np.unique(gene_array, return_counts=True) # howmany times the gene differtial expression in 16 tissue
gene_count_dict=dict(zip(unique, counts))  
specific_genes=[]
common_genes=[]


for x in gene_count_dict.items():
    if x[1]>=8:
        common_genes.append(x[0])
    else:
        specific_genes.append(x[0])  # find shared genes and specific genes
    
# write common genes to a file 
common_gene_file=open('common_gene_file.txt','w')
common_gene_file.write("\n".join(common_genes))
common_gene_file.close()



G=nx.Graph()
choice=combinations(tissue_genes.keys(),2)

# how many genes shared by each tissue
df=pd.DataFrame(index=mt,columns=mt)
for x in choice:
    xxx=set(list(tissue_genes[x[0]])).intersection(list(tissue_genes[x[1]]))
    G.add_edge(x[0],x[1],weight=len(xxx))
    df.loc[x[0],x[1]]=len(xxx)
    df.loc[x[1],x[0]]=len(xxx)
df.fillna(1000,inplace=True)

# clustering of the share gene_df
km = KMeans(n_clusters=3).fit(df)
km_labels=km.labels_.tolist()




# put the cluster network into a brain image 
img = Image.open('/Users/kkwang/Desktop/1.45.26.png')
plt.figure(figsize=(8,8))
plt.imshow(img,alpha=0.5)
plt.xticks([]), plt.yticks([])  # to hide tick values on X and Y axis

nodes=df.index.values


pos = {'DFC': (120, 200), 'VFC': (120, 230), 'MFC': (120,270),'OFC':(170,330),
       'MD':(370,300),'AMY':(350,330),'HIP':(350,360),'STR':(350,390),
       'STC':(370,270),'ITC':(300,370),'M1C':(350,100),'S1C':(410,100),
       'V1C':(600,300),'CBC':(500,400),'IPC':(500,125),'A1C':(430,240)} 

edges = G.edges()
nx.draw_networkx_nodes(G,pos,nodelist=[ nodes[ind] for ind,x in enumerate(km_labels) if x==0],node_color='pink')
nx.draw_networkx_nodes(G,pos,nodelist=[ nodes[ind] for ind,x in enumerate(km_labels) if x==1],node_color='cyan')
nx.draw_networkx_nodes(G,pos,nodelist=[ nodes[ind] for ind,x in enumerate(km_labels) if x==2],node_color='orange')
nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
nx.draw_networkx_labels(G,pos)                                        


tissue_left=list(set(mt)-set(['STR','AMY','HIP','MD','CBC']))
df=df.loc[tissue_left,:]
km = KMeans(n_clusters=3).fit(df)
km_labels=km.labels_.tolist()


plt.figure(figsize=(8,8))
plt.imshow(img,alpha=0.5)
plt.xticks([]), plt.yticks([]) 
G.remove_nodes_from(['STR','AMY','HIP','MD','CBC'])
nodes=df.index.values
nx.draw_networkx_nodes(G,pos,nodelist=[ nodes[ind] for ind,x in enumerate(km_labels) if x==0],node_color='magenta')
nx.draw_networkx_nodes(G,pos,nodelist=[ nodes[ind] for ind,x in enumerate(km_labels) if x==1],node_color='y')
nx.draw_networkx_nodes(G,pos,nodelist=[ nodes[ind] for ind,x in enumerate(km_labels) if x==2],node_color='teal')
nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
nx.draw_networkx_labels(G,pos)  




# clustermap of the tissue_share_genes_df
tissue_cla=['MD', 'CBC','AMY', 'HIP','STR','VFC', 'A1C','S1C','IPC','M1C', 'DFC','MFC', 'OFC' ,'ITC','V1C','STC']
tissue_cla_pal = sns.color_palette("Set2",5 )
cla_dict={0:2,1:3,2:5,3:4,4:2}
tissue_pal=[] # make the color bar
for index,x in cla_dict.items():
    while x>0:
        tissue_pal.append(tissue_cla_pal[index])
        x-=1
        
tissue_lut=dict(zip(tissue_cla, tissue_pal))
xxx=pd.Series(data=tissue_cla, index=tissue_cla)
row_colors_df = xxx.map(tissue_lut)
share_tissue_heatmap =sns.clustermap(df, cmap="coolwarm",robust=True,row_colors=row_colors_df,col_colors=row_colors_df)
share_tissue_heatmap.savefig('share_tissue_heatmap.png')









 
for mm in mt:   # tissue-data formation tissue: expression_df
    tissue_matrix=pd.read_csv('/Users/kkwang/mywork/{}_log2_fiter_matrix.tsv'.format(mm),index_col=None,header=0,sep='\t')
    tissue_matrix=tissue_matrix[tissue_matrix['gene_symbol'].isin(specific_genes)]
    tissue_data[mm]=tissue_matrix[tissue_matrix['gene_symbol'].isin(tissue_genes[mm])]
    del(tissue_matrix)  



period_info=pd.read_csv('/Users/kkwang/mywork/syllabus.2.csv',header=0)  # get the period information 
period_info_dict=period_info['col_numbers'].to_dict()        



for mm in mt:    # change the tissue_saple_df into tissue_period_df
    dddf=pd.DataFrame(index=tissue_data[mm].index.tolist(),columns=list(range(len(period_info_dict.keys()))))
    for ab in period_info_dict.items():
        col_range=ab[1].split(',')
        min_col=int(col_range[0])
        max_col=int(col_range[1])
        columns_want=[ col for col in tissue_data[mm].iloc[:,:-1].columns.tolist() if min_col<int(col)<=max_col]
        
        dddf[ab[0]]=tissue_data[mm].loc[:,columns_want].mean(axis=1).values
    dddf['gene_symbol']=tissue_data[mm]['gene_symbol']   
    dddf.fillna(1,inplace=True)
    dddf['tissue']=mm
    
    period_data[mm]=dddf
  

      
res = pd.concat([period_data[mm] for mm in period_data.keys()], axis=0,ignore_index=True) # add all tissue_together
# prepar  color bar for clustermap
res_pal = sns.color_palette("Set2",res['tissue'].unique().size )
#res_pal = sns.cubehelix_palette(n_colors=res['tissue'].unique().size,light=.9, dark=.1, reverse=True,start=1, rot=-2)
res_lut = dict(zip(map(str, res['tissue'].unique()), res_pal))
res_colors = pd.Series(res['tissue']).map(res_lut)





# clustermap of period_genes try_to look whether the gene expression of a tissue can be cluster back into its 
#original tissue 
cluster_map_res = sns.clustermap(res.iloc[:,:-2], cmap="winter",row_colors=res_colors,col_cluster=False)

for label in res['tissue'].unique():
    cluster_map_res.ax_col_dendrogram.bar(0, 0, color=res_lut[label],
                            label=label, linewidth=0)
cluster_map_res.ax_col_dendrogram.legend(loc="center", ncol=6)
cluster_map_res.savefig("period_genes.clustermap.png")





res_T=res.iloc[:,:-2].T
res_corr=res_T.corr()
#res_corr=res_corr.abs()
res_corr[res_corr==1]=0.999
res_corr[res_corr==-1]=-0.999
# fisher z transform z’ = 0.5[ln(1+r) – ln(1-r)]
res_corr_z=0.5 *(np.log(res_corr+1)-np.log(1-res_corr))





# clustermap of the genes_period_corr
cluster_map_res_corr = sns.clustermap(res_corr_z, cmap="coolwarm",row_colors=res_colors,col_cluster=False,col_colors=res_colors)

for label in res['tissue'].unique():
    cluster_map_res_corr.ax_col_dendrogram.bar(0, 0, color=res_lut[label],
                            label=label, linewidth=0)
cluster_map_res_corr.ax_col_dendrogram.legend(loc="center", ncol=6)
cluster_map_res_corr.ax_heatmap.set_yticklabels([])
#cluster_map_res_corr.ax_col_dendrogram.set_visible(False)
cluster_map_res_corr.savefig("period_genes_corr.clustermap.coolwarm.png")






cluster_map_res_corr = sns.clustermap(res_corr_z, cmap="coolwarm",row_colors=res_colors,col_cluster=False,col_colors=res_colors,row_cluster=False)

for label in res['tissue'].unique():
    cluster_map_res_corr.ax_col_dendrogram.bar(0, 0, color=res_lut[label],
                            label=label, linewidth=0)
cluster_map_res_corr.ax_col_dendrogram.legend(loc="center", ncol=6)
cluster_map_res_corr.ax_heatmap.set_yticklabels([])
#cluster_map_res_corr.ax_col_dendrogram.set_visible(False)
cluster_map_res_corr.savefig("period_genes_corr.not.clustermap.png")







# plot of the Z-transformed distribution
fig_pcc=plt.figure('pcc_distribution')
ax = fig_pcc.add_subplot(111)
ax.hist(res_corr_z.values.flatten(),density=True, bins='auto',facecolor='c', alpha=0.75,label='row pcc')
ax.set_xlabel('pearson correlation coefficient Z-trans')
ax.set_yticks([])
ax.set_title('Distribution of pcc')
#ax.set_facecolor('yellow')
ax.grid(False)
fig_pcc.savefig('period_genes_pcc_Z_distribution.png')
fig_pcc.clf()





# prepair for network analysis 
res_corr_z[res_corr_z <= -2]=1
res_corr_z[res_corr_z >= 2] =1
rows, cols = np.where(res_corr_z == 1)
edges = zip(res.index[rows.tolist()],res.index[cols.tolist()])

gr = nx.Graph()
gr.add_edges_from(edges)


degrees=dict(gr.degree())
values=sorted(set(list(degrees.values())))
hist_n=[list(degrees.values()).count(x) for x in values]
plt.figure()
plt.loglog(values,hist_n,'bv')
plt.xlabel('Degree')
plt.ylabel('Number of nodes')







