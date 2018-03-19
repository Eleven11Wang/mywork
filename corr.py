import pandas as pd 
import numpy as np
import networkx as nx


all_matrix=pd.read_csv('all_log2_fiter_matrix.tsv',sep='\t',header=0)

genes_symbol=all_matrix['gene_symbol']

all_matrix_T=all_matrix.iloc[:1000,:-1].T
all_log2_corr=all_matrix_T.corr()

all_log2_corr['gene_symbol']=genes_symbol[:1000].tolist()



disease_file=pd.read_csv('spci_df.tsv',sep='\t',index_col=0,header=0)

spci_corr=all_log2_corr[all_log2_corr['gene_symbol'].isin(disease_file['gene_symbol'].tolist())]


spci_gene_symbol=spci_corr['gene_symbol']
spci_corr=spci_corr.iloc[:,:-1]
spci_corr[spci_corr==1]=0.999
spci_corr[spci_corr==-1]=-0.999
# fisher z transform z’ = 0.5[ln(1+r) – ln(1-r)]
spci_corr_z=0.5 *(np.log(spci_corr+1)-np.log(1-spci_corr))

import matplotlib.pyplot as plt

fig_pcc=plt.figure('pcc_distribution')
ax = fig_pcc.add_subplot(111)
ax.hist(spci_corr_z.values.flatten().astype(np.float),density=True, bins='auto',facecolor='c', alpha=0.75,label='row pcc')

ax.set_xlabel('pearson correlation coefficient Z-trans')
ax.set_yticks([])
ax.set_title('Distribution of pcc')


fig_pcc.savefig('hhhhhhhhhhh.png')


spci_corr_z[spci_corr_z <= -2]=1
spci_corr_z[spci_corr_z >= 2] =1

spci_corr_z['gene_symbol']=spci_gene_symbol.tolist()

rows, cols = np.where(spci_corr_z == 1)
edges = zip(spci_corr_z.loc[spci_corr_z.index[rows.tolist()],'gene_symbol'].values,
                                              list(genes_symbol.loc[cols.tolist()]))

gr = nx.Graph()
gr.add_edges_from(edges)



hhhh=sorted(nx.connected_components(gr), key = len, reverse=True)
file=open('hhhhhhhhh.txt','w')
for x in hhhh:
    for xx in x:
        file.write(xx+'\t')
    file.write('\n')
    
file.close()




