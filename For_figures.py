#!/usr/bin/env python
# coding: utf-8

# In[1]:


import import_ipynb
import Scanpy_functions_v03262021 as sc_pipe
import scvelo as scv
scv.logging.print_version()
import warnings
import scirpy as ir
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
import bbknn
import logging
from sklearn.mixture import GaussianMixture
from scipy.stats     import norm
import glob
import os
import hvplot.pandas
import docx
from docx import Document
from docx.shared import Inches
from docx.shared import Pt
from scipy import sparse
import scanpy.external as sce
import holoviews as hv
import panel as pn
import bokeh
from bokeh.resources import INLINE
import scanorama
import gseapy


# In[2]:


# define sample metadata. Usually read from a file.
exclude_genes = ['Rpl', 'Rps', 'Trav', 'Traj', 'Trbj', 'Trbv','Mrp','Fau','Dap3','Uba52','Ighv', 'Igkv', 'Iglv']

samples = {
'Control-AM_3DE':{},
'Tumor-AM_3DE':{},
}

In_path ='/user/ifrec/liuyuchen/scRNASeq_DATA/scRNASeq_Tanikuchi/F2368_210212_105220_cellranger/'
out_path ='/user/ifrec/liuyuchen/Analysis_Reports/Tanikuchi_scRNA_0510'
multi_sample = True
sc.settings.verbosity = 3
sc.settings.figdir = out_path+'/figures'


# In[3]:


sc.set_figure_params(scanpy=True, dpi=400)
scv.set_figure_params( dpi=400)


# In[4]:


adata=sc.read_h5ad(out_path+'/Analysis0625.h5ad')


# In[5]:


sc.pl.umap(adata, color='leiden', legend_loc='on data',legend_fontsize=10) 


# In[6]:


sc.pl.umap(adata, color='leiden') 


# In[7]:


adata.obs['Sample'].value_counts()


# In[8]:


sc.pl.umap(adata,groups='Tumor-AM_3DE', color='Sample') 


# In[9]:


sc.pl.umap(adata,groups='Control-AM_3DE', color='Sample') 


# In[10]:


sc.tl.embedding_density(adata, groupby='Defined_Types')
sc.pl.embedding_density(adata, groupby='Defined_Types')


# In[11]:


scv.pl.velocity_embedding_stream(adata, basis='umap',color = 'leiden')


# In[12]:


markers= ['Crem', 'Rel', 'Jun']


# In[13]:


sc.pl.umap(adata, color=markers,cmap='OrRd') 


# In[14]:


adata.obs['defined_group']='others'
adata.obs.loc[adata.obs['leiden'].isin(['1','4','8']),'defined_group']='1,4,8'


# In[15]:


sc.tl.rank_genes_groups(adata,groupby='defined_group')


# In[16]:


sc.pl.rank_genes_groups(adata)


# In[17]:


markers= ['Junb','Fos','Smad3','Dach1','Cebpb']


# In[18]:


sc.pl.violin(adata,markers,groupby='defined_group')


# In[19]:


sc.pl.violin(adata[adata.obs['leiden'].isin(['1','4','6','7','8'])],markers,groupby='leiden')


# In[20]:


markers= ['Cd9', 'Ear1', 'S100a1', 'Ly6e', 'S100a6', 'Cd63']


# In[21]:


sc.pl.violin(adata[adata.obs['leiden'].isin(['1','4','6','7','8'])],markers,groupby='leiden')


# In[22]:


df_148 = sc.get.rank_genes_groups_df(adata, group="1,4,8")


# In[23]:


df_148.to_csv('Gene_list.csv')


# In[24]:


sc.tl.rank_genes_groups(adata,groupby='leiden')
sc.tl.filter_rank_genes_groups(adata)


# In[25]:


adata


# In[26]:


rank_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
top_5 = list(set(sum(rank_genes.head(15).values.tolist(),[])))


# In[27]:


top_5 


# In[28]:


sc.tl.dendrogram(adata,groupby='leiden')


# In[29]:


sc.pl.heatmap(adata,top_5, groupby='leiden', use_raw=False,log=True, standard_scale='var', dendrogram=True,show_gene_labels=True )


# In[30]:


sc.pl.heatmap(adata,top_5, groupby='leiden', use_raw=False,log=True, standard_scale='obs', dendrogram=True,show_gene_labels=True )


# In[31]:


sc.pl.rank_genes_groups_heatmap(adata,use_raw=False,log=True,show_gene_labels=True,standard_scale='var')


# In[32]:


adata


# In[33]:


scv.tl.paga(adata, groups='leiden',use_time_prior='latent_time')
scv.pl.paga(adata, basis='umap', size=50, threshold = 0,color = 'Defined_Types',
            min_edge_width=2, node_size_scale=2)


# In[34]:


scv.pl.paga(adata, basis='umap', size=50, dashed_edges=None,
            min_edge_width=2, node_size_scale=2)


# In[35]:


scv.pl.paga(adata, basis='umap', size=50, threshold = 0,color = 'Defined_Types',arrowsize=0,
            min_edge_width=2, node_size_scale=2)


# In[36]:


adata.obs.groupby(['leiden'])['Sample'].value_counts(normalize = True).unstack().plot.bar(stacked=True)
ax = plt.subplot(111)
plt.ylabel('Composition of cells')
plt.xlabel('Leiden cluster')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# In[37]:


markers =['Marco', 'Cd9', 'Ear1', 'S100a1',  'Ly6e', 'S100a6', 'Cd63']


# In[38]:


sc.tl.dendrogram(adata,var_names=markers,groupby='leiden')


# In[39]:


adata=adata[adata.obs['leiden'].isin(['1','4','6','7','8'])]


# In[40]:


sc.tl.dendrogram(adata,groupby='leiden')


# In[41]:


sc.pl.stacked_violin(adata,markers,groupby='leiden',cmap='Reds',standard_scale='var',dendrogram=True,swap_axes=True)


# In[42]:


sc.pl.stacked_violin(adata,markers,groupby='leiden',cmap='Reds',standard_scale='var',use_raw=False,swap_axes=True)


# In[43]:


sc.pl.stacked_violin(adata,markers,groupby='leiden',cmap='Reds',standard_scale='var',swap_axes=True)


# In[44]:


markers =['Marco', 'Cd9', 'Ear1', 'S100a1',  'Ly6e', 'S100a6', 'Cd63','Cd200','Lars2']


# In[45]:


sc.pl.stacked_violin(adata,markers,groupby='leiden',cmap='Reds',standard_scale='var',use_raw=False,swap_axes=True)


# In[46]:


sc.pl.stacked_violin(adata,markers,groupby='leiden',cmap='Reds',standard_scale='var',swap_axes=True)


# In[47]:


sc.pl.violin(adata,'Inhba',use_raw=False,groupby='leiden')


# In[48]:


sc.pl.violin(adata,'Inhba',groupby='leiden')


# In[49]:


for m in markers:
    sc.pl.violin(adata,m,groupby='leiden')


# In[50]:


for m in markers:
    sc.pl.violin(adata,m,use_raw=False,groupby='leiden')


# In[51]:


sc.pl.umap(adata,color ='Inhba',use_raw=False,color_map='OrRd')


# In[54]:


adata


# In[60]:


pd.DataFrame(adata[adata.obs['defined_group']=='1,4,8'][:,'Junb'].X.A).to_csv('Junb_148.csv')


# In[61]:


pd.DataFrame(adata[adata.obs['defined_group']=='others'][:,'Junb'].X.A).to_csv('Junb_others.csv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




