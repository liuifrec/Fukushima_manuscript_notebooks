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


adata_list = sc_pipe.data_IO(samples,In_path,IG=False)


# In[4]:


adata = sc_pipe.unify_value(adata_list)


# In[5]:


get_ipython().run_cell_magic('time', '', 'adata = sc_pipe.qc_and_preprocess(adata, out_path,feature_tag=  {\'mt\':\'mt-\',\'ribo\':("Rps","Rpl"),\'hb\':("^Hb[^(P)]")},multi_sample=True,organism = \'mmusculus\')\nadata = sc_pipe.feature_selection(adata,out_path, exclude_genes)\n')


# In[6]:


get_ipython().run_cell_magic('time', '', 'adata = sc_pipe.clustering(adata,samples,resol=1)\n')


# In[7]:


scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualizatio


# In[8]:


adata.obs['Sample'].value_counts()


# In[9]:


adata =  sc_pipe.differential(adata,out_path, multi_sample, top_n= 2,organism="mmusculus")


# In[10]:


for c in ['Clusters', '_X', '_Y']:
    del adata.obs[c]


# In[11]:


adata.obs['Condition'] = 'Tumor'
adata.obs.loc[adata.obs['Sample']=='Control-AM_3DE','Condition']='Control'


# In[12]:


Factors = ['Inhba']
adata.obs['Defined_Types'] =adata.obs['Condition']
for f in Factors:
    delta, var, cut_off,up_bound,low_bound= sc_pipe.Thrshold_by_Gaussian_with_plot(adata, f,out_path, 'All_samples', sparse = True)
    data = adata[:, f].X.A
    data =  np.interp(data, (data.min(), data.max()), (0, 10))
    adata.obs[f+'_binary'] = np.where(data > cut_off, f+'+', f+'-')
    adata.obs['Defined_Types'] =adata.obs['Defined_Types'].astype(str)+'_'+adata.obs[f+'_binary']


# In[13]:


sc.tl.embedding_density(adata, groupby='Inhba_binary')
sc.pl.embedding_density(adata, groupby='Inhba_binary')


# In[14]:


sc.tl.embedding_density(adata, groupby='Defined_Types')
sc.pl.embedding_density(adata, groupby='Defined_Types')


# In[15]:


All_annotation = list(adata.obs.columns)
adata.obs.index.name = 'cell'
sc.external.exporting.cellbrowser(adata, 'Tanikuchi_Inhba_0514', '0514_all_cells',  annot_keys=All_annotation, cluster_field='leiden',html_dir='Cell_browser_Tanikuchi')


# In[16]:


adata = adata[~adata.obs['leiden'].isin(['14','15'])]


# In[17]:


sub = adata[adata.obs['Inhba_binary']=='Inhba+'].copy()


# In[18]:


sub = sub[sub.obs['leiden'].isin(['6','7','8','4','1','9','0'])].copy()


# In[19]:


adata


# In[20]:


adata


# In[21]:


scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


# In[22]:


scv.tl.recover_dynamics(adata,n_jobs=30)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])
#scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=0.5)
scv.pl.velocity_embedding_stream(adata, basis='umap',color = 'leiden')


# In[23]:


scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80)
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='leiden', n_convolve=100)
scv.tl.paga(adata, groups='leiden',use_time_prior='latent_time')
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[24]:


sc.pl.umap(adata, color= 'Inhba',use_raw=False, color_map= 'OrRd')


# In[25]:


scv.pl.velocity_graph(adata,color='leiden', threshold=.1)


# In[26]:


scv.pl.velocity_graph(adata,color='latent_time', threshold=.1)


# In[27]:


sc.tl.embedding_density(adata, groupby='Defined_Types')
sc.pl.embedding_density(adata, groupby='Defined_Types')


# In[28]:


scv.tl.paga(adata, groups='leiden',use_time_prior='latent_time')
scv.pl.paga(adata, basis='umap', size=50, threshold = 0,color = 'Defined_Types',
            min_edge_width=2, node_size_scale=2)


# In[29]:


specific= adata[adata.obs['leiden'].isin(['6','7','8','4','1'])].copy()
scv.tl.paga(specific, groups='leiden',use_time_prior='latent_time')
scv.pl.paga(specific, basis='umap', size=50, threshold = 0,color = 'Defined_Types',
            min_edge_width=2, node_size_scale=2)


# In[30]:


sc.tl.embedding_density(specific, groupby='Defined_Types')
sc.pl.embedding_density(specific, groupby='Defined_Types')


# In[31]:


scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


# In[32]:


scv.pl.scatter(specific, color='latent_time', color_map='gnuplot', size=80)
top_genes = specific.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(specific, var_names=top_genes, sortby='latent_time', col_color='leiden', n_convolve=100)


# In[33]:


sc.tl.dendrogram(specific,groupby='leiden')


# In[34]:


specific =  sc_pipe.differential(specific,out_path, multi_sample, top_n= 2,organism="mmusculus")


# In[35]:


All_annotation = list(specific.obs.columns)
specific.obs.index.name = 'cell'
sc.external.exporting.cellbrowser(specific, 'Tanikuchi_sub1_0514', '0514_sub1',  annot_keys=All_annotation, cluster_field='leiden',html_dir='Cell_browser_Tanikuchi')


# In[36]:


specific


# In[37]:


sc.tl.rank_genes_groups(specific,groupby='leiden',reference='7')


# In[38]:


sc.pl.rank_genes_groups(specific, key='rank_genes_groups')
rank_genes = pd.DataFrame(specific.uns['rank_genes_groups']['names'])

rank_pavlues = pd.DataFrame(specific.uns['rank_genes_groups']['pvals'])
rank_pavlues.rename(columns=lambda x: x+'_pvals', inplace=True)
rank_fold = pd.DataFrame(specific.uns['rank_genes_groups']['logfoldchanges'])
rank_fold.rename(columns=lambda x: x+'_logfoldchanges', inplace=True)
pd.concat([rank_genes, rank_pavlues,rank_fold], axis=1, join='inner').to_csv(out_path+'/DE_genes_against_cluster7.csv')


# In[39]:


for c in ['1','4','6','8']:
    q= rank_genes[c].head(200).tolist()
    sc.queries.enrich(q, org="mmusculus",gprofiler_kwargs={'no_evidences':False}).to_csv(out_path+'/cluster_'+c+'_against_cluster7_enrichment.csv')
    


# In[40]:


sub3 =  adata[adata.obs['leiden'].isin(['1','4','7','8'])].copy()
sub3.obs['Groups'] = '1,4,8'
sub3.obs.loc[sub3.obs['leiden']=='7','Groups'] = '7'
sc.tl.rank_genes_groups(sub3,groupby='Groups')
sc.pl.rank_genes_groups(sub3, key='rank_genes_groups')

rank_genes = pd.DataFrame(sub3.uns['rank_genes_groups']['names'])

rank_pavlues = pd.DataFrame(sub3.uns['rank_genes_groups']['pvals'])
rank_pavlues.rename(columns=lambda x: x+'_pvals', inplace=True)
rank_fold = pd.DataFrame(sub3.uns['rank_genes_groups']['logfoldchanges'])
rank_fold.rename(columns=lambda x: x+'_logfoldchanges', inplace=True)
pd.concat([rank_genes, rank_pavlues,rank_fold], axis=1, join='inner').to_csv(out_path+'/DE_genes_two_groups.csv')


# In[41]:


for c in ['1,4,8','7']:
    q= rank_genes[c].head(200).tolist()
    sc.queries.enrich(q, org="mmusculus",gprofiler_kwargs={'no_evidences':False}).to_csv(out_path+'/Group_'+c+'_enrichment.csv')
    


# In[42]:


sub1 =  adata[adata.obs['leiden'].isin(['6','7'])].copy()
scv.tl.recover_dynamics(sub1,n_jobs=44)
scv.tl.paga(sub1, groups='leiden',use_time_prior='latent_time')
scv.pl.paga(sub1, basis='umap', size=50, threshold = 0,color = 'Defined_Types',
            min_edge_width=2, node_size_scale=2)

sub2 = adata[adata.obs['leiden'].isin(['6','8','4','1'])].copy()
scv.tl.recover_dynamics(sub2,n_jobs=44)
scv.tl.paga(sub2, groups='leiden',use_time_prior='latent_time')
scv.pl.paga(sub2, basis='umap', size=50, threshold = 0,color = 'Defined_Types',
            min_edge_width=2, node_size_scale=2)


# In[43]:


scv.pl.scatter(sub1, color='latent_time', color_map='gnuplot', size=80)
top_genes1 = sub1.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.pl.heatmap(sub1, var_names=top_genes1, sortby='latent_time', col_color='leiden', n_convolve=100)
scv.pl.scatter(sub2, color='latent_time', color_map='gnuplot', size=80)
top_genes2 = sub2.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.pl.heatmap(sub2, var_names=top_genes2, sortby='latent_time', col_color='leiden', n_convolve=100)


# In[44]:


common= list(set(top_genes1)&set(top_genes2))


# In[45]:


path67only = list(set(top_genes1)-set(top_genes2))


# In[46]:


path6841only  = list(set(top_genes2)-set(top_genes1))


# In[47]:


sc.queries.enrich(path67only, org="mmusculus",gprofiler_kwargs={'no_evidences':False}).to_csv(out_path+'/path_6_7_enrichment.csv')


# In[48]:


pd.DataFrame(path67only).to_csv(out_path+'/path_6_7_genes.csv')


# In[49]:


sc.queries.enrich(path6841only, org="mmusculus",gprofiler_kwargs={'no_evidences':False}).to_csv(out_path+'/path_6_8_4_1_enrichment.csv')


# In[50]:


pd.DataFrame(path6841only).to_csv(out_path+'/path_6_8_4_1_genes.csv')


# In[51]:


sc.queries.enrich(common, org="mmusculus",gprofiler_kwargs={'no_evidences':False}).to_csv(out_path+'/two_path_common_enrichment.csv')


# In[52]:


scv.pl.velocity(sub1, path67only)


# In[53]:


scv.pl.velocity(sub2,path6841only)


# In[54]:


pd.DataFrame(top_genes1).to_csv(out_path+'/genes1.csv')
pd.DataFrame(top_genes2).to_csv(out_path+'/genes2.csv')


# In[55]:


specific2= adata[adata.obs['leiden'].isin(['6','2','0','3','5','9','12'])].copy()
scv.tl.paga(specific2, groups='leiden',use_time_prior='latent_time')
scv.pl.paga(specific2, basis='umap', size=50, threshold = 0,color = 'Defined_Types',
            min_edge_width=2, node_size_scale=2)


# In[56]:


sc.tl.dendrogram(specific2,groupby='leiden')
specific2 =  sc_pipe.differential(specific2,out_path, multi_sample, top_n= 2,organism="mmusculus")


# In[57]:


sc.tl.embedding_density(specific2, groupby='Defined_Types')
sc.pl.embedding_density(specific2, groupby='Defined_Types')


# In[58]:


scv.pl.scatter(specific2, color='latent_time', color_map='gnuplot', size=80)
top_genes = specific2.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(specific2, var_names=top_genes, sortby='latent_time', col_color='leiden', n_convolve=100)


# In[59]:


All_annotation = list(specific2.obs.columns)
specific2.obs.index.name = 'cell'
sc.external.exporting.cellbrowser(specific2, 'Tanikuchi_sub2_0514', '0514_sub2',  annot_keys=All_annotation, cluster_field='leiden',html_dir='Cell_browser_Tanikuchi')


# In[60]:


cp = list(adata.uns['Defined_Types_colors'])
adata.obs.groupby('leiden')['Defined_Types'].value_counts(normalize = True).unstack().plot.bar(stacked=True,color =cp)
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# In[61]:


scv.pl.proportions(sub,groupby='Sample',layers=['ambiguous', 'spliced', 'unspliced'])


# In[62]:


sub.obs['Condition'].value_counts()


# In[63]:


scv.pp.moments(sub, n_pcs=30, n_neighbors=30)


# In[64]:


sc.tl.leiden(sub,key_added='re-clustered')
sc.tl.paga(sub, groups='re-clustered')
sc.pl.paga(sub,node_size_power= 1.0,threshold=0.1 ,edge_width_scale = 0.1) 
sc.tl.umap(sub,  init_pos=sc.tl._utils.get_init_pos_from_paga(sub), maxiter=100,min_dist=0.1, spread=4.0)
sc.pl.umap(sub, color='leiden', legend_loc='on data',legend_fontsize=10) 
sc.pl.umap(sub, color='re-clustered', legend_loc='on data',legend_fontsize=10) 


# In[65]:


sc.tl.rank_genes_groups(sub, 're-clustered',use_raw=False)
sc.pl.rank_genes_groups(sub)
sc.tl.dendrogram(sub,groupby='re-clustered')
sc.pl.rank_genes_groups_dotplot(sub,n_genes=2,use_raw=False)
sc.pl.rank_genes_groups_stacked_violin(sub,n_genes=2,use_raw=False)
sc.pl.rank_genes_groups_heatmap(sub,n_genes=2,use_raw=False)


# In[66]:


sc.tl.embedding_density(sub, groupby='Defined_Types')
sc.pl.embedding_density(sub, groupby='Defined_Types')


# In[67]:


scv.tl.recover_dynamics(sub,n_jobs=30)
scv.tl.velocity(sub, mode='dynamical')
scv.tl.velocity_graph(sub)
scv.tl.score_genes_cell_cycle(sub)
scv.pl.scatter(sub, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])
scv.tl.velocity_confidence(sub)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(sub, c=keys, cmap='coolwarm', perc=[5, 95])


# In[68]:


scv.pl.scatter(sub, color_gradients=[ 'G2M_score'], smooth=False)


# In[69]:


sc.pl.umap(sub,color = 'phase')


# In[70]:


scv.pl.velocity_embedding_stream(sub, basis='umap',color = 're-clustered')


# In[71]:


scv.pl.velocity_embedding_stream(sub, basis='umap',color = 'leiden')


# In[72]:


scv.tl.latent_time(sub,root_key='0')
scv.pl.scatter(sub, color='latent_time', color_map='gnuplot', size=80)
top_genes = sub.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(sub, var_names=top_genes, sortby='latent_time', col_color='re-clustered', n_convolve=100)


# In[73]:


scv.tl.paga(sub, groups='re-clustered',use_time_prior='latent_time')
scv.pl.paga(sub, basis='umap', size=50, threshold = 0,color = 'Defined_Types',
            min_edge_width=2, node_size_scale=2)


# In[74]:


scv.tl.paga(sub, groups='re-clustered',use_time_prior='latent_time',minimum_spanning_tree =False)
scv.pl.paga(sub, basis='umap', size=50, threshold = 0,color = 'Defined_Types',
            min_edge_width=2, node_size_scale=2)


# In[75]:


cp = list(sub.uns['Defined_Types_colors'])
sub.obs.groupby('re-clustered')['Defined_Types'].value_counts(normalize = True).unstack().plot.bar(stacked=True,color =cp)
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# In[76]:


adata.write(out_path+'/Analysis0625.h5ad')


# In[ ]:




