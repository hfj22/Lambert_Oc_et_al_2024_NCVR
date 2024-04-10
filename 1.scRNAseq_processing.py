#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


# In[3]:


get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")
get_ipython().run_line_magic('matplotlib', 'inline')

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300


# In[4]:


#D5 data including raw data and metadata
adata = sc.read('d5_raw.h5ad')
adata


# In[5]:


sc.pp.filter_genes(adata, min_counts=1)
adata


# In[6]:


sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')


# In[7]:


filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                              flavor='seurat',
                                              n_top_genes=4000,
                                              log=True)


# In[8]:


adata_TF = adata[:, ["Id3", "Eno1", "Mef2c", "Hif1a", "Nfia", "Hmga1", "Ebf1", "Rora", "Klf9", "Jund", "Foxp1", "Thra", "Prrx1", "Nfix", "Stat3","Twist1","Runx2"]]
adata_TF


# In[9]:


# Subset the genes
adata_HVG = adata[:, filter_result.gene_subset]
adata_HVG


# In[10]:


adata = ad.concat([adata_TF, adata_HVG], axis=1, merge="same", uns_merge="same")
adata


# In[11]:


sc.pp.normalize_per_cell(adata)


# In[12]:


adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()

sc.pp.log1p(adata)
sc.pp.scale(adata)


# In[13]:


sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata)


# In[14]:


sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.diffmap(adata)

sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')


# In[15]:


# Current cluster name
cluster_list = adata.obs.seurat_clusters.unique()
cluster_list


# In[16]:


# Cluster anottation dictionary
annotation = {"path1":[4,3],
              "path2": [6,10],
              "transition":[2],
              "base":[0, 5, 7, 8],
             "stress": [1],
             "proliferation": [9]}

annotation_rev = {}
for i in cluster_list:
    for k in annotation:
        if int(i) in annotation[k]:
            annotation_rev[i] = k
            
annotation_rev


# In[17]:


adata.obs["paths"] = [annotation_rev[i] for i in adata.obs.seurat_clusters]


# In[18]:


sc.tl.paga(adata, groups='paths')


# In[19]:


plt.rcParams["figure.figsize"] = [6, 4.5]


# In[20]:


sc.pl.paga(adata)


# In[21]:


sc.tl.draw_graph(adata, init_pos='paga', random_state=123)


# In[22]:


sc.pl.draw_graph(adata, color='paths', legend_loc='on data')


# In[23]:


sc.pl.draw_graph(adata, color='seurat_clusters', legend_loc='on data')


# In[24]:


plt.rcParams["figure.figsize"] = [4.5, 4.5]


# In[25]:


markers = {"proliferation":["Mki67"],
           "synthetic state":["Spp1", "Tnfrsf11b", "Col8a1", "Mgp", "Timp1", "Ly6a"],
            "contractile state":["Myh11", "Tagln"],
            "stress":["Atf3", "Fos", "Cryab","Hspa1b"],
           "other":["Stat3", "Twist1", "Runx2", "Runx1"]
            }

for cell_type, genes in markers.items():
    print(f"marker genes of {cell_type}")
    sc.pl.draw_graph(adata, color=genes, cmap="viridis", use_raw=False, ncols=4)
    plt.show()


# In[26]:


adata.write_h5ad("chosen_d5_preprocessed new.h5ad")
adata

