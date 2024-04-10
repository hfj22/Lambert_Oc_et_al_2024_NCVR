#!/usr/bin/env python
# coding: utf-8

# In[1]:


import copy
import glob
import time
import os
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from tqdm.notebook import tqdm


# In[2]:


import celloracle as co
from celloracle.applications import Pseudotime_calculator
co.__version__


# In[3]:


plt.rcParams["figure.figsize"] = [5,5]
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")
plt.rcParams["savefig.dpi"] = 300

get_ipython().run_line_magic('matplotlib', 'inline')


# In[4]:


oracle = co.load_hdf5("chosen_d5base.atac.union_sim_pdt.celloracle.oracle")
oracle.adata.uns['paths_colors'] = ['#a9561e', '#ff6600', '#bbbbbb', '#ff0000', '#dddddd', '#ffcc00']


# In[5]:


#pseudotime calculation
pt = Pseudotime_calculator(oracle_object=oracle)


# In[6]:


print("Clustering name: ", pt.cluster_column_name)
print("Cluster list", pt.cluster_list)


# In[7]:


# Check data
pt.plot_cluster(fontsize=8)


# In[8]:


# These cluster can be classified into path1

clusters_in_path1 = ['stress','base', 'transition','path1','proliferation']

# Make dictionary
lineage_dictionary = {"path1": clusters_in_path1}

# Inpur lineage information into pseudotime object
pt.set_lineage(lineage_dictionary=lineage_dictionary)

# Visualize lineage information
pt.plot_lineages()


# In[9]:


#find the cell with highest myh11 expression
max_expression_index = pt.adata[:, "Myh11"].X.argmax()
cell_name_with_max_expression = pt.adata.obs_names[max_expression_index]
print(cell_name_with_max_expression)


# In[10]:


#show the cluster in that
print(pt.adata.obs.loc[cell_name_with_max_expression,'seurat_clusters'])


# In[11]:


pt.set_root_cells(root_cells={"path1": cell_name_with_max_expression, "path2": cell_name_with_max_expression})
# Check root cell and lineage
pt.plot_root_cells()


# In[12]:


# Check diffusion map data.
"X_diffmap" in pt.adata.obsm


# In[13]:


pt.obsm_key


# In[14]:


# Calculate pseudotime
pt.get_pseudotime_per_each_lineage()

# Check results
pt.plot_pseudotime(cmap="rainbow")


# In[15]:


# Check result
pt.adata.obs[["Pseudotime"]].head()


# In[16]:


# Add calculated pseudotime data to the oracle object
oracle.adata.obs = pt.adata.obs


# In[17]:


# Save updated oracle object
oracle.to_hdf5("chosen_d5base.atac.union_sim_pdt_new.celloracle.oracle")


# In[ ]:




