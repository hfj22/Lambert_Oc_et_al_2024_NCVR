#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


# In[3]:


import celloracle as co
co.__version__


# In[4]:


get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")
get_ipython().run_line_magic('matplotlib', 'inline')

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300


# In[5]:


save_folder = "GRNplots"
os.makedirs(save_folder, exist_ok=True)


# In[6]:


adata = sc.read("chosen_d5_preprocessed new.h5ad")
adata


# In[7]:


df=pd.read_parquet("./baseGRN/base_GRN_dataframe.parquet")

base_GRN = df
# Check data
base_GRN.head()


# In[8]:


oracle = co.Oracle()


# In[9]:


print("metadata columns :", list(adata.obs.columns))
print("dimensional reduction: ", list(adata.obsm.keys()))


# In[10]:


adata.X = adata.layers["raw_count"].copy()

oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="paths",
                                   embedding_name="X_draw_graph_fa")


# In[11]:


oracle.import_TF_data(TF_info_matrix=base_GRN)


# In[12]:


oracle.perform_PCA()

plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
print(n_comps)
n_comps = min(n_comps, 50)


# In[13]:


n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")


# In[14]:


k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")


# In[15]:


oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)


# In[16]:


oracle.to_hdf5("chosen_d5base.atac.union.celloracle.oracle")


# In[17]:


#GRN calculation
sc.pl.draw_graph(oracle.adata, color="paths")


# In[18]:


get_ipython().run_cell_magic('time', '', 'links = oracle.get_links(cluster_name_for_GRN_unit="paths", alpha=10,\n                         verbose_level=10, test_mode=False)')


# In[19]:


links.links_dict.keys()


# In[24]:


links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)


# In[25]:


links.filtered_links["path1"].to_csv(f"chosen_filtered_GRN_atac_for_path1.csv")
links.filtered_links["path2"].to_csv(f"chosen_filtered_GRN_atac_for_path2.csv")
links.filtered_links["base"].to_csv(f"chosen_filtered_GRN_atac_for_base.csv")
links.filtered_links["stress"].to_csv(f"chosen_filtered_GRN_atac_for_stress.csv")
links.filtered_links["proliferation"].to_csv(f"chosen_filtered_GRN_atac_for_proliferation.csv")
links.filtered_links["transition"].to_csv(f"chosen_filtered_GRN_atac_for_transition.csv")


# In[26]:


plt.rcParams["figure.figsize"] = [9, 4.5]


# In[27]:


save_folder = "GRNplots"
os.makedirs(save_folder, exist_ok=True)


# In[28]:


links.plot_degree_distributions(plot_model=True,save=f"{save_folder}/degree_distribution/")


# In[29]:


plt.rcParams["figure.figsize"] = [6, 4.5]


# In[30]:


links.get_network_score()


# In[32]:


links.merged_score.to_csv(f"chosen_network_scores.csv")


# In[33]:


links.to_hdf5(file_path="chosen_links.atac.union.celloracle.links")


# In[40]:


links.plot_scores_as_rank(cluster="path1", n_gene=30, save=f"{save_folder}/ranked_score")


# In[43]:


links.plot_score_comparison_2D(value="degree_centrality_all",
                               cluster1="path1", cluster2="base",
                               percentile=98,
                               save=f"{save_folder}/score_comparison")

