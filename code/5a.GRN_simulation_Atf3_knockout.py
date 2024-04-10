#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import sys

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


# In[3]:


import celloracle as co
co.__version__


# In[4]:


plt.rcParams["figure.figsize"] = [6,10]
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")
plt.rcParams["savefig.dpi"] = 600

get_ipython().run_line_magic('matplotlib', 'inline')


# In[5]:


save_folder = "GRNplots"
os.makedirs(save_folder, exist_ok=True)


# In[6]:


oracle = co.load_hdf5("chosen_d5base.atac.union.celloracle.oracle")
oracle


# In[7]:


links = co.load_hdf5(file_path="chosen_links.atac.union.celloracle.links")


# In[8]:


links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)


# In[9]:


goi = "Atf3"
sc.pl.draw_graph(oracle.adata, color=[goi, oracle.cluster_column_name],
                 layer="imputed_count", use_raw=False, cmap="viridis")


# In[10]:


sc.get.obs_df(oracle.adata, keys=[goi], layer="imputed_count").hist()
plt.show()


# In[11]:


oracle.simulate_shift(perturb_condition={goi: 0.0},
                      n_propagation=3)


# In[12]:


oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)

oracle.calculate_embedding_shift(sigma_corr = 0.05)


# In[13]:


fig, ax = plt.subplots(1, 2,  figsize=[15, 7])

scale = 40
# Show quiver plot
oracle.plot_quiver(scale=scale, ax=ax[0])
ax[0].set_title(f"Perturbation simulation results: {goi} KO")

# Show quiver plot that was calculated with randomized GRN.
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Perturbation simulation with randomized GRNs")

plt.show()


# In[14]:


n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)


# In[15]:


oracle.suggest_mass_thresholds(n_suggestion=12)


# In[16]:


min_mass = 0.004
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)


# In[17]:


min_mass = 0.005
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)


# In[18]:


min_mass = 0.006
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)


# In[19]:


fig, ax = plt.subplots(1, 2,  figsize=[15, 7])

scale_simulation = 0.5
# Show quiver plot
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Perturbation simulation results: {goi} KO")

# Show quiver plot that was calculated with randomized GRN.
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Perturbation simulation with randomized GRNs")

plt.show()


# In[20]:


fig, ax = plt.subplots(figsize=[8, 8])

oracle.plot_cluster_whole(ax=ax, s=10)
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)


# In[21]:


oracle.to_hdf5("chosen_d5base.atac.union_sim.celloracle.oracle")


# In[22]:


import copy
import glob
import time
import shutil
from tqdm.notebook import tqdm


# In[23]:


import celloracle as co
from celloracle.applications import Pseudotime_calculator
co.__version__


# In[24]:


plt.rcParams["figure.figsize"] = [5,5]
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")
plt.rcParams["savefig.dpi"] = 300

get_ipython().run_line_magic('matplotlib', 'inline')


# In[25]:


#pseudotime calculation
pt = Pseudotime_calculator(oracle_object=oracle)


# In[26]:


print("Clustering name: ", pt.cluster_column_name)
print("Cluster list", pt.cluster_list)


# In[27]:


pt.plot_cluster(fontsize=8)


# In[28]:


# These cluster can be classified into path1 or path2

clusters_in_path1 = ['stress','base', 'transition','path1','proliferation']
clusters_in_path2 = ['stress','base','transition', 'path2']

lineage_dictionary = {"path1": clusters_in_path1,
           "path2": clusters_in_path2}
pt.set_lineage(lineage_dictionary=lineage_dictionary)
pt.plot_lineages()


# In[29]:


pt.adata.obs_names[np.flatnonzero(pt.adata.obs['seurat_clusters']  == '0')[0]]


# In[30]:


pt.set_root_cells(root_cells={"path1": "AAACCCACATCTAACG", "path2": "AAACCCACATCTAACG"})
# Check root cell and lineage
pt.plot_root_cells()


# In[31]:


"X_diffmap" in pt.adata.obsm


# In[32]:


pt.obsm_key


# In[33]:


pt.get_pseudotime_per_each_lineage()
pt.plot_pseudotime(cmap="rainbow")


# In[34]:


pt.adata.obs[["Pseudotime"]].head()


# In[35]:


# Add calculated pseudotime to the object
oracle.adata.obs = pt.adata.obs


# In[36]:


oracle.to_hdf5("chosen_d5base.atac.union_sim_pdt.celloracle.oracle")


# In[37]:


fig, ax = plt.subplots(figsize=[6,6])
sc.pl.embedding(adata=oracle.adata, basis=oracle.embedding_name, ax=ax, cmap="rainbow",
                color=["Pseudotime"])


# In[38]:


from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")


# In[39]:


gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)


# In[40]:


gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":4}, plot=True)


# In[41]:


gradient.calculate_gradient()
scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)


# In[42]:


fig, ax = plt.subplots(figsize=[6, 6])
gradient.plot_dev_flow_on_grid(scale=scale_dev, ax=ax)


# In[43]:


gradient.to_hdf5("oracle.celloracle.gradient")


# In[44]:


from celloracle.applications import Oracle_development_module

dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)


# In[45]:


from celloracle.visualizations.config import CONFIG
CONFIG["cmap_ps"] = "coolwarm"


# In[46]:


dev.visualize_development_module_layout_0(s=5,
                                          scale_for_simulation=scale_simulation,
                                          s_grid=50,
                                          scale_for_pseudotime=scale_dev,
                                          vm=0.02)


# In[47]:


# Show inner product score
fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax)


# In[48]:


fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)

