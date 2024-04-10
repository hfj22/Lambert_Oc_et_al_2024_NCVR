#!/usr/bin/env python
# coding: utf-8

# In[1]:


import copy
import glob
import importlib
import time
import os
import shutil
import sys
from importlib import reload

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from tqdm.notebook import tqdm


# In[2]:


import celloracle as co
from celloracle.applications import Oracle_development_module, Oracle_systematic_analysis_helper
co.__version__


# In[3]:


#plt.rcParams["font.family"] = "arial"
plt.rcParams["figure.figsize"] = [5,5]
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")
plt.rcParams["savefig.dpi"] = 300
plt.rcParams['pdf.fonttype']=42

get_ipython().run_line_magic('matplotlib', 'inline')


# In[4]:


# Make Oracle development module class to load data
file_path = "Systematic_OE_results_all.celloracle.hdf5"
dev = Oracle_development_module()
dev.set_hdf_path(path=file_path)


# In[5]:


# If we use the function below, we can see information of the saved data
info = dev.get_hdf5_info()

print("Genes\n", info["gene_list"])

print("\nSimulation conditions\n", info["misc_list"])


# In[6]:


# Load one results
dev.load_hdf5(gene="Runx1", misc="Whole_cells")

# Visualize result
dev.visualize_development_module_layout_0(s=5, scale_for_simulation=0.5, s_grid=50,
                                          scale_for_pseudotime=40, vm=0.02)


# In[7]:


# Load data with Oracle_systematic_analysis_helper.
helper = Oracle_systematic_analysis_helper(hdf5_file_path="Systematic_OE_results_all.celloracle.hdf5")


# In[8]:


ps_path1_pos=helper.sort_TFs_by_positive_ip(misc="Path1")
ps_path1_pos.to_csv("OE_all_ps_path1_pos.csv")
ps_path1_pos


# In[9]:


ps_path1_neg=helper.sort_TFs_by_neagative_ip(misc="Path1")
ps_path1_neg.to_csv("OE_all_ps_path1_neg.csv")
ps_path1_neg

