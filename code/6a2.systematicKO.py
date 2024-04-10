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
from celloracle.applications import Oracle_development_module


# In[3]:


oracle = co.load_hdf5("./chosen_d5base.atac.union_sim_pdt_new.celloracle.oracle")
links = co.load_hdf5(file_path="./chosen_links.atac.union.celloracle.links")
gradient = co.load_hdf5("./oracle.celloracle.gradient")

assert((oracle.adata.obsm[oracle.embedding_name] == gradient.embedding).all())


# In[4]:


links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)


# In[5]:


genes = oracle.active_regulatory_genes
genes


# In[6]:


clusters = ['base','transition','path1', 'proliferation']
cell_idx_path1 = np.where(oracle.adata.obs["paths"].isin(clusters))[0]

index_dictionary = {"Whole_cells": None,
                    "Path1": cell_idx_path1}


# In[7]:


# 0. Define parameters
n_propagation = 3
n_neighbors=200

file_path = "./Systematic_KO_results_all.celloracle.hdf5" 


def pipeline(gene_for_KO):
     
    # 1. Simulate KO
    oracle.simulate_shift(perturb_condition={gene_for_KO: 0},
                                 ignore_warning=True,
                                 n_propagation=3)
    oracle.estimate_transition_prob(n_neighbors=n_neighbors, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    # Do simulation for all conditions.
    for lineage_name, cell_idx in index_dictionary.items():
        
        dev = Oracle_development_module()
        # Load development flow
        dev.load_differentiation_reference_data(gradient_object=gradient)
        # Load simulation result
        dev.load_perturb_simulation_data(oracle_object=oracle, cell_idx_use=cell_idx, name=lineage_name)
        # Calculate inner product
        dev.calculate_inner_product()
        dev.calculate_digitized_ip(n_bins=10)
        
        # Save results in a hdf5 file.
        dev.set_hdf_path(path=file_path) 
        dev.dump_hdf5(gene=gene_for_KO, misc=lineage_name)


# In[8]:


get_ipython().run_cell_magic('time', '', 'for gene in tqdm(genes):\n    pipeline(gene_for_KO=gene)')

