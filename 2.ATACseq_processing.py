#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

import seaborn as sns

import os, sys, shutil, importlib, glob
from tqdm import tqdm_notebook as tqdm

get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300


# In[2]:


from celloracle import motif_analysis as ma


# In[3]:


file_path_of_bed_file = "All_VSMC_MACS2PE_Peaks.bed"
bed = ma.read_bed(file_path_of_bed_file)
print(bed.shape)
bed.head()


# In[4]:


peaks = ma.process_bed_file.df_to_list_peakstr(bed)
peaks


# In[5]:


tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="mm10")

tss_annotated.tail()


# In[6]:


peak_id_tss = ma.process_bed_file.df_to_list_peakstr(tss_annotated)
tss_annotated = pd.DataFrame({"peak_id": peak_id_tss,
                              "gene_short_name": tss_annotated.gene_short_name.values})
tss_annotated = tss_annotated.reset_index(drop=True)
print(tss_annotated.shape)
tss_annotated.head()


# In[7]:


tss_annotated.to_csv("processed_peak_file.csv")

