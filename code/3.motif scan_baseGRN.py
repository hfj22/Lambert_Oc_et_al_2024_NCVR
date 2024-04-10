#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import seaborn as sns

import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm


# In[3]:


from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object


# In[4]:


get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")
get_ipython().run_line_magic('matplotlib', 'inline')

plt.rcParams['figure.figsize'] = (15,7)
plt.rcParams["savefig.dpi"] = 600


# In[5]:


ref_genome = "mm10"
genome_installation = ma.is_genome_installed(ref_genome=ref_genome)
print(ref_genome, "installation: ", genome_installation)


# In[6]:


if not genome_installation:
    import genomepy
    genomepy.install_genome(ref_genome, "UCSC")
else:
    print(ref_genome, "is installed.")


# In[7]:


peaks = pd.read_csv("processed_peak_file.csv", index_col=0)
peaks.head()


# In[8]:


def decompose_chrstr(peak_str):
    """
    Args:
        peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'

    Returns:
        tuple: chromosome name, start position, end position
    """

    *chr_, start, end = peak_str.split("_")
    chr_ = "_".join(chr_)
    return chr_, start, end

from genomepy import Genome

def check_peak_foamat(peaks_df, ref_genome):
    """
    Check peak fomat.
     (1) Check chromosome name.
     (2) Check peak size (length) and remove sort DNAs (<5bp)

    """

    df = peaks_df.copy()

    n_peaks_before = df.shape[0]

    # Decompose peaks and make df
    decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
    df_decomposed = pd.DataFrame(np.array(decomposed))
    df_decomposed.columns = ["chr", "start", "end"]
    df_decomposed["start"] = df_decomposed["start"].astype(np.int)
    df_decomposed["end"] = df_decomposed["end"].astype(np.int)

    # Load genome data
    genome_data = Genome(ref_genome)
    all_chr_list = list(genome_data.keys())


    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])


    # Filter peaks with invalid chromosome name
    n_threshold = 5
    df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]

    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

    # Data counting
    n_invalid_length = len(lengths[lengths < n_threshold])
    n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
    n_peaks_after = df.shape[0]

    #
    print("Peaks before filtering: ", n_peaks_before)
    print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
    print("Peaks with invalid length: ", n_invalid_length)
    print("Peaks after filtering: ", n_peaks_after)

    return df


# In[9]:


peaks = check_peak_foamat(peaks, ref_genome)


# In[10]:


tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome)


# In[11]:


get_ipython().run_cell_magic('time', '', '\ntfi.scan(fpr=0.02,\n         motifs=None,           verbose=True)\n\n# Save tfinfo object\ntfi.to_hdf5(file_path="test1.celloracle.tfinfo")')


# In[12]:


tfi.scanned_df.head()


# In[13]:


tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=10)
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)


# In[14]:


df = tfi.to_dataframe()
df.head()


# In[15]:


folder = "baseGRN"
os.makedirs(folder, exist_ok=True)

# Save result as a dataframe
df = tfi.to_dataframe()
df.to_parquet(os.path.join(folder, "base_GRN_dataframe.parquet"))


# In[16]:


td = tfi.to_dictionary(dictionary_type="targetgene2TFs")
save_as_pickled_object(td, os.path.join(folder, "TFinfo_targetgene2TFs.pickled"))


# In[ ]:




