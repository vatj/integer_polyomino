#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys, os
os.nice(0)


# In[2]:


import pandas as pd
import numpy as np
import glob


# In[3]:


import integer_polyomino.assembly as ipa
sys.path.append(os.path.join(os.getcwd(), "..", "src", "integer_polyomino", "scripts"))
import hdf


# In[4]:


data_dir = os.path.join(os.getcwd(), "..", "data", "V" + ipa.__version__)
if not os.path.exists(data_dir):
    raise ValueError("Specify an existing directory")


# In[5]:


# file_name = 'GenomeMetrics_N3_C6_T25_B150_Cx8_J10_D1_S0.txt'
current_files = glob.glob(os.path.join(data_dir,'*.txt'))
file_names = [file.rsplit('/')[-1] for file in current_files]
genome_files = [file_name for file_name in file_names if ("GenomeMetrics_" in file_name)]
duplicate_files = [file_name for file_name in file_names if ("GenomeMetricsDuplicate" in file_name)]
file_hdf = 'Processed_GenomeMetrics.h5'


# In[6]:


genome_files = ['GenomeMetrics_N4_C6_T25_B200_Cx8_J1000_D1_S0.txt',
 'GenomeMetrics_N3_C8_T25_B150_Cx10_J1000_D1_S0.txt',
 'GenomeMetrics_N3_C6_T25_B150_Cx8_J1000_D1_S0.txt',
 'GenomeMetrics_N2_C6_T25_B100_Cx8_J1000_D1_S0.txt']
duplicate_files = ['GenomeMetricsDuplicate_N4_C8_T25_B200_Cx10_J1000_D1_S0.txt',
 'GenomeMetricsDuplicate_N4_C6_T25_B200_Cx8_J1000_D1_S0.txt',
 'GenomeMetricsDuplicate_N3_C6_T25_B150_Cx8_J1000_D1_S0.txt']


# In[7]:


pd.set_option('io.hdf.default_format', 'table')


# In[8]:


with pd.HDFStore(os.path.join(data_dir, file_hdf)) as store:
    hdf.write_to_hdf(data_dir, genome_files, store, 'genome', False)


# In[11]:


with pd.HDFStore(os.path.join(data_dir, file_hdf)) as store:
    hdf.write_to_hdf(data_dir, duplicate_files, store, 'genome', True)


# In[ ]:




