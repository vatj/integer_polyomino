
# coding: utf-8

# In[ ]:


import sys, os
os.nice(15)


# In[ ]:


import pandas as pd
import numpy as np
import glob


# In[ ]:


import integer_polyomino.assembly as ipa
sys.path.append(os.path.join(os.getcwd(), "..", "src", "integer_polyomino", "scripts"))
import hdf


# In[ ]:


data_dir = os.path.join(os.getcwd(), "..", "data", "V" + ipa.__version__)
if not os.path.exists(data_dir):
    raise ValueError("Specify an existing directory")


# In[ ]:


# file_name = 'GenomeMetrics_N3_C6_T25_B150_Cx8_J10_D1_S0.txt'
current_files = glob.glob(os.path.join(data_dir,'*.txt'))
file_names = [file.rsplit('/')[-1] for file in current_files]
file_hdf = 'Processed_GenomeMetrics.h5'


# In[ ]:


pd.set_option('io.hdf.default_format', 'table')


# In[ ]:


genome_files = [file_name for file_name in file_names if ("GenomeMetrics_" in file_name)]
duplicate_files = [file_name for file_name in file_names if ("GenomeMetricsDuplicate" in file_name)]


# In[ ]:


with pd.HDFStore(os.path.join(data_dir, file_hdf)) as store:
    hdf.write_to_hdf(data_dir, genome_files, store, 'genome', False)


# In[ ]:


with pd.HDFStore(os.path.join(data_dir, file_hdf)) as store:
    hdf.write_to_hdf(data_dir, duplicate_files, store, 'genome', False)

