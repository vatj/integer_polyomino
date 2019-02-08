
# coding: utf-8

# In[ ]:


import sys, os
os.nice(15)


# In[ ]:


import integer_polyomino.assembly as ipa
import integer_polyomino.gpmap as gp
sys.path.append(os.path.join(os.getcwd(), "..", "src", "integer_polyomino", "scripts"))
import graph_topo
import plotly_utilities as pu


# In[ ]:


import plotly.graph_objs as go
from plotly.offline import init_notebook_mode, iplot
import plotly.figure_factory as ff
import plotly
import pandas as pd
import numpy as np
import seaborn as sns


# In[ ]:


get_ipython().magic('matplotlib inline')
init_notebook_mode(connected=True)
sns.set(style="dark", context="talk")


# In[ ]:


data_dir = os.path.join(os.getcwd(), "..", "data", "V" + ipa.__version__)
if not os.path.exists(data_dir):
    raise ValueError("Specify an existing directory")


# In[ ]:


p = dict()

p["n_genes"] = 3
p["low_colour"] = 0
p["gen_colour"] = 6
p["high_colour"] = 8
p["threshold"] = 25
p["phenotype_builds"] = p["n_genes"] * 50
p["fixed_table"] = False
p["determinism"] = 1
p["n_jiggle"] = 1000
p["table_file"] = os.path.join(data_dir, "PhenotypeTable_D{determinism}.txt".format(**p))
set_metric_file = os.path.join(data_dir, "SetMetrics_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.txt".format(**p))

genome_metric_file = "GenomeMetrics_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}".format(**p)
file_hdf = 'Processed_GenomeMetrics.h5'


# In[ ]:


with pd.HDFStore(data_dir + file_hdf, mode='r') as store:
    dfg = store[genome_metric_file]


# In[ ]:


iplot(pu.pretty_genome_pie(dfg, row))


# In[ ]:


iplot(pu.pretty_genome_hbar(dfg, row))


# In[ ]:


iplot(pu.metric_subplots(dfg, [row], p, 'group'))


# In[ ]:


iplot(pu.distribution_metrics_phenotype(pID='{(3, 0)}', file_name = genome_metric_file, hdf=os.path.join(data_dir, file_hdf)))

