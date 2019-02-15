
# coding: utf-8

# In[1]:


import sys, os
os.nice(15)


# In[2]:


import integer_polyomino.assembly as ipa
import integer_polyomino.gpmap as gp
sys.path.append(os.path.join(os.getcwd(), "..", "src", "integer_polyomino", "scripts"))
import graph_topo
import plotly_utilities as pu


# In[3]:


import plotly.graph_objs as go
from plotly.offline import init_notebook_mode, iplot
import plotly.figure_factory as ff
import cufflinks as cf
import plotly
import pandas as pd
import numpy as np
import seaborn as sns


# In[4]:


get_ipython().magic('matplotlib inline')
init_notebook_mode(connected=True)
cf.set_config_file(offline=True)
sns.set(style="dark", context="talk")


# In[5]:


data_dir = os.path.join(os.getcwd(), "..", "data", "V" + ipa.__version__)
if not os.path.exists(data_dir):
    raise ValueError("Specify valid data directory")


# In[9]:


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


# In[10]:


p["n_genes"] += 1
p["phenotype_builds"] = p["n_genes"] * 50
duplicate_set_metric_file = os.path.join(data_dir, "SetMetricsDuplicate_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.txt".format(**p))


# In[11]:


dfs = pd.read_csv(set_metric_file, sep=" ")
dfd = pd.read_csv(duplicate_set_metric_file, sep=" ")


# In[12]:


metrics = ['srobustness', 'irobustness', 'evolvability',  'complex_diversity', 'diversity', 'robust_evolvability', 'complex_evolvability', 'rare', 'unbound']


# In[13]:


df = dfs.merge(dfd, on=("pIDs"), suffixes=('_sim', '_dup'))
for metric in metrics:
    df[metric + '_rel'] = (df[metric + '_dup'] - df[metric + '_sim'])
    df[metric + '_norm'] = df[metric + '_rel'] / (df[metric + '_sim'] + df[metric + '_dup'])

df["max_size"] = df.pIDs.apply(lambda x: max(eval(x))[0])


# In[14]:


iplot(pu.scatter_metric_size(df, 'srobustness_sim', max_size=15, multi=True, title='srobustness_sim'))


# In[15]:


iplot(pu.dual_metric_size(df, metrics=['srobustness_sim', 'srobustness_dup'], max_size=15, title='Comparison srobustness'))


# In[16]:


iplot(pu.scatter_metrics(df, xMetric='srobustness_norm', yMetric='evolvability_norm', max_size=15, multi=True, title='Comparison srobustness'))


# In[17]:


fig = pu.all_metrics_size(df, metrics, suffixe='_sim', max_size=15, title='Effect of Gene Duplication on GP-map', colors=DEFAULT_PLOTLY_COLORS[0])
fig2 = pu.all_metrics_size(df, metrics, suffixe='_dup', max_size=15, title='Effect of Gene Duplication on GP-map',colors=DEFAULT_PLOTLY_COLORS[1])
for index in range(len(fig2.data)):
    fig.add_trace(fig2.data[index])
iplot(fig)

