#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys, os
import pandas as pd
import numpy as np
os.nice(15)


# In[2]:


import integer_polyomino.assembly as ipa
import integer_polyomino.gpmap as gp
sys.path.append(os.path.join(os.getcwd(), "..", "src", "integer_polyomino", "scripts"))
import graph_topo


# In[3]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[4]:


data_dir = os.path.join(os.getcwd(), "..", "data", "V" + ipa.__version__)
if not os.path.exists(data_dir):
    os.makedirs(data_dir)


# In[26]:


parameters = dict()

parameters["threshold"] = 0.25
parameters["phenotype_builds"] = 120
parameters["fixed_table"] = False
parameters["determinism"] = 1
parameters["table_file"] = os.path.join(data_dir, "PhenotypeTable_D{determinism}.txt".format(**parameters))


# In[27]:


ipa.AssemblePlasticGenotype([0,0,1,2,0,0,0,2], **parameters)


# In[32]:


ipa.AssemblePlasticGenotypeFrequency([0,0,1,2,0,0,0,2], **parameters);


# In[30]:


ipa.AssemblePlasticGenotypes([[0,0,0,1,0,0,2,2]], **parameters);


# In[35]:


p = dict()

p["n_genes"] = 2
p["low_colour"] = 0
p["high_colour"] = 6
p["threshold"] = 25
p["phenotype_builds"] = p["n_genes"] * 50
p["fixed_table"] = False
p["determinism"] = 1
p["genome_file"] = os.path.join(data_dir, "Genomes_N{n_genes}_C{high_colour}_T{threshold}_B{phenotype_builds}.txt".format(**p))
p["table_file"] = os.path.join(data_dir, "table_N{n_genes}_C{high_colour}_T{threshold}_B{phenotype_builds}.txt".format(**p))
p["gpmap_file"] = os.path.join(data_dir, "gpmap_N{n_genes}_C{high_colour}_T{threshold}_B{phenotype_builds}.txt".format(**p))
p["threshold"] /= 100


# In[36]:


minimal_genomes = gp.MinimalGenomes(**p);


# In[37]:


gpmap = gp.MinimalMap(**p);


# In[38]:


gp.GenotypeNeighbourhood([0,0,1,2], low_colour=-1, high_colour=2);


# In[39]:


uniques = [];
for pIDs, genomes in polyo.MinimalMap(**p).items():
    uniques.extend(graph_topo.TrimTopo(genomes))


# In[40]:


p["n_genes"] = 2
p["gen_colour"] = 6
p["high_colour"] = 8
p["n_jiggle"] = 100
p["dup_aware"] = False
p.pop("genome_file", None)
p["genomes"] = uniques
p["genome_metric_file"] = os.path.join(data_dir, "GenomeMetrics_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}.txt".format(**p))
p["set_metric_file"] = os.path.join(data_dir, "SetMetrics_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}.txt".format(**p))


# In[41]:


gp.MetricSampling(**p)


# In[43]:


metric_struct = gp.MetricNeighbourhood([0,0,0,1,0,0,0,2], **p)
neighbourhood = pd.Series([str(x) for x in gp.PhenotypeNeighbourhood([0,0,0,1,0,0,0,2], **p)], dtype="category")


# In[45]:


neighbourhood.value_counts()


# In[24]:





# In[ ]:




