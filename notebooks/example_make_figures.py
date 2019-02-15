
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
import plotly.io as pio
import pandas as pd
import numpy as np
from plotly.colors import DEFAULT_PLOTLY_COLORS
import cufflinks as cf
import itertools


# In[ ]:


init_notebook_mode(connected=True)
cf.set_config_file(offline=True)

if not os.path.exists(os.getcwd() + '/../figures'):
    os.mkdir(os.getcwd() + '/../figures')
figure_dir = os.getcwd() + '/../figures/'

if not os.path.exists(os.getcwd() + '/../data/V1'):
    os.mkdir(os.getcwd() + '/../data/V1')
data_dir = os.getcwd() + "/../data/V1/"


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
p["table_file"] = data_dir + "PhenotypeTable_D{determinism}.txt".format(**p)
genome_metric_file = "GenomeMetrics_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}".format(**p)
set_metric_file = data_dir + "SetMetrics_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.txt".format(**p)
file_hdf = 'Processed_GenomeMetrics.h5'
gpmap_dir = figure_dir + 'gpmap_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}/'.format(**p)


# In[ ]:


p["n_genes"] += 1
p["phenotype_builds"] = p["n_genes"] * 50
max_size = 10
duplicate_set_metric_file = data_dir + "SetMetricsDuplicate_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.txt".format(**p)
duplicate_genome_metric_file = "GenomeMetricsDuplicate_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}".format(**p)
p["n_genes"] -= 1
p["phenotype_builds"] = p["n_genes"] * 50


# In[ ]:


dfs = pd.read_csv(set_metric_file, sep=" ")
dfs_d = pd.read_csv(duplicate_set_metric_file, sep=" ")


# In[ ]:


metrics = ['srobustness', 'irobustness', 'evolvability',  'complex_diversity', 'diversity', 'robust_evolvability', 'complex_evolvability', 'rare', 'unbound']
suffixes = ['_sim', '_dup', '_norm']


# In[ ]:


df = dfs.merge(dfs_d, on=("pIDs"), suffixes=('_sim', '_dup'))
for metric in metrics:
    df[metric + '_rel'] = (df[metric + '_dup'] - df[metric + '_sim'])
    df[metric + '_norm'] = df[metric + '_rel'] / (df[metric + '_sim'] + df[metric + '_dup'])

df["max_size"] = df.pIDs.apply(lambda x: max(eval(x))[0])


# In[ ]:


if not os.path.exists(gpmap_dir):
    os.mkdir(gpmap_dir)


# In[ ]:


current_path = gpmap_dir + 'deterministic_vs_multi/'
if not os.path.exists(current_path):
    os.mkdir(current_path)

for metric, suffixe in itertools.product(metrics, suffixes):
    title = 'Deterministic v.s. Multi ' + metric + suffixe + ' metric in GPmap N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}'.format(**p)
    filename = 'deterministic_vs_multi_' + metric + suffixe + '_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.pdf'.format(**p)
    pio.write_image(pu.scatter_metric_size(df, metric + suffixe, max_size=max_size, multi=True, title=title), file=current_path + filename)


# In[ ]:


current_path = gpmap_dir + 'simple_vs_duplicate/'
if not os.path.exists(current_path):
    os.mkdir(current_path)

for metric in metrics:
    title = 'Simple v.s. Duplicate ' + metric + ' metric in GPmap N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}'.format(**p)
    filename = 'simple_vs_duplicate_' + metric + '_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.pdf'.format(**p)
    pio.write_image(pu.dual_metric_size(df, metrics=[metric + '_sim', metric + '_dup'], max_size=max_size, title=title, symbol=['circle', 'x']), file=current_path + filename)


# In[ ]:


distribution_dir = gpmap_dir + 'phenotypes_distribution/'
if not os.path.exists(distribution_dir):
    os.mkdir(distribution_dir)


# In[ ]:


with pd.HDFStore(data_dir + file_hdf, mode='r') as store:
    pIDs_simple = store.select(genome_metric_file, 'columns == pIDs')
    pIDs_dup = store.select(duplicate_genome_metric_file, 'columns == pIDs')


# In[ ]:


current_path = distribution_dir + 'simple/'
if not os.path.exists(current_path):
    os.mkdir(current_path)

for pIDs in pIDs_simple.pIDs.unique():
    if not os.path.exists(current_path + pIDs):
        os.mkdir(current_path + pIDs)
    for metric in metrics:
        title = 'Distribution of ' + metric + ' metric for phenotype ' + pIDs + '<br> in GPmap N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}'.format(**p)
        filename = '/distribution_' + metric + '_for_phenotype_' + pIDs + '_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.pdf'.format(**p)
        layout = go.Layout(title=title)
        pio.write_image(go.Figure(data=pu.distribution_metric_phenotype_trace(pID=pIDs, file_name=genome_metric_file, metric=metric, hdf=data_dir+file_hdf), layout=layout), file=current_path + pIDs + filename)


# In[ ]:


current_path = distribution_dir + 'duplicate/'
if not os.path.exists(current_path):
    os.mkdir(current_path)

for pIDs in pIDs_dup.pIDs.unique():
    if not os.path.exists(current_path + pIDs):
        os.mkdir(current_path + pIDs)
    for metric in metrics:
        title = 'Distribution of ' + metric + ' metric for phenotype ' + pIDs + '<br>in GPmap-duplicate N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}'.format(**p)
        filename = '/distribution_' + metric + '_for_phenotype_' + pIDs + '_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.pdf'.format(**p)
        layout = go.Layout(title=title)
        pio.write_image(go.Figure(data=pu.distribution_metric_phenotype_trace(pID=pIDs, file_name=duplicate_genome_metric_file, metric=metric, hdf=data_dir+file_hdf), layout=layout), file=current_path + pIDs + filename)


# In[ ]:


current_path = distribution_dir + 'comparison/'
if not os.path.exists(current_path):
    os.mkdir(current_path)

for pIDs in pIDs_dup.pIDs.unique():
    if not os.path.exists(current_path + pIDs):
        os.mkdir(current_path + pIDs)
    for metric in metrics:
        title = 'Distribution of ' + metric + ' metric for phenotype ' + pIDs + '<br>in GPmap N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}'.format(**p)
        filename = '/distribution_' + metric + '_for_phenotype_' + pIDs + '_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.pdf'.format(**p)
        layout = go.Layout(title=title)
        traces = pu.distribution_metric_phenotype_dup_trace(pID=pIDs, genome_file=genome_metric_file, duplicate_file=duplicate_genome_metric_file, metric=metric, hdf=data_dir+file_hdf)
        pio.write_image(go.Figure(data=traces, layout=layout), file=current_path + pIDs + filename)

