
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


# %matplotlib inline
init_notebook_mode(connected=True)
cf.set_config_file(offline=True)
sns.set(style="dark", context="talk")


# In[5]:


data_dir = os.path.join(os.getcwd(), "..", "data", "V" + ipa.__version__)
if not os.path.exists(data_dir):
    raise ValueError("Specify valid data directory")

figure_dir = os.path.join(os.getcwd(), '..', 'figures')
if not os.path.exists(figure_dir):
    os.mkdir(figure_dir)


# In[6]:


metrics = ['srobustness', 'irobustness', 'evolvability',  'complex_diversity', 'diversity', 'robust_evolvability', 'complex_evolvability', 'rare', 'unbound']
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
p["n_genes"] += 1
p["phenotype_builds"] = p["n_genes"] * 50
duplicate_set_metric_file = os.path.join(data_dir, "SetMetricsDuplicate_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.txt".format(**p))
dfs = pd.read_csv(set_metric_file, sep=" ")
dfd = pd.read_csv(duplicate_set_metric_file, sep=" ")


# In[7]:


df = dfs.merge(dfd, on=("pIDs"), suffixes=('_sim', '_dup'))
for metric in metrics:
    df[metric + '_rel'] = (df[metric + '_dup'] - df[metric + '_sim'])
    df[metric + '_norm'] = df[metric + '_rel'] / (df[metric + '_sim'] + df[metric + '_dup'])

df["max_size"] = df.pIDs.apply(lambda x: max(eval(x))[0])


# In[8]:


import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from itertools import product

def VisualiseSingleShape(shape,ax=None,corner=(0,1),extent=1,add_direction=False):
#     cols=['darkgreen','royalblue','firebrick','goldenrod','mediumorchid']
    cols = sns.color_palette('muted')
    hatchs=['//','.','\\\\','O']
    ar_offsets={0:(0,-.25,0,.25),1:(-.25,0,.25,0),2:(0,.25,0,-.25),3:(.25,0,-.25,0)}
    rescale=False
    if ax is None:
        rescale=True
        fig = plt.figure()
        ax = fig.gca()

    dx=shape[0]
    dy=shape[1]
    max_d=max(dx,dy)
    extentS=extent/max_d
    for i,j in product(range(dx),range(dy)):
        if(shape[2+i+j*dx]):
            
            coords=corner[0]+(i/max_d)*extent, corner[1]-((j+1)/max_d)*extent

            ax.add_patch(Rectangle(coords, extentS, extentS, facecolor=cols[(shape[2+i+j*dx]-1)//4], edgecolor='black', fill=True,lw=5,transform=ax.transAxes))
#             ax.add_patch(Rectangle(coords, extentS, extentS , edgecolor='black',fill=False,lw=4.5,transform=ax.transAxes))
            theta=(shape[2+i+j*dx]-1)%4;
            if add_direction:  
                ax.arrow(coords[0]+extentS/2, coords[1]+extentS/2, ar_offsets[theta][2]*extentS, ar_offsets[theta][3]*extentS, head_width=0.01*extent, head_length=0.025*extent, fc='k', ec='k',transform=ax.transAxes)

    if rescale:
        maxB=max(dx,dy)
        #ax.set_xlim([-maxB/2-0.5,maxB/2+0.5])
        #ax.set_ylim([-maxB/2-0.5,maxB/2+0.5])
        ax.set_aspect('equal')
        ax.set_axis_off()
        plt.show(block=False)
        
    return fig


# In[9]:


fig_format = 'pdf'


# In[10]:


fig = VisualiseSingleShape([1,2,5,1])
fig.savefig(os.path.join(figure_dir, 'heterodimer') + '.' + fig_format, format=fig_format)


# In[11]:


fig = VisualiseSingleShape([2,2,1,5,5,1])
fig.savefig(os.path.join(figure_dir, 'heterotetramer') + '.' + fig_format, format=fig_format)


# In[12]:


fig = VisualiseSingleShape([2,2,1,1,1,1])
fig.savefig(os.path.join(figure_dir, 'homotetramer')+ '.' + fig_format, format=fig_format)


# In[13]:


fig = VisualiseSingleShape([2,2,1,0,5,1])
fig.savefig(os.path.join(figure_dir, 'heterotrimer')+ '.' + fig_format, format=fig_format)


# In[14]:


# VisualiseSingleShape([3,2,0,1,0,1,5,1])


# In[15]:


# VisualiseSingleShape([3,3,0,1,0,1,5,1,0,1,0])


# In[16]:


fig = VisualiseSingleShape([4,4,0,1,0,0,0,5,5,1,1,5,5,0,0,0,1,0])
fig.savefig(os.path.join(figure_dir, 'heteroctomer')+ '.' + fig_format, format=fig_format)


# In[17]:


fig = VisualiseSingleShape([1,1,1])
fig.savefig(os.path.join(figure_dir, 'bluemonomer') + '.' + fig_format, format=fig_format)


# In[18]:


fig = VisualiseSingleShape([1,1,5])
fig.savefig(os.path.join(figure_dir, 'orangemonomer') + '.' + fig_format, format=fig_format)

