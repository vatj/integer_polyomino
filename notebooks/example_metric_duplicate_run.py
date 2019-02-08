
# coding: utf-8

# In[ ]:


import sys, os
os.nice(15)


# In[ ]:


import integer_polyomino.assembly as ipa
import integer_polyomino.gpmap as gp
sys.path.append(os.path.join(os.getcwd(), "..", "src", "integer_polyomino", "scripts"))
import graph_topo


# In[ ]:


data_dir = os.path.join(os.getcwd(), "..", "data", "V" + ipa.__version__)
if not os.path.exists(data_dir):
    os.makedirs(data_dir)


# In[ ]:


p = dict()

p["n_genes"] = 4
p["low_colour"] = 0
p["high_colour"] = 8
p["threshold"] = 25
p["phenotype_builds"] = p["n_genes"] * 50
p["fixed_table"] = False
p["determinism"] = 1
p["genome_file"] = os.path.join(data_dir, "Genomes_N{n_genes}_C{high_colour}_T{threshold}_B{phenotype_builds}_D{determinism}_S{low_colour}.txt".format(**p))
p["table_file"] = os.path.join(data_dir, "PhenotypeTable_D{determinism}.txt".format(**p))
p["gpmap_file"] = os.path.join(data_dir, "Gpmap_N{n_genes}_C{high_colour}_T{threshold}_B{phenotype_builds}_D{determinism}_S{low_colour}.txt".format(**p))
p["threshold"] /= 100


# In[ ]:


minimal_genomes = gp.MinimalGenomes(**p);
uniques = [];
for pIDs, genomes in gp.MinimalMap(**p).items():
    uniques.extend(graph_topo.TrimTopo(genomes))
    
print(len(uniques), len(minimal_genomes))


# In[ ]:


dup_uniques = []
p.pop("genome_file", None)
p.pop("gpmap_file", None)

for pIDs, genomes in gp.MinimalMap(genomes=gp.GenomesDuplication(uniques), **p).items():
    dup_uniques.extend(graph_topo.TrimTopo(genomes))


# In[ ]:


p["n_genes"] += 1
p["phenotype_builds"] = p["n_genes"] * 50
p["gen_colour"] = p["high_colour"]
p["high_colour"] += 2
p["threshold"] = int(p['threshold'] * 100)
p["n_jiggle"] = 1000
p["dup_aware"] = True
p.pop("genome_file", None)
p["genomes"] = dup_uniques
p["genome_metric_file"] = os.path.join(data_dir, "GenomeMetricsDuplicate_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.txt".format(**p))
p["set_metric_file"] = os.path.join(data_dir, "SetMetricsDuplicate_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}.txt".format(**p))
p["threshold"] /= 100


# In[ ]:


gp.MetricSampling(**p)

