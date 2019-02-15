
# coding: utf-8

# In[ ]:


import sys, os
os.nice(15)


# In[ ]:


from shutil import copyfile


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

p["n_genes"] = 2
p["low_colour"] = 0
p["high_colour"] = 10
p["threshold"] = 25
p["phenotype_builds"] = p["n_genes"] * 100
p["fixed_table"] = False
p["determinism"] = 1
p["genome_file"] = os.path.join(data_dir, "Genomes_N{n_genes}_C{high_colour}_T{threshold}_B{phenotype_builds}_D{determinism}_S{low_colour}.txt".format(**p))
p["threshold"] /= 100


# In[ ]:


minimal_genomes = gp.MinimalGenomes(**p);
uniques = [];
for pIDs, genomes in gp.MinimalMap(**p).items():
    uniques.extend(graph_topo.TrimTopo(genomes))
    
print(len(uniques), len(minimal_genomes))


# In[ ]:


p["threshold"] = int(p["threshold"] * 100)
p.pop("genome_file", None)
p["table_file"] = os.path.join(data_dir, "PhenotypeTable_N{n_genes}_T{threshold}_D{determinism}.txt".format(**p))
gp.PrintTableFromMap(genomes=uniques, **p)


# In[ ]:


copyfile(os.path.join(data_dir, "PhenotypeTable_N{n_genes}_T{threshold}_D{determinism}.txt".format(**p)), os.path.join(data_dir, "PhenotypeTable_D{determinism}.txt".format(**p)))

