import networkx as nx
from collections import defaultdict, deque, Counter
import itertools

import numpy as np

def Transform_Graph_From_List(tile_kit):
    graph_kit=nx.MultiDiGraph()
    graph_kit.add_nodes_from(range(len(tile_kit)))

    ## Add edges for internal structure in clockwise orientation
    for i_edge in range(int(len(tile_kit)/4)):
        slices=deque([0,1,2,3])
        for _ in range(4):
            graph_kit.add_edge(i_edge*4+slices[0],i_edge*4+slices[1])
            slices.rotate(-1)

    ## Add edges for graph connections
    for index,face in enumerate(tile_kit):
        i=0
        if face:
            while(Interaction_Matrix(face) in tile_kit[index+i:]):
                i+=tile_kit[index+i:].index(Interaction_Matrix(face))
                graph_kit.add_edge(index,index+i)#,color='r')
                graph_kit.add_edge(index+i,index)#,color='r')
                i+=1

    return graph_kit

def PartitionPhenotype(genotypes):
    """
    genotypes: list of genotypes that all have the same phenotype
    returns a dict of lists, where the list contain isomorphic genotypes
    """
    index = 0
    graph_clusters = dict()
    graph_clusters[str(genotypes[index])] = str(index)
    network_graphs = dict()
    network_graphs[str(index)] = Transform_Graph_From_List(genotypes[index])

    for genotype in genotypes[1:]:
        ref_graph = Transform_Graph_From_List(genotype)
        
        for key, comp_graph in network_graphs.items():
            if nx.is_isomorphic(ref_graph, comp_graph):
                graph_clusters[str(genotype)] = key
                break
        else:
            index += 1
            graph_clusters[str(genotype)] = str(index)
            network_graphs[str(index)] = ref_graph

    return graph_clusters


def TrimTopo(genotypes):
    """
    genotypes: list of genotypes that all have the same phenotype
    returns a single representant for each non-isomorphic assembly_graph
    """
    uniques = [genotypes[0]]
    network_graphs = [Transform_Graph_From_List(genotypes[0])]

    for genotype in genotypes[1:]:
        ref_graph=Transform_Graph_From_List(genotype)
<<<<<<< HEAD
        ref_graph_chiral = Transform_Graph_From_List(genotype[::-1])
        for i,comp_graph in enumerate(network_graphs):
            if (nx.is_isomorphic(ref_graph,comp_graph)  || nx.is_isomorphic(ref_graph_chiral, comp_graph)):
=======
        ref_graph_chiral=Transform_Graph_From_List(genotype[::-1])
        for i,comp_graph in enumerate(network_graphs):
            if (nx.is_isomorphic(ref_graph,comp_graph) or nx.is_isomorphic(ref_graph_chiral, comp_graph)):
>>>>>>> 64d799422f5d92d3d503ca153681faf7537cd776
                break
        else:
            uniques.append(genotype)
            network_graphs.append(ref_graph)

    return uniques


def Interaction_Matrix(colour):
    return colour if colour <=0 else  (1-colour%2)*(colour-1)+(colour%2)*(colour+1)


def StripIsomorphisms(file_in):
    tile_kits=[[int(i) for i in line.rstrip('\n').split()] for line in open(file_in)]
    assembly_graphs=list(zip(list(range(len(tile_kits))),[Transform_Graph_From_List(tile_kit) for tile_kit in tile_kits]))

    unique_assembly_graphs_indices=[]

    while len(assembly_graphs)>1:
        unique_assembly_graphs_indices.append(assembly_graphs[0][0])
        assembly_graphs[:]=[assembly_graph for assembly_graph in assembly_graphs[1:] if not nx.is_isomorphic(assembly_graph[1],assembly_graphs[0][1])]

    with open(file_in.rstrip('.txt')+'_Out.txt', 'w') as outfile:
        for index in unique_assembly_graphs_indices:
            genotype= ' '.join(map(str,tile_kits[index]))+'\n'
            outfile.write(genotype)

from multiprocessing import Pool
def StripInParallel(runs):
    pool = Pool(processes=4)
    pool.map_async(StripIsomorphisms,["A2_T20_C200_N500_Mu0.003125_O25_K15000_Run{}_Genotype".format(r) for r in runs])
    pool.close()


def Load_Tiles(fileN):
    Raw_Topologies_Input = [line.rstrip('\n') for line in open(fileN)]
    Raw_Topologies = [[int(face) for face in line.split()] for line in Raw_Topologies_Input]
    topology_dict=defaultdict(list)

    for tile_kit in Raw_Topologies:
        zeroes=Order_By_Zeroes(tile_kit)
        topology_dict[zeroes].append(tile_kit)

    return topology_dict


def Trim_Topologies(fin):
    TD=Load_Tiles(fin)
    unique_C=0
    for key in reversed(sorted(TD.keys())):
        del_count=0
        print("on ",key, len(TD[key]))
        for i,comp1 in enumerate(TD[key]):
            G1=Transform_Graph_From_List(comp1)
            del_count=0
            for j,comp2 in enumerate(TD[key][i+1:]):
                G2=Transform_Graph_From_List(comp2)
                if nx.is_isomorphic(G1,G2):
                    #print i+1+j-del_count,len(value),i,j,del_count
                    del TD[key][i+1+j-del_count]
                    del_count+=1
            unique_C+=1
    print(unique_C)
    with open(fin.rstrip('.txt')+'_Iso.txt', 'w') as outfile:
        for zero_set in list(TD.values()):
            for genotype in zero_set:
                genotype_str= ' '.join(map(str,genotype))+'\n'
                outfile.write(genotype_str)




def Order_By_Zeroes(tile_kit):
    z_Map=defaultdict(int)
    for i in range(int(len(tile_kit)/4)):
        tile=tile_kit[i*4:i*4+4]
        if 0 in tile:
            z_Map[tile.count(0)]+=1
        else:
            z_Map[0]+=1

    return (z_Map[0],z_Map[1],z_Map[2],z_Map[3],z_Map[4])


    while(len(z_Map)>0):
        new_tile_kit.extend(tile_kit[max(z_Map, key=z_Map.get)*4:max(z_Map, key=z_Map.get)*4+4])
        del z_Map[max(z_Map, key=z_Map.get)]

    return new_tile_kit


def Enumerate_Topology(tile_kit):
    queue_kit=[deque(tile) for tile in tile_kit]


    for tile_rot in range(4**len(tile_kit)):
        for tile_index in range(len(tile_kit)-1,-1,-1):
            if tile_rot%(4**tile_index)==0:
                queue_kit[tile_index].rotate(-1)
                yield itertools.permutations(queue_kit)
                break


def cycle_tile(tile):
    for rotation in range(4):
        tile.rotate(-1)
        yield tile
