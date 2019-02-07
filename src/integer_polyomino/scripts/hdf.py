import pandas as pd
import numpy as np
import graph_topo as gt

def genome_metric_preprocess(df):
    iso_groups(df)
    df['Iso_index'] = df['Iso_index'].astype(np.uint8)
    df['original'] = df['original'].astype('category')
    df['pIDs'] = df['pIDs'].apply(lambda x: str(eval(x))).astype('category')
    df['diversity'] = df['diversity'].astype(np.uint16)
    df['complex_diversity'] = df['complex_diversity'].astype(np.uint16)
    df['neutral_weight'] = df['neutral_weight'].astype('category')
    # df['frequencies'] = df['frequencies'].apply(lambda x: np.array(eval(x), dtype=np.uint8))
    metrics = ['srobustness','irobustness','evolvability', 'robust_evolvability', 'complex_evolvability', 'rare','unbound']
    for metric in metrics:
        df[metric] = df[metric].astype(np.float16)

def iso_groups(df):
    df['Iso_index'] = df['diversity'].astype(np.int16)
    for pID in df.pIDs.unique():
        small_df = df[df.pIDs == pID]
        partition = gt.PartitionPhenotype(list(map(eval, small_df.original.unique())))
        df['Iso_index'].update(small_df['original'].apply(lambda x: partition[str(eval(x))]))

def set_metric_preprocess(df):

#     df['diversity_tracker'] = df['diversity_tracker'].apply(lambda x: np.array(eval(x), dtype=np.int16))
    df['analysed'] = df['analysed'].astype(np.int16)
    df['misclassified'] = df['misclassified'].astype(np.int16)
    df['diversity'] = df['diversity'].astype(np.int16)
#     df['originals'] = df['originals'].apply(lambda x: list(eval(x)))
    df['neutral_size'] = df['neutral_size'].astype(np.int32)
    metrics = ['srobustness', 'irobustness', 'evolvability', 'robust_evolvability', 'complex_evolvability', 'rare','unbound', 'complex_diversity']
    for metric in metrics:
        df[metric] = df[metric].astype(np.float16)

def write_to_hdf(file_path, files, store, kind, overwrite):
    for file_name in files:
        if((overwrite) or not(('/' + file_name[:-4]) in store.keys())):
            df = pd.read_csv(file_path + file_name, sep=' ')
            if(kind == 'set'):
                set_metric_preprocess(df)
            elif(kind == 'genome'):
                genome_metric_preprocess(df)
            store.append(file_name[:-4], df, format='table', data_columns=True)
