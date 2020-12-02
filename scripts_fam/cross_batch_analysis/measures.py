import pandas as pd
import sys
import os
import glob
from collections import defaultdict
from sklearn.preprocessing import normalize
from sklearn.metrics import pairwise_distances
import numpy as np


def exact_mc_perm_test(xs, ys, nmc):
    n, k = len(xs), 0
    diff = np.abs(np.mean(xs) - np.mean(ys))
    zs = np.concatenate([xs, ys])
    for j in range(nmc):
        np.random.shuffle(zs)
        k += diff <= np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
    return k / nmc


def counts_distance(files):
    distance_in_family = []
    distance_others = []
    for file in files:
        df = pd.read_csv(file, index_col=0)
        df_copy = df.drop('pop_average', axis=1)
        columns = []
        for col, source, fam in zip(df_copy.columns.values, df_copy.loc['seq_source'].values, df_copy.loc['family'].values):
            columns.append(f'{col} - {source} - {fam}')
        df_copy = df_copy.drop(['seq_source', 'family', 'relationship'], axis=0)

        cols = df_copy.columns
        df_copy[cols] = df_copy[cols].apply(pd.to_numeric)

        contig_df = pd.read_csv('/scratch/groups/dpwall/personal/chloehe/unmapped_reads/ref_genome/contig_df.csv', index_col=0)
        index = []
        for idx in df_copy.index.values:
            index.append(contig_df.loc[idx].short_des)

        # threshold = 10
        # df_copy = df_copy[df_copy.ge(threshold).any(axis=1)]
        df_scaled = normalize(df_copy)

        df_scaled = pd.DataFrame(df_scaled, columns=columns, index=index)
        
        distance = pairwise_distances(df_scaled.T)

        family = defaultdict(list)
        for idx, col in enumerate(columns):
            family[col.split(' - ')[2]].append(idx)

        fam_samples = set()
        for i in range(len(columns)):
            fam_samples.add(i)

        for idx, sample in enumerate(columns):
            fam = sample.split(' - ')[2]
            rel = family[fam].copy()
            others = list(fam_samples - set(rel))
            rel.remove(idx)
            for r in rel:
                distance_in_family.append(distance[idx, r])
            for o in others:
                distance_others.append(distance[idx, o])
                
    return distance_in_family, distance_others
