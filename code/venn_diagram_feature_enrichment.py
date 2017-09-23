#from matplotlib import pyplot as plt
import numpy as np
#from matplotlib_venn import venn3, venn3_circles
import seaborn
import pandas as pd
import itertools
import networkx as nx
import argparse
from networkx.readwrite.gml import literal_stringizer
import pickle
from common import parse_config, get_chrom_size, get_clusters
import os
import sys



def get_enrichment(features):
    if (os.path.isfile(config["ENRICHMENT_DIR"] + 'enrichment_all.csv')):
        df = pd.read_csv(config["ENRICHMENT_DIR"] + 'enrichment_all.csv', index_col=0)
    else:
        df = pd.DataFrame()
        for feature in features:
            df_feat =  pd.read_csv(config["ENRICHMENT_DIR"] + 'enrichment_' + feature + '.csv', index_col=0)
            df = pd.concat([df, df_feat], axis=1)
    return df

def get_clust_ind_enriched_feature(df):
    # for each feature, get all indices that are enriched (>1)
    # make a dictionary {feature: indices}
    features = df.columns.values
    dict_features_indices = {}
    for feature in features:
        # get all enriched indices for this feature
        ind = set(df[df[feature] > 1].index.values)
        dict_features_indices[feature] = ind
    return dict_features_indices

def get_edges_from_node(combo, node_name):
    edge_list = []
    for feature in combo:
        edge = tuple((node_name, feature))
        edge_list.append(edge)
    return edge_list


def get_combinations(dict_features_indices, features, df):
    # get combinations of different lengths
    node_list = []
    edge_list = []
    
    for i in range(1,len(features)+1):
        combos = list(itertools.combinations(features,i))
        for combo in combos:
            print combo
            # go through each member of the combo
            # start with all indices and intersect them in loop
            set_combo = set(df.index.values)
            for feature in combo:
                ind = dict_features_indices[feature]
                set_combo = set_combo & ind
            # make sure the intersection doesn't contain any of the remaining features
            remaining_features = set(features) - set(combo)
            union_remaining = get_union_features(dict_features_indices, remaining_features)
            set_combo = set_combo - union_remaining

            size = len(set_combo)
            print size
            if (len(combo) == 1):
                node_name = combo[0]
            else:
                node_name = str(combo)
            node_list.append((node_name, {'size' : size}))
            edges = get_edges_from_node(combo, node_name)
            edge_list.append(edges)
    # get rid of list within list
    edge_list = list(itertools.chain.from_iterable(edge_list))
    return node_list, edge_list

def make_graph(node_list, edge_list, ENRICHMENT_NEW_DIR):
    print node_list
    G = nx.Graph()
    scale_factor = 5
    G.add_edges_from(edge_list)
    node_sizes = [ n[1]['size']*scale_factor for n in node_list]
    nodes = [n[0] for n in node_list]
    nx.draw_networkx(G, pos=nx.spring_layout(G), nodelist = nodes, node_size = node_sizes)
    #nx.write_gml(G, ENRICHMENT_NEW_DIR + 'all_nodes_edges.adj1Mb.gml', stringizer = literal_stringizer)
    plt.axis('off')
    plt.show()

def venn_diagram(node_list):
    plt.figure(figsize=(4,4))
    subset = (31, 13, 32, 2, 0, 11, 178)
    v = venn3(subsets=subset, set_labels = ('H3K27ac', 'H3K4me3', 'RNAPII'))
    plt.show()
    subset = (244, 267, 224, 212, 201, 202, 192)
    v = venn3(subsets=subset, set_labels = ('H3K4me3', 'H3K4me1', 'POLR2A'))
    plt.show()
    subset = (195, 212, 118, 171, 66, 26, 24)
    v = venn3(subsets=subset, set_labels = ('H3K27me3', 'POLR2A', 'H3K9me3'))
    plt.show()

def filter_clusters(filename, filename_new, df):
    # load clusters
    clusters = get_clusters(filename)
    # active clusters
    enriched_regions = df[(df['RNAPII'] > 1) & ((df['H3K9me3'] <= 1))].dropna().index
    print enriched_regions
    s = pd.Series(enriched_regions)
    s.to_csv(config["ENRICHMENT_DIR"] + 'active.txt', header=False, index= False)    
    clusters_filt = [clust for i,clust in enumerate(clusters) if i in enriched_regions]
    f = open(filename_new, 'wb')
    pickle.dump(clusters_filt, f)
    f.close()

def filter_clusters_inactive(filename, filename_new, df, ENRICHMENT_NEW_DIR):
    # load clusters
    clusters = get_clusters(filename)
    # active clusters
    enriched_regions = df[(df['POLR2A'] <= 1) & ((df['H3K9me3'] > 1))].dropna().index    
    clusters_filt = [clust for i,clust in enumerate(clusters) if i in enriched_regions]
    s = pd.Series(enriched_regions)
    s.to_csv(ENRICHMENT_NEW_DIR + 'inactive.txt', header=False, index= False) 
    f = open(filename_new, 'wb')
    pickle.dump(clusters_filt, f)
    f.close()

def get_union_features(dict_features_indices, features):
    # return union set of all passed features
    set_union = set()
    for feature in features:
        set_union = set.union(set_union, dict_features_indices[feature])
    return set_union

def get_num_onefeature_only(dict_features_indices, features):
    for feature in features:
        # get all indices that only belong to this feature
        ind_feature = dict_features_indices[feature]
        # get union of the remaining features
        remaining_feat = [f for f in features if f != feature]
        set_union = get_union_features(dict_features_indices, remaining_feat)
        diff = ind_feature - set_union
        print feature
        print len(diff)

def get_intersection_features(dict_features_indices, features):
    set_intersection = get_union_features(dict_features_indices, features)
    for feature in features:
        set_intersection = set_intersection.intersection(dict_features_indices[feature])
    #print len(set_intersection)
    fname = config["ENRICHMENT_DIR"] + 'intersection.txt'
    f = open(fname, 'w')
    f.write('\n'.join(map(str, list(set_intersection))))
    f.close()
    return set_intersection

def get_intersection_inactive_union(set_intersect, dict_features_indices, feature):
    set_union = set.union(set_intersect, dict_features_indices[feature])
    return len(set_union)

def get_intersection_mutiple_feat(dict_features_indices, features):
    # generate combinations of features

    for feature in features:
        # get all indices that only belong to this feature
        ind_feature = dict_features_indices[feature]
        # get union of the remaining features
        remaining_feat = [f for f in features if f != feature]
        set_union = get_union_features(dict_features_indices, remaining_feat)
        diff = ind_feature - set_union
        print feature
        print len(diff)

def write_enriched_feature_indices(dict_features_indices, features):
    # write tab separated list of indices for each feature
    for feature in features:
        ind = list(dict_features_indices[feature])
        f = open(config["ENRICHMENT_DIR"] + feature + '.txt', 'w')
        f.write('\n'.join(map(str, ind)))
        f.close()


def main():
    global config
    config_fn = sys.argv[1]
    config = parse_config(config_fn)

    filename = config["GRAPH_NEW_DIR"] + config["filename"]

    features = np.array(['H3K27ac', 'H3K4me3', 'H3K27me3', 'H3K9me3', 'RNAPII'])
    df = get_enrichment(features)
    features = df.columns.values
    dict_features_indices = get_clust_ind_enriched_feature(df)
    #features = np.array(['H3K27ac', 'H3K4me3' , 'RNAPII'])
    #features = np.array(['H3K9ac', 'POLR2A', 'H3K4me1', 'H3K4me3' , 'H3K36me3', 'H3K9me3'])

    #features = np.array(['H3K9ac', 'POLR2A', 'H3K9me3', 'H3K36me3'])
    #features = np.array(['H3K4me3', 'H3K4me1', 'POLR2A'])
    #features = np.array(['H3K27me3', 'POLR2A', 'H3K9me3'])
    #features = np.array(['POLR2A', 'H3K9me3'])
    write_enriched_feature_indices(dict_features_indices, features)
    node_list, edge_list = get_combinations(dict_features_indices, features, df)
    #make_graph(node_list, edge_list, ENRICHMENT_NEW_DIR)

    set_intersect = get_intersection_features(dict_features_indices, features)
    l = get_intersection_inactive_union(set_intersect, dict_features_indices, 'H3K9me3')
    ##print l

    filename_new = config["GRAPH_ACTIVE_DIR"] + config["filename"]
    filter_clusters(filename, filename_new, df)
    #filter_clusters_inactive(filename, filename_new, df, ENRICHMENT_NEW_DIR)

    ##print len(set_intersect), len(df)

if __name__ == "__main__":
    main()
