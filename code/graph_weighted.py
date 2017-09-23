import matplotlib.pyplot as plt
import seaborn
import numpy as np
import itertools
import networkx as nx 
import pickle 
from collections import defaultdict
import os
import subprocess
import scipy
import argparse
import datetime
import community 
from networkx.readwrite.gml import literal_stringizer
import sys
from common import parse_config, get_chrom_size

# make the fonts pretty
plt.rc('font', family='serif')

# make the time stamp the same
TIMESTAMP = str(datetime.datetime.now()).replace(' ','-')

def make_weighted_graph(data):
    # remove all rows that have weight of 0
    data = data[data[:,6].astype('float') != 0]
    G = nx.Graph()

    edges1 = data[:,0:2]
    edges2 = data[:,3:5]
    weights = data[:,6].astype('float')
    # make edges hashable
    edges1_h = map(tuple, edges1)
    edges2_h = map(tuple, edges2)

    edge_connections = [(e1,) + (e2,) + (w,) for e1, e2, w in zip(edges1_h, edges2_h, weights)]
    # this will get rid off repeated connections
    G.add_weighted_edges_from(edge_connections)

    #save_adjacency(G, GRAPH_DIR)
    nx.write_gpickle(G, config["GRAPH_NEW_DIR"] + 'all_nodes_edges.adj1Mb.pickle')
    nx.write_gexf(G, config["GRAPH_NEW_DIR"] + 'all_nodes_edges.adj1Mb.gexf')
    nx.write_gml(G, config["GRAPH_NEW_DIR"]  + 'all_nodes_edges.adj1Mb.gml', stringizer = literal_stringizer)
    return G

def write_correlation_distr_file(sparseMat, nodes):
    filename = config["GRAPH_NEW_DIR"] + 'correlation-distr.txt'
    outMat = open(filename, 'w')
    outMat.write("%s\n" %nodes)
    print nodes
    for ii in range(nodes):
        for jj in range(ii + 1, nodes):
            if (jj, ii) in sparseMat:
                outMat.write("%s\n" %sparseMat[(jj, ii)])
            else:
                outMat.write("%s\n" %.0)

    return filename

def weighted_correlation_clustering(G):
    # compute the correlation
    corr_dok = nx.to_scipy_sparse_matrix(G, format = 'dok')    
    keys = corr_dok.keys()
    values = corr_dok.values()
    # create a dictionary
    d = dict(itertools.izip(keys,values))
    # write the dictionary to file
    num_nodes = G.number_of_nodes()
    infile = write_correlation_distr_file(d, num_nodes)
    outname = config["GRAPH_NEW_DIR"] + 'correlation-distr_clusters.txt'

    print infile
    print outname
    # INSERT VOTE/BOEM running function here
    subprocess.call("~/correlation-distr/bin64/chainedSolvers vote boem stats print " + infile + " > " + outname, shell =True)
    # delete infile
    subprocess.call("rm " + infile, shell =True)

    # read in the results
    labels = np.loadtxt(outname, dtype = 'int')

    # need to keep TIMESTAMP so that different version of clusters are not overwritten
    #TIMESTAMP = str(datetime.datetime.now()).replace(' ','-')
    filename = config["GRAPH_NEW_DIR"] + 'weighted_correlation_cluster_members.' + TIMESTAMP + '.pickle'
    clusters =  write_clusters_from_labels(labels, G, filename)

    # write clusters to file, input file is the pickled clusters
    output_file = config["GRAPH_NEW_DIR"] + 'weighted_correlation_cluster_members.' + TIMESTAMP + '.txt'
    write_clusters(clusters, tad_size, output_file)    

def write_clusters_from_labels(labels, G, filename):

    clusters = []
    for clust_num in range(np.min(labels), np.max(labels)):
        # get indices that are labeled this way
        node_ind = np.where(labels == clust_num)[0]
        # get node names that correspond to the indices
        node_names = [G.nodes()[n] for n in node_ind]
        clusters.append(node_names)
    f = open(filename, 'wb')
    pickle.dump(clusters, f)
    f.close()
    return clusters

def write_clusters(clusters, tad_size, output_file):
    
    file = open(output_file, "w")
    for clust in clusters:                                                                     
        if ((len(clust) > 1)):
            line = ''
            chr_list = []
            for node in clust:
                line = line + str(int(node[0])) + '\t' + str(int(node[1])) + '\t' + str(int(node[1]) + tad_size) + '\t'
                chr_list.append(node[0])
            line = line + '\n'
            # check that the cluster is regions between different chromosomes
            if (len(set(chr_list)) != 1):
                file.write(line)
    return


def main():
    global config
    config_fn = sys.argv[1]
    config = parse_config(config_fn)

    data_weighted = np.loadtxt(config["GRAPH_NEW_DIR"]  + 'data.weighted.txt', dtype= 'S12')

    G = make_weighted_graph(data_weighted)

    print 'Number of nodes: ', G.number_of_nodes()
    print 'Number of edges: ', G.number_of_edges()

    weighted_correlation_clustering(G)


if __name__ == "__main__":
    main()
