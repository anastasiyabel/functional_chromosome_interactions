import json
import os
import pickle
from collections import defaultdict
import pandas as pd
import numpy as np
import itertools

def parse_config(config_fn):
    config = json.load(open(config_fn))
    if 'HIC_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['HIC_DIR']) + '/'
        config['HIC_DIR'] = directory
        create_dir(directory)
    if 'HIC_FILT_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['HIC_FILT_DIR']) + '/'
        config['HIC_FILT_DIR'] = directory
        create_dir(directory)
    if 'HIC_OUT_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['HIC_OUT_DIR']) + '/'
        config['HIC_OUT_DIR'] = directory
        create_dir(directory)
    if 'HIC_NEW_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['HIC_NEW_DIR']) + '/'
        config['HIC_NEW_DIR'] = directory
        create_dir(directory)
    if 'INTRACHROM_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['INTRACHROM_DIR']) + '/'
        config['INTRACHROM_DIR'] = directory
        create_dir(directory)
    if 'CHIPSEQ_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['CHIPSEQ_DIR']) + '/'
        config['CHIPSEQ_DIR'] = directory
        create_dir(directory)
    if 'CHIPSEQ_OUT_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['CHIPSEQ_OUT_DIR']) + '/'
        config['CHIPSEQ_OUT_DIR'] = directory
        create_dir(directory)
    if 'GENOME_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['GENOME_DIR']) + '/'
        config['GENOME_DIR'] = directory
        create_dir(directory)
    if 'GRAPH_NEW_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['GRAPH_NEW_DIR']) + '/'
        config['GRAPH_NEW_DIR'] = directory
        create_dir(directory)
    if 'GRAPH_ACTIVE_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['GRAPH_ACTIVE_DIR']) + '/'
        config['GRAPH_ACTIVE_DIR'] = directory
        create_dir(directory)
    if 'ENRICHMENT_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['ENRICHMENT_DIR']) + '/'
        config['ENRICHMENT_DIR'] = directory
        create_dir(directory)
    if 'TFBS_DIR' in config:
        directory = os.path.join(os.path.dirname(config_fn), config['TFBS_DIR']) + '/'
        config['TFBS_DIR'] = directory
        create_dir(directory)
    return config

def create_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_chrom_size(chrom, directory):
        # get chrom size from known file
        lengths_filename = directory + 'chrom_hg19.sizes'
        with open(lengths_filename) as f:
                for line in f.readlines():
                        d = line.rstrip().split('\t')
                        if (d[0] == ('chr' + str(chrom))):
                                chrom_size = int(d[1])
        return chrom_size

def get_clusters(filename):
    f = open(filename, 'rb')
    clusters = pickle.load(f)    

    clusters_filt = []
    for i, clust in enumerate(clusters):
        dic_chr = defaultdict(float)
        if (len(clust) > 1):

            for node in clust:
                if (len(clust) != 1):
                    chr = node[0]
                    dic_chr[chr] = dic_chr[chr] + 1
            
        # check that all chromosomes are different
        if (len(dic_chr.keys()) > 1):
            clusters_filt.append(clust)
    return clusters_filt

def get_num_bases_genome(config):
    chrom_size_file =  config["GENOME_DIR"] + 'chrom_hg19.sizes'
    chrom_df = pd.read_csv(chrom_size_file, header = None, names = ['chr', 'size'], sep = '\t')
    num_bases_genome = chrom_df.iloc[0:len(config["chrs"]),1].sum()
    return num_bases_genome

def get_clustered_regions(filename, config):
    clusters_nonsingleton = get_clusters(filename)
    clustered = np.asarray(list(itertools.chain.from_iterable(clusters_nonsingleton)))
    df_clustered = pd.DataFrame(clustered, columns = ['chr', 'start'])
    df_clustered['start'] = df_clustered['start'].astype('int')
    df_clustered['stop'] = df_clustered['start'] + config["HIC_RESOLN"]
    return df_clustered