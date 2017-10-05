import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import sys
import numpy as np
import seaborn
import pandas as pd
import pickle
import itertools
import os.path
import pybedtools
import multiprocessing as mp
from joblib import Parallel, delayed
from common import parse_config, get_chrom_size, get_clusters

plt.rc('font', family='serif')

""" This script computes fold enrichment of a certain feature within a cluster. """

def get_peakfilename_feature(feature):
    mapping_file = config["CHIPSEQ_DIR"] + 'feature_filenames.txt'
    #mapping = np.loadtxt(mapping_file, dtype = 'str', delimiter = ' ')
    mapping = pd.read_csv(mapping_file, sep = '\t')
    # using feature name return the appropriate file with all the data peaks
    feature_file = mapping[mapping['feature'] == feature]['filename'].values[0]
    feature_file = config["CHIPSEQ_DIR"] + feature_file
    return feature_file


def get_num_bases_genome():
    chrom_size_file =  config["GENOME_DIR"] + 'chrom_hg19.sizes'
    chrom_df = pd.read_csv(chrom_size_file, header = None, names = ['chr', 'size'], sep = '\t')
    num_bases_genome = chrom_df.iloc[0:len(config["chrs"]),1].sum()
    return num_bases_genome

def transform_clust2bed(clust):
    clust_df = pd.DataFrame(clust, columns = ['chr', 'start'], dtype = 'S12')
    clust_df['chr'] = clust_df['chr'].apply(lambda x: 'chr' + str(int(float(x))) if ((x != 'X') | (x != 'Y')) else 'chr' + x)
    clust_df['stop'] = clust_df['start'].astype('int') + config["HIC_RESOLN"]
    clust_bed = pybedtools.BedTool.from_dataframe(clust_df)
    return clust_bed

def get_bed_intersect_features(features_active, name):
    features_intersect = pybedtools.BedTool(get_peakfilename_feature(features_active[0]))
    for feature in features_active[1:]:
        feature_file = get_peakfilename_feature(feature)
        feature_bed = pybedtools.BedTool(feature_file)      
        features_intersect = features_intersect.intersect(feature_bed)
    # bases overlap feature
    feature_bed_sorted = features_intersect.sort()
    bases_feature = feature_bed_sorted.total_coverage()

    fold_enrichment_list = []
    for clust in clusters:
        fold_enrichment = fold_enrichment_clust(clust, features_intersect, bases_feature, num_bases_genome)
        fold_enrichment_list.append(fold_enrichment)
    df_enrichment[feature] = fold_enrichment_list
    print feature
    print float(np.sum(np.asarray(fold_enrichment_list) > 1.0))/len(fold_enrichment_list)*100
    
    plt.figure()
    plt.hist(fold_enrichment_list, bins = 40)
    plt.xlabel('Fold enrichment of ' + feature + ' in clusters', fontsize = 12)
    plt.ylabel('# clusters', fontsize = 12)
    plt.savefig(ENRICHMENT_NEW_DIR + 'fold_enrichment_perclust_' + name + '.png')
    plt.close()
    df_enrichment.to_csv(ENRICHMENT_NEW_DIR + 'enrichment_' + name + '.csv')


def fold_enrichment_clust(clust, feature_bed, bases_feature, num_bases_genome):
    # compute fold enrichment of feature in particular cluster
    clust_bed = transform_clust2bed(clust)
    # bases in cluster and overlap feature
    cluster_and_feature = feature_bed.coverage(clust_bed).to_dataframe()
    # get columns with number of overlapping bases (-3 index is based on bedtools help menu)
    bases_cluster_and_feature = cluster_and_feature.iloc[:,-3].sum()
    # bases overlap cluster
    bases_cluster = config["HIC_RESOLN"]*len(clust)
    # numerator
    numerator = float(bases_cluster_and_feature)/num_bases_genome
    denominator = (float(bases_feature)/num_bases_genome)*(float(bases_cluster)/num_bases_genome)
    fold_enrichment = numerator/denominator
    return fold_enrichment


def run_one_feature(feature, clusters):
    df_enrichment = pd.DataFrame()
    num_bases_genome = get_num_bases_genome()
    
    # get appropriate file
    feature_file = get_peakfilename_feature(feature)
    feature_bed = pybedtools.BedTool(feature_file)
    # bases overlap feature
    feature_bed_sorted = feature_bed.sort()
    bases_feature = feature_bed_sorted.total_coverage()

    fold_enrichment_list = []
    for clust in clusters:
        fold_enrichment = fold_enrichment_clust(clust, feature_bed, bases_feature, num_bases_genome)
        fold_enrichment_list.append(fold_enrichment)
    df_enrichment[feature] = fold_enrichment_list
    print feature
    print float(np.sum(np.asarray(fold_enrichment_list) > 1.0))/len(fold_enrichment_list)*100
    
    #plt.figure()
    #plt.hist(fold_enrichment_list, bins = 40)
    #plt.xlabel('Fold enrichment of ' + feature + ' in clusters', fontsize = 12)
    #plt.ylabel('# clusters', fontsize = 12)
    #plt.savefig(config["ENRICHMENT_DIR"] + 'fold_enrichment_perclust_' + feature + '.png')
    #plt.close()
    #df_enrichment.to_csv(config["ENRICHMENT_DIR"] + 'enrichment_' + feature + '.csv')
 
    return df_enrichment

def main():
    global config
    config_fn = sys.argv[1]
    config = parse_config(config_fn)

    # get clusters
    filename = config["GRAPH_NEW_DIR"] + config["filename"]
    clusters = get_clusters(filename)
    print len(clusters)
    features = ['H3K36me3', 'POLR2A', 'H3K9me3', 'H3K9ac', 'H3K4me3', 'H3K4me1', 'H3K27me3']

    enrichment_list = Parallel(n_jobs=config["NUM_PROC"])(delayed(run_one_feature)(feature, clusters) for feature in features)

    df = pd.DataFrame()
    for df_feat in enrichment_list:
        df = pd.concat([df, df_feat], axis=1)
    df.to_csv(config["ENRICHMENT_DIR"] + 'enrichment_all.csv')


if __name__ == "__main__":
    main()
