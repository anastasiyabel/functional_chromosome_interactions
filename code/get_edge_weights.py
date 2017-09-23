import matplotlib.pyplot as plt
import seaborn
import numpy as np
import itertools
import pickle 
import pandas as pd
from collections import defaultdict
import argparse
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
import os
import sys
from common import parse_config, get_chrom_size
# make the fonts pretty
plt.rc('font', family='serif')



def bin_intermingling_regions_hic_resoln(df_intermingling, chr1, chr2):
    # break up large intermingling regions into 250kb bins (or resolution of HiC data)
    combos_list = []
    for row in df_intermingling.iterrows():
        row = row[1]
        start_row = range(int(row['start row']), int(row['stop row']) + config["HIC_RESOLN"], config["HIC_RESOLN"])
        start_col = range(int(row['start col']), int(row['stop col']) + config["HIC_RESOLN"], config["HIC_RESOLN"])
        combos = list(itertools.product(start_row, start_col))
        combos_list.append(combos)

    combos_list = np.asarray(list(itertools.chain.from_iterable(combos_list)))

    df_intermingling_binned = pd.DataFrame(combos_list, columns=['start row', 'start col'])
    df_intermingling_binned['chr1'] = [chr1]*len(combos_list)
    df_intermingling_binned['chr2'] = [chr2]*len(combos_list)
    df_intermingling_binned['stop row'] = df_intermingling_binned['start row'] + config["HIC_RESOLN"]
    df_intermingling_binned['stop col'] = df_intermingling_binned['start col'] + config["HIC_RESOLN"]

    df_intermingling_binned.to_csv(config["HIC_NEW_DIR"]  + 'intermingling_regions.chr' + str(chr1) + '_chr' + str(chr2) + '.binned.csv') 
    return df_intermingling_binned



def get_edge_weights_intermingling(method):
    df_inter = pd.DataFrame()
    chr_pairs = list(itertools.combinations(config["chrs"], 2))
    for pair in chr_pairs:
        print pair
        chr1, chr2 = pair
        fname = config["HIC_NEW_DIR"] + 'intermingling_regions.chr' + str(chr1) + '_chr' + str(chr2) + '.avg_filt.csv'
        if (os.path.isfile(fname) == True):
            df_intermingling = pd.read_csv(fname, index_col = 0)
            # check if dataframe is empty
            if (len(df_intermingling) > 0):

                df_intermingling_binned = bin_intermingling_regions_hic_resoln(df_intermingling, chr1, chr2)

                # get correlation values for each intermingling pair
                df_intermingling_binned_corr = get_correlation(df_intermingling_binned, chr1, chr2)

                # concatenate to previous edge weights
                df_inter = pd.concat([df_inter, df_intermingling_binned_corr])
        else:
            print 'Missing: ', pair
    # add intra-chromosomal interactions for each chr
    df_intra = add_intrachromosomal_edges(df_inter, method = method)
    # get rid off any duplicates
    df_intra = df_intra.drop_duplicates()
    print df_intra.shape

    # concatenate inter- and intra-chromosomal regions
    df_edges = pd.concat([df_inter, df_intra])
    return df_edges
    

def get_correlation(df_intermingling_binned, chr1, chr2):
    df_chipseq1 = pd.read_csv(config["CHIPSEQ_OUT_DIR"] + 'features_matrix_chr' + str(chr1) + 'norm.csv', index_col = 0)
    df_chipseq2 = pd.read_csv(config["CHIPSEQ_OUT_DIR"] + 'features_matrix_chr' + str(chr2) + 'norm.csv', index_col = 0)
    # get data only for specific locations that are intermingling
    pos_chr1 = np.unique(df_intermingling_binned['start row'].astype('str'))
    pos_chr2 = np.unique(df_intermingling_binned['start col'].astype('str'))

    features1 = df_chipseq1[pos_chr1]
    features2 = df_chipseq2[pos_chr2]

    features = pd.concat([features1.transpose(), features2.transpose()])
    #corr = np.corrcoef(features)
    corr = spearmanr(features1, features2)[0]

    # get correlation for every intermingling region
    corr_list = []
    for row in df_intermingling_binned.iterrows():
        row = row[1]
        # get corresponding row in the correlation matrix
        start_row = str(int(row['start row']))
        start_col = str(int(row['start col']))

        # get inddices of which correlation you should look up
        ind_row = np.where(pos_chr1 == start_row)[0][0]
        ind_col = np.where(pos_chr2 == start_col)[0][0] + len(pos_chr1)

        if (isinstance(corr, np.ndarray)):
            corr_val = corr[ind_row, ind_col]
            corr_list.append(corr_val)
        else:
            corr_list.append(corr)

    df_intermingling_binned['correlation'] = corr_list
    # keep only positively correlated values
    df_intermingling_binned['correlation'][df_intermingling_binned['correlation'] < 0] = 0

    return df_intermingling_binned

def add_intrachromosomal_edges(df_inter, method = 'tad'):
    df_intra = pd.DataFrame(columns = df_inter.columns)
    if (method == 'tad'):
        print 'Arrowhead intrachromsomal edges...'
        filename = config["INTRACHROM_DIR"] + 'GSE63525_IMR90_Arrowhead_domainlist.txt'
        domains = pd.read_csv(filename, sep = '\t')
        chr_list = config["chrs"]
        for chr in chr_list:
            print chr
            # get domains for this chromosome specifically
            domains_chr = domains[(domains['chr1'] == str(chr))]
            # get intermingling regions for this chr
            df_inter_pos1 = df_inter[(df_inter['chr1'] == chr)]['start row']
            df_inter_pos2 = df_inter[(df_inter['chr2'] == chr)]['start col']
            df_inter_pos = np.unique(pd.concat([df_inter_pos1, df_inter_pos2]))


            interactions_chr_list = []
            for domain in domains_chr.iterrows():
                domain = domain[1]
                start_pos = (domain['x1']/config["HIC_RESOLN"])*config["HIC_RESOLN"]
                stop_pos = (domain['x2']/config["HIC_RESOLN"])*config["HIC_RESOLN"]
                pos = range(start_pos, stop_pos+config["HIC_RESOLN"], config["HIC_RESOLN"])
                pos_list = [p for p in pos if p in df_inter_pos]

                # make combinations of all possible interactions in the loop
                interactions = list(itertools.combinations(pos_list, 2))
                if (len(interactions) != 0):
                    interactions_chr_list.append(interactions)

            interactions_chr_list = np.asarray(list(itertools.chain.from_iterable(interactions_chr_list)))
            df_intra_chr = pd.DataFrame(interactions_chr_list, columns=['start row', 'start col'])
            df_intra_chr['chr1'] = [chr]*len(interactions_chr_list)
            df_intra_chr['chr2'] = [chr]*len(interactions_chr_list)
            df_intra_chr['stop row'] = df_intra_chr['start row'] + config["HIC_RESOLN"]
            df_intra_chr['stop col'] = df_intra_chr['start col'] + config["HIC_RESOLN"]

            df_intra_corr = get_correlation(df_intra_chr, chr, chr)

            df_intra = pd.concat([df_intra, df_intra_corr])

    
    if (method == 'las'):
        print 'LAS blocks intrachromosomal edges...'
        chr_pairs = list(itertools.combinations(range(1, 23), 2))
        for pair in chr_pairs:
            print pair
            chr1, chr2 = pair
            # read in LAS domains
            fname = config["HIC_NEW_DIR"] + 'intermingling_regions.chr' + str(chr1) + '_chr' + str(chr2) + '.avg_filt.csv'
            df_intermingling = pd.read_csv(fname, index_col = 0)
            # iterate through LAS blocks that are detected
            for row in df_intermingling.iterrows():
                row = row[1]
                start_row, stop_row = row['start row'], row['stop row']
                interactions = list(itertools.combinations(np.arange(start_row, stop_row + config["HIC_RESOLN"], config["HIC_RESOLN"]),2))
                # make dataframe to get correlations
                interactions = np.asarray(interactions)
                # check that interactions are non-zero
                if (len(interactions) > 0):
                    df_intra_chr = make_df_intra_chr(interactions, chr1)
                    # get correlation
                    df_intra_corr = get_correlation(df_intra_chr, chr1, chr1, CHIPSEQ_DIR)
                    # append to existing list
                    df_intra = pd.concat([df_intra, df_intra_corr])

                start_col = row['start col']
                stop_col = row['stop col']
                interactions = list(itertools.combinations(np.arange(start_col, stop_col + config["HIC_RESOLN"], config["HIC_RESOLN"]),2))
                # make dataframe to get correlations
                interactions = np.asarray(interactions)
                if (len(interactions) > 0):
                    df_intra_chr = make_df_intra_chr(interactions, chr2)
                    # get correlation
                    df_intra_corr = get_correlation(df_intra_chr, chr2, chr2)
                    # append to existing list
                    df_intra = pd.concat([df_intra, df_intra_corr])
    return df_intra

def make_df_intra_chr(interactions, chr):
    df_intra_chr = pd.DataFrame(interactions, columns = ['start row', 'start col'])
    df_intra_chr['chr1'] = [chr]*len(interactions)
    df_intra_chr['chr2'] = [chr]*len(interactions)
    df_intra_chr['stop row'] = df_intra_chr['start row'] + config["HIC_RESOLN"]
    df_intra_chr['stop col'] = df_intra_chr['start col'] + config["HIC_RESOLN"]
    df_intra_chr = df_intra_chr.astype('int')
    return df_intra_chr

def main():
    global config
    config_fn = sys.argv[1]
    config = parse_config(config_fn)

    method_intra = 'tad'
    df_edges = get_edge_weights_intermingling(method_intra)
    df_edges = df_edges[['chr1', 'start row', 'stop row', 'chr2', 'start col', 'stop col', 'correlation']]
    print df_edges.head()
    data_weighted = df_edges.as_matrix()
    np.savetxt(config["GRAPH_NEW_DIR"] + 'data.weighted.txt', data_weighted, fmt = ['%s', '%d', '%d', '%s', '%d',  '%d',  '%f'])
    print data_weighted.shape

if __name__ == "__main__":
    main()