import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import random
from scipy import sparse
from scipy.special import comb
from scipy.special import gammaln
from scipy.special import erfcx
from scipy.stats import norm
import scipy.stats
import seaborn
import csv
import pandas as pd
import pickle
from collections import defaultdict
import operator
from scipy.sparse import csr_matrix
import itertools
import os.path
import math
import pybedtools
from common import parse_config, get_chrom_size
from large_average_submatrix_hic_avgcutoff_iter import map_rownum2pos, map_colnum2pos

""" This script filters out centromeric regions, row and column sum outliers, and repeat regions from the HiC contact matrix.
Each matrix is corrected by mean and standard devation based on the mean and std of HiC across the entire genome."""

def get_hic_matrix(hic_filename, chr1, chr2):
        # construct matrix where each axis corresponds to position along the chromosome
        # this matrix will be rectangular b/c sizes of interchromosomal matrices are not the same
        # returns scipy sparse matrix
        data = np.loadtxt(hic_filename, delimiter = '\t')
        row_ind = data[:,0]/config["HIC_RESOLN"]
        col_ind = data[:,1]/config["HIC_RESOLN"]
        contact_values = data[:,2]
        # obtain chromosome sizes
        chr1_size = get_chrom_size(chr1, config["GENOME_DIR"])/config["HIC_RESOLN"] + 1
        chr2_size = get_chrom_size(chr2, config["GENOME_DIR"])/config["HIC_RESOLN"] + 1

        hic_matrix = csr_matrix((contact_values, (row_ind, col_ind)), shape = (chr1_size, chr2_size))
        hic_dense = np.asarray(hic_matrix.todense())
        row_labels = np.arange(chr1_size)*config["HIC_RESOLN"]
        col_labels = np.arange(chr2_size)*config["HIC_RESOLN"]
        df = pd.DataFrame(hic_dense, index = row_labels, columns = col_labels)
        # get rid of nans
        df = df.fillna(0)
        return df

def plot_hic_matrix(df, chr1, chr2, label, plotname):
    data = df.as_matrix()
    plt.figure()
    plt.imshow(data, cmap = 'Reds')
    cbar = plt.colorbar()
    #cbar.set_label('log(1+x) transformed rescaled HiC observed contacts')
    #cbar.set_label('Transformed HiC contacts', fontsize = 12)
    cbar.set_label(label)
    cbar.solids.set_rasterized(True)
    
    # label ticks with genomic position (Mb)
    xaxis = range(0, df.shape[1], 100)
    xlabels = [str(map_colnum2pos(df, x)/1000000.0) for x in xaxis]
    plt.xticks(xaxis, xlabels)
    yaxis = range(0, df.shape[0], 100)
    ylabels = [str(map_rownum2pos(df, y)/1000000.0) for y in yaxis]
    plt.yticks(yaxis, ylabels)

    plt.xlabel('chr' + str(chr2) + ' (Mb)', fontsize = 14)
    plt.ylabel('chr' + str(chr1) + ' (Mb)', fontsize = 14)

    #plt.savefig(config["HIC_DIR"] + 'hic_transformed_rescaled.chr' + str(chr1) + '_chr' + str(chr2) + '.png')
    plt.savefig(plotname)
    plt.close()

def transform(data):
    data_trans = np.log(1 + data)
    return data_trans

def filter_centromere(df, chrom, row_or_col, filter_size = 1000000):
    # get centromere locations
    centrom_filename = config["GENOME_DIR"] + 'chrom_hg19.centromeres'
    df_centrom = pd.read_csv(centrom_filename, sep = '\t')
    chr_data = df_centrom[df_centrom['chrom'] == 'chr' + str(chrom)]
    centrom_start = int(math.floor(float(chr_data['chromStart'])/config["HIC_RESOLN"])*config["HIC_RESOLN"])
    centrom_end = int(math.ceil(float(chr_data['chromEnd'])/config["HIC_RESOLN"])*config["HIC_RESOLN"])

    centrom_start = centrom_start - filter_size
    centrom_end = centrom_end + filter_size

    if (row_or_col == 'row'):
        df.loc[centrom_start:centrom_end, :] = 0
    if (row_or_col == 'col'):
        df.loc[:, centrom_start:centrom_end] = 0
    return df

def find_repeat_locations_wg():
    # get repeats regions
    filename = config["GENOME_DIR"] +  'rmsk.txt'
    df_repeats = pd.read_csv(filename, header=None , sep= '\t', usecols = (5,6,7), names = ['chr', 'start', 'stop'])
    df_repeats_bed = pybedtools.BedTool.from_dataframe(df_repeats)

    # make a dictionary of all repeat regions for each chromosome - these regions should be taken out from hic
    # this threshold is based on the 95th percentile of repeat coverage histgram over 250kb regions - PLEASE CHANGE THIS IF YOU WANT
    threshold = 0.63009760000000004
    dic_repeats_tofilter = {}
    df_coverage = pd.DataFrame()
    #chr_list = range(1,23)
    chr_list = config["chrs"]
    for chr in chr_list:
        print chr
        # construct a dataframe of all possible start sites for this chr
        chr_size = get_chrom_size(chr, config["GENOME_DIR"])
        start_regions = range(0, chr_size, config["HIC_RESOLN"])
        df_chr = pd.DataFrame({'chr':['chr' + str(chr)]*len(start_regions), 'start':start_regions})
        df_chr['stop'] = df_chr['start'] + config["HIC_RESOLN"]
        df_chr_bed = pybedtools.BedTool.from_dataframe(df_chr)
        coverage_chr = df_chr_bed.coverage(df_repeats_bed).to_dataframe(names = ['chr', 'start', 'end', 'bases covered', 'length A', 'length B', 'coverage'])
        tofilter = coverage_chr[coverage_chr['coverage'] >= threshold]['start'].values
        dic_repeats_tofilter[chr] = tofilter
        df_coverage = pd.concat([df_coverage, coverage_chr])
    return dic_repeats_tofilter

def filter_repeats(df, chr, dic_repeats_tofilter, row_or_col):
    regions2filter = dic_repeats_tofilter[chr]
    if (row_or_col == 'row'):
        df.loc[regions2filter,:] = 0
    if (row_or_col == 'col'):
        df.loc[:,regions2filter] = 0
    return df

def filter_hic(df, row, col, threshold_row, threshold_col):
    # take hic matrix and zero out columns or rows that are outliers
    ind_row = np.arange(len(row))[row>threshold_row]
    df.iloc[ind_row,:] = 0
    ind_col = np.arange(len(col))[col>threshold_col]
    df.iloc[:,ind_col] = 0
    return df

def raw_matrices():
    chr_pairs = list(itertools.combinations(config["chrs"], 2))
    for pair in chr_pairs:
            print pair
            chr1, chr2 = pair
            hic_filename =  config["HIC_DIR"] +'chr' + str(chr1) + '_chr' + str(chr2) + '.txt'
            df = get_hic_matrix(hic_filename, chr1, chr2)
            # make the minimum value 0
    df = df - np.min(np.min(df))
    # transdorm
    df = transform(df)
    # plot the matrix
    plot_hic_matrix(df, chr1, chr2, 'log-transformed Obs/Expected', config["HIC_DIR"] + 'chr' + str(chr1) + '_' + 'chr' + str(chr2) + '.png')   


def row_col_sums(dic_repeats_tofilter):
    # record nonzero values over all hic matrices
    nonzero_entries = []

    #chr_pairs = list(itertools.combinations(range(1, 23), 2))
    chr_pairs = list(itertools.combinations(config["chrs"], 2))
    #chr_pairs = list(itertools.product(range(1,23), ['X']))
    for pair in chr_pairs:
        print pair
        chr1, chr2 = pair
        hic_filename =  config["HIC_DIR"] +'chr' + str(chr1) + '_chr' + str(chr2) + '.txt'
        df = get_hic_matrix(hic_filename, chr1, chr2)
        df = df - np.min(np.min(df))
        # plot the matrix
        print 'Before filtering'
        print np.count_nonzero(df.sum(axis=0))
        print np.count_nonzero(df.sum(axis=1))

        # FILTER OUT CENTROMERIC REGIONS (2Mb)
        df = filter_centromere(df, chr1, 'row', filter_size = 2000000)
        df = filter_centromere(df, chr2, 'col', filter_size = 2000000)
        print 'After centromere filtering'
        print np.count_nonzero(df.sum(axis=0))
        print np.count_nonzero(df.sum(axis=1))
        # filter out repeats
        df = filter_repeats(df, chr1, dic_repeats_tofilter, 'row')
        df = filter_repeats(df, chr2, dic_repeats_tofilter, 'col')
        print 'After repeat filtering'
        print np.count_nonzero(df.sum(axis=0))
        print np.count_nonzero(df.sum(axis=1))
        # LOG TRANSFORM
        df_transformed = transform(df)
        
        # get all the row sums
        row_orig = np.sum(df_transformed, axis = 1).as_matrix()
        col_orig = np.sum(df_transformed, axis = 0).as_matrix()
        # get rid of nonzero
        row = row_orig[np.nonzero(row_orig)]
        col = col_orig[np.nonzero(col_orig)]

        threshold_row = detect_upper_outliers(row)
        threshold_col = detect_upper_outliers(col)

        # plt.figure()
        # seaborn.violinplot(row, orient = 'v')
        # plt.ylabel('HiC contacts log(1+x) transformed row sums')
        # plt.title('chr' + str(chr1) + ' - ' + str(chr2))
        # plt.savefig(config["HIC_FILT_DIR"] + 'row_sums.chr' + str(chr1) + '_chr' + str(chr2) + '.png')
        # plt.close()

        # plt.figure()
        # seaborn.violinplot(col, orient = 'v')
        # plt.ylabel('HiC contacts log(1+x) transformed row sums')
        # plt.title('chr' + str(chr1) + ' - ' + str(chr2))
        # plt.savefig(config["HIC_FILT_DIR"] + 'col_sums.chr' + str(chr1) + '_chr' + str(chr2) + '.png')
        # plt.close()

        print np.sum(row>threshold_row)
        print np.sum(col>threshold_col)

        df_filt = filter_hic(df_transformed, row_orig, col_orig, threshold_row, threshold_col)
        # save new dataframe
        df_filt.to_csv(config["HIC_FILT_DIR"] + 'chr' + str(chr1) + '_chr' + str(chr2) + '.txt')

        # record nonzero entries
        data = df_filt.as_matrix()
        data_nonzero = data[np.nonzero(data)]
        nonzero_entries.append(data_nonzero)

    # save the nonzero entries
    nonzero_entries = np.asarray(list(itertools.chain.from_iterable(nonzero_entries)))
    np.savetxt(config["HIC_FILT_DIR"] + 'whole_genome_nonzero.logtrans.txt', nonzero_entries)
    return nonzero_entries

def whole_genome_mean_std(nonzero_entries):
    mean = np.mean(nonzero_entries)
    std = np.std(nonzero_entries)
    print 'Mean = ', mean
    print 'St. dev. = ', std
    return mean, std

def detect_upper_outliers(arr):
    p25 = np.percentile(arr, 25)
    p75 = np.percentile(arr, 75)
    upper = p75 + 1.5*(p75-p25)
    return upper

def z_score_hic_matrix(mean, std):
    chr_pairs = list(itertools.combinations(config["chrs"], 2))
    #chr_pairs = list(itertools.product(range(1,23), ['X']))
    for pair in chr_pairs:
        print pair
        chr1, chr2 = pair
        hic_filename =  config["HIC_FILT_DIR"] +'chr' + str(chr1) + '_chr' + str(chr2) + '.txt'
        df = pd.read_csv(hic_filename, index_col = 0)

        # zscore matrix
        df  = (df - mean)/std
        # save new matrix
        df.to_csv(config["HIC_FILT_DIR"] + 'chr' + str(chr1) + '_chr' + str(chr2) + '.zscore.txt')


def output_blacklist_locations(mean, std):
    # output a BED file with locations that have been removed (these locations will have an observed value of 0)
    print "Generating a blacklist of locations..."
    chr_pairs = list(itertools.combinations(config["chrs"], 2))
    blacklist = defaultdict(list)

    for pair in chr_pairs:
        print pair
        chr1, chr2 = pair
        # read in
        filename = config["HIC_FILT_DIR"] + 'chr' + str(chr1) + '_chr' + str(chr2) + '.txt'
        df = pd.read_csv(filename, index_col=0)
        # find out which columns have all zeroes in them
        zero_cols = df.columns[(df == 0).all()]
        blacklist[chr2].append(zero_cols)
        # flip the dataframe
        df = df.T
        zero_rows = df.columns[(df == 0).all()]
        blacklist[chr1].append(zero_rows)
    # process the list
    for chr in blacklist.keys():
        values_list = blacklist[chr]
        blacklist[chr] = set(map(int, list(itertools.chain.from_iterable(values_list))))
    # pickle this dictionary
    f = open(config["HIC_FILT_DIR"] + 'blacklist.pickle', 'wb')
    pickle.dump(blacklist, f)


def main():
    global config
    config_fn = sys.argv[1]
    config = parse_config(config_fn)
        
    #raw_matrices()
    # construct a dictionary of regions that have more than >50% repeat coverage - these should be filtered out in HiC matrices
    dic_repeats_tofilter = find_repeat_locations_wg()

    nonzero_entries = row_col_sums(dic_repeats_tofilter)
    ##nonzero_entries = np.loadtxt(config["HIC_FILT_DIR"] + 'whole_genome_nonzero.logtrans.txt')
    mean, std = whole_genome_mean_std(nonzero_entries)
    z_score_hic_matrix(mean, std)
    output_blacklist_locations(mean, std)
    
if __name__ == "__main__":
    main()
