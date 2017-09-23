import sys
import argparse
import numpy as np 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import seaborn
import pandas as pd
import pybedtools
from joblib import Parallel, delayed
from common import parse_config, get_chrom_size


def count_peaks(bed, bed_chrom):
    out = pybedtools.bedtool.BedTool.map(bed_chrom, bed, c = 3, o = 'count_distinct')
    # just get the counts
    counts = out.to_dataframe()['name'].values
    return counts

def read_peakfilename(filename):
    bed = pybedtools.BedTool(filename)
    # sort
    return bed.sort()

def make_chrom_bed(chrom):
    # divide the chromosome into segments of HIC_RESOLN length
    chrom_size = get_chrom_size(chrom, config["GENOME_DIR"])
    # include the very last bit of the incomplete 250kb region
    stop_pos = np.arange(config["HIC_RESOLN"], chrom_size + config["HIC_RESOLN"], config["HIC_RESOLN"], dtype = 'int')
    df_chrom = pd.DataFrame()
    df_chrom['chrom'] = ['chr' + str(chrom)]*len(stop_pos)
    df_chrom['start'] = stop_pos - config["HIC_RESOLN"]
    df_chrom['stop'] = stop_pos
    # convert to bed file
    bed_chrom = pybedtools.BedTool.from_dataframe(df_chrom)
    return bed_chrom


def get_matrix_chr(chrom, df):
    bed_chrom = make_chrom_bed(chrom)
    bed_chrom_df = bed_chrom.to_dataframe()
    # make a dataframe to store results into
    feature_matrix = pd.DataFrame(index = df['feature'].values, columns = bed_chrom_df['start'].values)
    for i in range(len(df)):
        f = df.loc[i, 'filename']
        feature = df.loc[i, 'feature']
        # get bed file of the feature
        bed = read_peakfilename(config["CHIPSEQ_DIR"] + f)
        # get counts for this feature and this chromosome
        counts = count_peaks(bed, bed_chrom)
        # store results into matrix
        feature_matrix.loc[feature, :] = counts
    # write feature matrix to file
    feature_matrix.to_csv(config['CHIPSEQ_OUT_DIR'] + 'features_matrix_chr' + str(chrom) + '.csv')


def main():
    global config
    config_fn = sys.argv[1]
    config = parse_config(config_fn)

    # read in file with all names and file names
    df = pd.read_csv(config["CHIPSEQ_DIR"] + 'feature_filenames.txt', sep = '\t')

    chr_list = config["chrs"]
    Parallel(n_jobs = config['NUM_PROC'])(delayed(get_matrix_chr)(chrom, df) for chrom in chr_list)


if __name__ == "__main__":
    main()
