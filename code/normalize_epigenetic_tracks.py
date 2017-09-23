import sys
import numpy as np 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import seaborn
import pandas as pd
import pickle
from common import parse_config

def get_filtered_chipseq(chr, blacklist):
    df_chipseq = pd.read_csv(config["CHIPSEQ_OUT_DIR"] + 'features_matrix_chr' + str(chr) + '.csv', index_col = 0)
    # get all blacklisted loccations
    blacklist_chr = blacklist[chr]
    # get a list of columns to keep
    allcols = set(map(int,df_chipseq.columns))
    cols2keep = allcols - blacklist_chr
    df_chipseq_filt = df_chipseq[list(map(str,cols2keep))]
    return df_chipseq_filt

def get_mean_std():
    blacklist = pickle.load(open(config["HIC_FILT_DIR"] + 'blacklist.pickle', 'rb'))
    # collect chipseq data across all chromosomes into one dataframe
    df_all = pd.DataFrame()
    chr_list = config["chrs"]
    for chr in chr_list:
        df_chipseq_filt = get_filtered_chipseq(chr, blacklist)
        df_all = pd.concat([df_all, df_chipseq_filt],axis=1)
    # transform
    df_all = np.log(df_all + 1)
    # find mean and standard dev
    mean_features = np.mean(df_all, axis =1)
    std_features = np.std(df_all, axis=1)
    return mean_features, std_features

def normalize_chipseq(mean_features, std_features):
    chr_list = config["chrs"]
    for chr in chr_list:
        # get chipseq data
        df_chipseq = pd.read_csv(config["CHIPSEQ_OUT_DIR"] + 'features_matrix_chr' + str(chr) + '.csv', index_col = 0)
        # transform
        df_chipseq = np.log(df_chipseq + 1)
        # normalize
        df_norm = (df_chipseq.T - mean_features)/std_features
        # transpose back
        df_norm = df_norm.T
        # save
        df_norm.to_csv(config["CHIPSEQ_OUT_DIR"] + 'features_matrix_chr' + str(chr) + 'norm.csv')



def main():
    global config
    config_fn = sys.argv[1]
    config = parse_config(config_fn)

    mean_features, std_features = get_mean_std()
    normalize_chipseq(mean_features, std_features)

if __name__ == "__main__":
    main()