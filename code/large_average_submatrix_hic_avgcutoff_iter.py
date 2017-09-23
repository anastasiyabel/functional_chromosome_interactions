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
import argparse
import csv
import pandas as pd
import pickle
from collections import defaultdict
import operator
from scipy.sparse import csr_matrix
import itertools
import os.path
import os
from joblib import Parallel, delayed
from common import parse_config, get_chrom_size

def simulate_data():
	data = np.random.randn(20, 40)
	x = 2.5*np.random.randn(5,3) + 10
	data[2:7, 3:6] = x
	return data

def transform(data):
	data_trans = np.log(1 + data)
	return data_trans

def rescale(data, mean, std):
	# convert to z scores
	z = (data-mean)/std
	return z

def check_submatrix_below_threshold(sub_matrix, threshold):
	# check whether score is below threshold based on the the size of submatrix
	num_rows, num_cols = sub_matrix.shape
	avg_sub_matrix = np.sqrt(num_rows*num_cols)*np.average(sub_matrix)
	
	#threshold = scipy.stats.norm.ppf(percentile_threshold)
	print avg_sub_matrix, threshold
	if (avg_sub_matrix > threshold):
		check = True
	else:
		check = False
	return check

def residual(u,data, rows, cols):
	avg = np.mean(u)
	# only subtract avg for submatrix
	data[np.ix_(rows,cols)] = data[np.ix_(rows,cols)] - avg
	return data

def large_average_submatrix_adj(data, chr1, chr2, threshold_new):
	# store some data on iterations
	dir = config["HIC_NEW_DIR"] + str(chr1) + '_' + str(chr2) + '/'
	if (os.path.isdir(dir) == False):
		os.makedirs(dir)

	# algorithm until score falls below threshold
	continue_search = True
	iter = 0
	# store matrices
	start_rows, stop_rows, start_cols, stop_cols, best_score_list, avg_list = [], [], [], [], [], []

	while (continue_search):
		rows, cols, sub_matrix, best_score = search_main(data, dir, iter)
		data = residual(sub_matrix, data, rows, cols)
		# check whether score is below threshold based on the the size of submatrix
		continue_search = check_submatrix_below_threshold(sub_matrix, threshold_new)
		if (continue_search == True):
			start_rows.append(rows[0])
			stop_rows.append(rows[-1])
			start_cols.append(cols[0])
			stop_cols.append(cols[-1])
			best_score_list.append(best_score)
			avg_list.append(np.average(sub_matrix))
			iter = iter + 1
			print 'Best score = ', best_score
			print 'Average = ', np.average(sub_matrix)
			print rows[0], rows[-1], cols[0], cols[-1]

	return start_rows, stop_rows, start_cols, stop_cols, best_score_list, avg_list


def search_main(data, dir, iter):
	num_iter = 100
	# keep track of submatrix params that you get out
	search_attributes = np.empty((num_iter, 5))
	for iteration in range(num_iter):
		start_row, k, start_col, l, curr_score = search(data)
		search_attributes[iteration] = start_row, k, start_col,l ,curr_score

	# save the iterations
	np.savetxt(dir + 'sub_matrix' + str(iter) + '.txt', search_attributes)


	best_start_row, best_k, best_start_col, best_l, best_score =  search_attributes[np.argmax(search_attributes[:,4])]
	rows = np.arange(best_start_row, best_start_row + best_k, dtype = 'int')
	cols = np.arange(best_start_col, best_start_col + best_l, dtype = 'int')
	sub_matrix = data[np.ix_(rows,cols)]
	return rows, cols, sub_matrix, best_score

def score(u, data):
	
	m,n = data.shape
	k,l = u.shape
	tau = np.mean(u)
	sc = comb(m,k)*comb(n,l)*norm.cdf(-tau*np.sqrt(k*l))
	sc = -np.log(sc)
	return sc

def score_sum(sum_u, k,l, data):
	m,n = data.shape
	cnr = gammaln(m + 1) - gammaln(k+1) - gammaln(m-k+1)
	cnc = gammaln(n + 1) - gammaln(l+1) - gammaln(n-l+1)
	ar = sum_u/np.sqrt(k*l)
	rest2 = -(ar*ar)/2.0 + np.log(erfcx(ar/np.sqrt(2))*0.5)
	sc = -rest2 -cnr - cnc
	return sc

def grouped_sum(array, N):
	length = len(array) - N + 1
	# initialize the array
	adj_sum = np.zeros((length))
	for i in range(0,N):
		adj_sum = adj_sum + array[i:length+i]
	return adj_sum


def search(data):
	# run search procedure with fixed k, l first
	max_num_rows = int(10000000.0/config["HIC_RESOLN"])
	max_num_cols = int(10000000.0/config["HIC_RESOLN"])
	k = random.randint(1, max_num_rows)
	l = random.randint(1, max_num_cols)
	row_set, col_set = search_fixed_k_l(data, k,l)

	# allow k and l to vary
	# initialize the running average
	pre_score = -1000000
	curr_score = 0
	# iterate until convergence
	while(pre_score != curr_score):
		# sum across columns
		row_summed = np.sum(col_set, axis =1)
		start_row, k, score_rows = enumerate_adj_submatrix_scores(data, row_summed, max_num_rows, k, l, 'row')
		# make a row set
		row_set = data[start_row:start_row+k, :]

		# columns
		col_summed = np.sum(row_set, axis =0)
		start_col, l, score_cols = enumerate_adj_submatrix_scores(data, col_summed, max_num_cols, k, l, 'col')
		# make a col set
		col_set = data[:,start_col:start_col+l]

		# update scores
		pre_score = curr_score
		curr_score = score_cols
		#print 'Score = ', pre_score, curr_score

	return start_row, k, start_col, l, curr_score

def enumerate_adj_submatrix_scores(data, row_summed, max_num_rows, k, l, row_or_col):
	if (row_or_col == 'row'):
			start_row_best_list = []
			start_row_best_ind_list = []
			# let the number of rows to include vary (+1 to make the range inclusive)
			possible_num_rows = range(1, max_num_rows + 1)
			for i in possible_num_rows:
				# make all possible submatrices by summing adjacent rows
				adj_row_sum = grouped_sum(row_summed, i)
				score_list = [score_sum(sum_u, i, l, data) for sum_u in adj_row_sum]
				# find best starting row
				start_row_best_ind, start_row_best = max(enumerate(score_list), key=operator.itemgetter(1))
				start_row_best_ind_list.append(start_row_best_ind)
				start_row_best_list.append(start_row_best)
	if (row_or_col == 'col'):
			start_row_best_list = []
			start_row_best_ind_list = []
			possible_num_rows = range(1, max_num_rows + 1)
			for i in possible_num_rows:
				# make all possible submatrices by summing adjacent rows
				adj_row_sum = grouped_sum(row_summed, i)
				# LINE BELOW is THE ONLY DIFFENECE BETWEEN ROW AND COL CODE
				score_list = [score_sum(sum_u, k, i, data) for sum_u in adj_row_sum]
				# find best starting row
				start_row_best_ind, start_row_best = max(enumerate(score_list), key=operator.itemgetter(1))
				start_row_best_ind_list.append(start_row_best_ind)
				start_row_best_list.append(start_row_best)
	# choose the best scoring 
	ind, score_rows = max(enumerate(start_row_best_list), key=operator.itemgetter(1))
	start_row = start_row_best_ind_list[ind]
	k = possible_num_rows[ind]
	return start_row, k, score_rows

def search_fixed_k_l(data, k, l):
	# initialize (select l adjacent columns at random)
	num_rows = data.shape[0]
	num_cols = data.shape[1]
	# choose a random starting position for column
	start_col = random.randint(0, num_cols-l)
	col_set = data[:,start_col:start_col+l]

	# initialize the running average
	pre_avg = -1000000
	curr_avg = 0
	# iterate until convergence
	while(pre_avg != curr_avg):
		# get k rows with the largest sum over l columns (adjacent rows)
		# make another matrix that is the sum of k adjacent columnns
		row_summed_data = np.asarray([np.sum(col_set[i:i+k,:]) for i in range(0, col_set.shape[0]-k+1)])
		# choose starting row that gave the largest sum
		start_row = np.argmax(row_summed_data)
		row_set = data[start_row:start_row+k, :]

		# get l rows with the largest sum over k rows (adjacent columns)
		# make another matrix that is the sum of l adjacent rows
		col_summed_data = np.asarray([np.sum(row_set[:,j:j+l]) for j in range(0, row_set.shape[1]-l+1)])
		# choose starting row that gave the largest sum
		start_col = np.argmax(col_summed_data)
		col_set = data[:,start_col:start_col+l]

		# compute the new average of the submatrix
		sub_matrix = data[np.ix_(range(start_row, start_row+k), range(start_col, start_col+l))]
		# update averages
		pre_avg = curr_avg
		curr_avg = np.mean(sub_matrix)
		#print curr_avg, pre_avg
	return row_set, col_set

def df_remove_zeros_rows_cols(df):
	# remove indices of matrix that are all zero col or row wise
	# drop rows that are all 0
	df = df[(df.T != 0).any()]
	# drop columns that are all 0
	df = df[df.columns[(df != 0).any()]]
	return df

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

def map_pos2rownum(df, row_pos):
	return np.where(df.index.values == row_pos)[0][0]
def map_pos2colnum(df, row_pos):
	return np.where(df.columns.values == str(int(row_pos)))[0][0]

def map_rownum2pos(df, row_num):
	positions = df.index.values
	return positions[row_num]

def map_colnum2pos(df, col_num):
	positions = df.columns.values
	return float(positions[col_num])

def map_num2pos(df, start_rows, stop_rows, start_cols, stop_cols):
	# for each row and column number figure out what position on the genome does it correspond to
	start_row_pos = [map_rownum2pos(df, row_num) for row_num in start_rows]
	stop_row_pos = [map_rownum2pos(df, row_num) for row_num in stop_rows]
	start_col_pos = [map_colnum2pos(df, col_num) for col_num in start_cols]
	stop_col_pos = [map_colnum2pos(df, col_num) for col_num in stop_cols]
	return start_row_pos, stop_row_pos, start_col_pos, stop_col_pos

def numclust_avg(pair):
	chr1, chr2 = pair
	fname = config["HIC_NEW_DIR"]  + 'intermingling_regions.chr' + str(chr1) + '_chr' + str(chr2) + '.csv'
	# check if the file exists
	if (os.path.isfile(fname) == True):
		

		df_intermingling = pd.read_csv(fname, index_col = 0)
		plt.figure()
		plt.plot(xrange(df_intermingling.shape[0]), df_intermingling['score'], 'o-')
		plt.xlabel('Cluster #')
		plt.ylabel('Score')
		plt.savefig(config["HIC_NEW_DIR"] + 'cluster_score.chr' + str(chr1) + '_chr' + str(chr2) + '.png')
		plt.close()

                plt.figure()
                plt.plot(xrange(df_intermingling.shape[0]), df_intermingling['avg'], 'o-')
                plt.xlabel('Cluster #')
                plt.ylabel('Average')
                plt.savefig(config["HIC_NEW_DIR"] + 'cluster_average.chr' + str(chr1) + '_chr' + str(chr2) + '.png')
                plt.close()

def determine_min_max_hic():
    # get the minimum and maximum transformed Hi-C contact for plotting
    chr_pairs = list(itertools.combinations(config["chrs"], 2))
    min_list = []
    max_list = []
    for pair in chr_pairs:
        chr1, chr2 = pair
        fname = config["HIC_NEW_DIR"]  + 'intermingling_regions.chr' + str(chr1) + '_chr' + str(chr2) + '.csv'
        if (os.path.isfile(fname) == True):
                # read in hic matrix
                hic_filename =  config["HIC_FILT_DIR"] +'chr' + str(chr1) + '_chr' + str(chr2) + '.zscore.txt'
                df = pd.read_csv(hic_filename, index_col = 0)
                data = df.as_matrix()
		min_chr_pair = np.min(data)
		max_chr_pair = np.max(data)
		min_list.append(min_chr_pair)
		max_list.append(max_chr_pair)
    minl = min(min_list)
    maxl = max(max_list)
    return minl, maxl

def draw_identified_LASregions(pair, minl, maxl):
	chr1, chr2 = pair
	plt.rc('font', family='serif')
	# no gridlines
	seaborn.set_style("dark", {'axes.grid':False})

	numclust = 50
	fname = config["HIC_NEW_DIR"]  + 'intermingling_regions.chr' + str(chr1) + '_chr' + str(chr2) + '.csv'
	# check if the file exists
	if (os.path.isfile(fname) == True):

		hic_filename =  config["HIC_FILT_DIR"] +'chr' + str(chr1) + '_chr' + str(chr2) + '.zscore.txt'
		df = pd.read_csv(hic_filename, index_col = 0)
		data = df.as_matrix()

		#plt.figure(figsize = (100, 100))
		plt.figure()
		plt.imshow(data, cmap = 'Reds', vmin = minl, vmax = maxl)
		cbar = plt.colorbar()
		#cbar.set_label('log(1+x) transformed rescaled HiC observed contacts')
		cbar.set_label('Transformed Hi-C contacts', fontsize = 12)
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

		#plt.savefig(config["HIC_NEW_DIR"] + 'hic_transformed_rescaled.chr' + str(chr1) + '_chr' + str(chr2) + '.png')
		plt.savefig(config["HIC_NEW_DIR"] + 'hic_transformed_rescaled.chr' + str(chr1) + '_chr' + str(chr2) + 'commonscale.png')

		df_intermingling = pd.read_csv(fname, index_col = 0)
		# iterate over all las regions found
		for num in range(0, len(df_intermingling)):
			region = df_intermingling.iloc[num]
			start_row = map_pos2rownum(df, region['start row'])
			stop_row = map_pos2rownum(df, region['stop row'])
			start_col = map_pos2colnum(df, region['start col'])
			stop_col = map_pos2colnum(df, region['stop col'])

			# draw vertical lines - columns are same
			plt.plot([start_col, start_col], [start_row, stop_row], 'k-', lw = 0.8)
			plt.plot([stop_col, stop_col], [start_row, stop_row], 'k-', lw = 0.8)
			# draw horizontal lines - rows are same
			plt.plot([start_col, stop_col], [start_row, start_row], 'k-', lw = 0.8)
			plt.plot([start_col, stop_col], [stop_row, stop_row], 'k-', lw = 0.8)

		#plt.savefig(fname.split('.csv')[0] + '.png')
		plt.savefig(fname.split('.csv')[0] + 'commonscale.png')
		#plt.savefig(fname.split('.csv')[0] + '.pdf', format = 'pdf', dpi = 1000)
		plt.close()



def run_LAS(pair, threshold_new):
	chr1, chr2 = pair
	fname = config["HIC_NEW_DIR"] + 'intermingling_regions.chr' + str(chr1) + '_chr' + str(chr2) +'.avg_filt.csv'
	if (os.path.isfile(fname) == False):
		print chr1, chr2
		# read in hic matrix
		hic_filename =  config["HIC_FILT_DIR"] +'chr' + str(chr1) + '_chr' + str(chr2) + '.zscore.txt'
		df = pd.read_csv(hic_filename, index_col = 0)
		data = df.as_matrix()	

		# run LAS algorithm
		start_rows, stop_rows, start_cols, stop_cols, best_score_list, avg_list = large_average_submatrix_adj(data, chr1, chr2, threshold_new)

		# convert indices to positions
		start_row_pos, stop_row_pos, start_col_pos, stop_col_pos = map_num2pos(df, start_rows, stop_rows, start_cols, stop_cols)

		# store results in dataframe
		dic = {'start row' : pd.Series(start_row_pos), 'stop row' : pd.Series(stop_row_pos), 'start col' : pd.Series(start_col_pos), 'stop col': pd.Series(stop_col_pos), 'score' : pd.Series(best_score_list), 'avg' : pd.Series(avg_list)}
		df_intermingling = pd.DataFrame(dic, columns=dic.keys())
		df_intermingling.to_csv(fname)		

def main():
	global config
	config_fn = sys.argv[1]
	config = parse_config(config_fn)

	#percentile_new =  0.999999999999999 #7.941

	iters = 100
	threshold_new = scipy.stats.norm.ppf(config["pvalue_threshold"])
	print "LAS z-score threshold = ", threshold_new

	chr_pairs = list(itertools.combinations(config["chrs"], 2))
	Parallel(n_jobs = config['NUM_PROC'])(delayed(run_LAS)(pair, threshold_new) for pair in chr_pairs)  

	#minl, maxl = determine_min_max_hic()
	#Parallel(n_jobs = config['NUM_PROC'])(delayed(draw_identified_LASregions)(pair, minl, maxl) for pair in chr_pairs)  
	#Parallel(n_jobs = config['NUM_PROC'])(delayed(numclust_avg)(pair) for pair in chr_pairs)
	

if __name__ == "__main__":
    main()
