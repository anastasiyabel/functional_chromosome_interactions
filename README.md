# Network analysis for identification of functional interchromosomal interactions

This pipeline identifies clusters of genomic regions that are both spatially colocalized in 3D and coregulated. We present an intergrative approach that leverages 1D functional genomic features (e.g. epigenetic marks) with 3D interactions from high-throughput chromosome conformation capture (Hi-C) data to identify functional interchromosomal interactions. In the paper "Network analysis identifies chromosome intermingling regions as regulatory hotspots for transcription" (to appear), we:

* identify highly interacting domains in interchromosomal Hi-C maps by determining large average submatrices
* superimpose regulatory marks on the interacting domains
* construct network of interacting regions with edges weighted by the correlation of the superimposed marks
* cluster network to obtain spatially colocalized and coregulated domains


In order to identify interacting domains, [filter](https://github.com/anastasiyabel/functional_chromosome_interactions/blob/master/code/filter_hic_contacts.py) interchromosomal [Hi-C matrices](https://github.com/anastasiyabel/functional_chromosome_interactions/tree/master/inter_chromosome_250kb_imr90) (250kb resolution) and run [large average submatirx (LAS) finding algorithm](https://github.com/anastasiyabel/functional_chromosome_interactions/blob/master/code/large_average_submatrix_hic_avgcutoff_iter.py). These will be the nodes of the network.
```python
python filter_hic_contacts.py ../run_params.json
python large_average_submatrix_hic_avgcutoff_iter.py ../run_params.json
```

To obtain genomic features, count the number of peaks for each feature (at 250kb resolution). The files for each feature must be provided in the [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html) and specified [here](https://github.com/anastasiyabel/functional_chromosome_interactions/blob/master/peaks/feature_filenames.txt).

```python
python get_feature_matrix_perchr.py ../run_params.json
python normalize_epigenetic_tracks.py ../run_params.json
```
In order to quantify coregulation between all regions that were found to be interacting via LAS, compute the correlation between feature profiles for each interacting pair of regions - this will correspond to edge weights in the network.

```python
python get_edge_weights.py ../run_params.json
```
Construct a network of interacting regions, weighted by correlation between each pair of regions and cluster the network.
```python
python graph_weighted.py ../run_params.json
```
Perform analysis on the clusters of the network, for example, find fold enrichment of active and repressive regulatory marks in each cluster.
```python
python feature_fold_enrichment.py ../run_params.json
```
The directories for input/output, resolution of Hi-C data, and specification of threshold for LAS algorithm are located in the configuration file [run_params.json](https://github.com/anastasiyabel/functional_chromosome_interactions/blob/master/run_params.json).

Additional programs are required to run code in this repository:
* bedtools
* pybedtools
* pandas
* numpy
* matplotlib
* pickle
* itertools
* scipy
* joblib
* sys
* collections
* networkx
* seaborn
* subprocess
* installation of [weighted correlation clustering](http://www.ling.ohio-state.edu/~elsner.14/resources/correlation-readme.html)
