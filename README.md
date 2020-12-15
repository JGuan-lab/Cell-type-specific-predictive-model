# Cell-type-specific-predictive-model
Cell-type-specific predictive model for prioritizing disease-associated genes and gene sets

The codes were used for the analyses of single-nucleus gene expression data of ASD and controls for identifying cell type-specific ASD-associated genes and gene sets. The analyzed single-nucleus gene expression data was deposited at: https://doi.org/10.5281/zenodo.4319554

## About computation time and memory usage

The codes except MultiLabelModeling are run independently for each cell type. 

Take the cell type of L2/3 as an example, we provide the computation time and memory usage. We run the codes on a computer with the specs: CPU @ 2.10GHz 2.10GHz (2 processors), each processor with 8 cores and 16 threads.

For the code of AllFeatureModeling, it takes about 14 minutes and 199 GB of RAM when 12 threads are using.

For the code of RFEModeling, it takes about 8 hours and 128 GB of RAM when 12 threads are using.

For the code of GenesetModeling, it takes about 13 hours and 18 GB of RAM when 12 threads are using.

For the code of MultiLabelModeling, it takes about 7 days and 200 GB of RAM when 1 thread is using.
