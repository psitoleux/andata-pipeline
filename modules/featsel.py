import scanpy as sc

import anndata2ri
import logging
import numpy as np

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects.packages import importr


def get_topk_indices(X : np.array, k : int) -> np.array:
    
    
    return np.argsort(X)[-k:][::-1]

def set_var_key(adata : sc.AnnData, key : str, value : np.array) -> sc.AnnData:
    
    adata.var[key] = value
    
    return adata

def get_boolean_array_from_idx(k : int, idx : np.array) -> np.array:
    
    bool_arr = np.zeros(k, dtype = bool)
    
    bool_arr[idx] = True
    
    return bool_arr
    

def highly_expressed_genes(adata : sc.AnnData, n_top_genes : int = 1800) -> sc.AnnData:
    
    n_counts_per_gene = adata.X.sum(0).A1
    n_counts_sorted = np.sort(n_counts_per_gene)
    
    cutoff = n_counts_sorted[-n_top_genes]
    
    sc.pp.filter_genes(adata,  min_counts = cutoff)

    print('highly expressed genes selected')
    
    return adata

def most_often_expressed_genes(adata : sc.AnnData, n_top_genes : int = 1800) -> sc.AnnData:
    
    nb_cells_nonzero_per_gene = (adata.X > 0).sum(0).A1
    nb_cells_nonzero_sorted = np.sort(nb_cells_nonzero_per_gene)
    
    cutoff = nb_cells_nonzero_sorted[-n_top_genes]
    
    sc.pp.filter_genes(adata,  min_cells = cutoff)
    
    print('most often expressed genes selected')
    
    return adata

def highly_variable_genes(adata : sc.AnnData, n_top_genes : int = 1800) -> sc.AnnData:
    sc.pp.highly_variable_genes(adata, inplace = True) 
    
    return adata 

def deviance(adata : sc.AnnData, raw : sc.AnnData, n_genes : int = 4_000, ) -> sc.AnnData:
    
    ro.pandas2ri.activate()
    anndata2ri.activate()

    # Import R library
    scry = importr("scry")

    # Save the AnnData object in the R environment
    ro.globalenv["adata"] = adata

    # Perform feature selection with deviance on the non-normalized counts matrix
    ro.r('sce = devianceFeatureSelection(adata, assay="X")')

    # Export the binomial deviance values as a vector
    binomial_deviance = ro.r("rowData(sce)$binomial_deviance")
    binomial_deviance = np.array(binomial_deviance)

    # Sort the vector and select the top n_comps highly deviant genes
    idx = binomial_deviance.argsort()[-n_genes:]
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True

    # Save results in the AnnData object
    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = binomial_deviance

    # Compute the mean and dispersion for each gene across all cells
    sc.pp.highly_variable_genes(adata, layer="scran_normalization")
    
    return adata
