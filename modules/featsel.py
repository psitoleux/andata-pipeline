import scanpy as sc

import anndata2ri
import logging
import numpy as np

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

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
