from abmodule import Module
from scanpy import sc

import scanpy as sc
import anndata2ri
import logging
import numpy as np

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

rcb.logger.setLevel(logging.ERROR)


def soupx(adata: sc.AnnData, raw: sc.AnnData, min_cells: int = 20) -> sc.AnnData:
    """
    Run SoupX on the input AnnData object.

    Args:
        adata (AnnData): Input AnnData object with normalized counts. 
        raw (AnnData): Input AnnData object containing raw counts.
        min_cells (int): Minimum number of cells a gene must be detected in. Default is 20.

    Returns:
        AnnData: AnnData object with corrected counts stored in the 'soupX_counts' layer.
    """
    
    ro.pandas2ri.activate()
    anndata2ri.activate()

    # Extract the transposed droplet table from the raw AnnData object
    data_tod = raw.X.T
    del raw

    # Import R library
    soupx = importr("SoupX")
    
    adata_pp = adata.copy()

    # Compute leiden clusters
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")

    # Extract the leiden clusters
    soupx_groups = adata_pp.obs["soupx_groups"]

    # Free memory by deleting the processed copy
    del adata_pp

    # Prepare the input data for SoupX
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T
    
    # Run SoupX
    # Set up R environment
    ro.globalenv["data"] = data
    ro.globalenv["data_tod"] = data_tod
    ro.globalenv["genes"] = genes
    ro.globalenv["cells"] = cells
    ro.globalenv["soupx_groups"] = soupx_groups

    # Run SoupX in R
    ro.r('''
        library(SoupX)
        # Specify row and column names of data
        rownames(data) <- genes
        colnames(data) <- cells
        # Ensure correct sparse format for table of counts and table of droplets
        data <- as(data, "sparseMatrix")
        data_tod <- as(data_tod, "sparseMatrix")
        
        # Generate SoupChannel Object for SoupX 
        sc <- SoupChannel(data_tod, data, calcSoupProfile = FALSE)
        
        # Add extra meta data to the SoupChannel object
        soupProf <- data.frame(row.names = rownames(data), 
                               est = rowSums(data) / sum(data), 
                               counts = rowSums(data))
        sc <- setSoupProfile(sc, soupProf)
        # Set cluster information in SoupChannel
        sc <- setClusters(sc, soupx_groups)
        
        # Estimate contamination fraction
        sc <- autoEstCont(sc, doPlot = FALSE)
        # Infer corrected table of counts and round to integer
        out <- adjustCounts(sc, roundToInt = TRUE)
    ''')
    
    # Retrieve the output
    out = np.array(ro.globalenv["out"])

    # Store the corrected counts in the AnnData object
    adata.layers["counts"] = adata.X
    adata.layers["soupX_counts"] = out.T
    adata.X = adata.layers["soupX_counts"]

    # Filter out genes not detected in at least min_cells cells
    sc.pp.filter_genes(adata, min_cells=min_cells)

    print(f"Total number of genes after filtering: {adata.n_vars}")

    return adata
