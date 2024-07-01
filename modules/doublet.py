import scanpy as sc
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import numpy as np
import scipy.sparse as sp

def scdbl(adata: sc.AnnData, seed: int = 123) -> sc.AnnData:
    # Import R libraries
    scDblFinder = importr("scDblFinder")
    SingleCellExperiment = importr("SingleCellExperiment")
    Matrix = importr("Matrix")
    
    # Prepare the data matrix
    data_mat = adata.X.T  # Transpose the matrix
    
    if sp.issparse(data_mat):
        # Convert the sparse matrix to COO format
        data_coo = data_mat.tocoo()
        
        # Transfer data to the R environment
        ro.globalenv["data_row"] = ro.IntVector(data_coo.row + 1)  # R is 1-indexed
        ro.globalenv["data_col"] = ro.IntVector(data_coo.col + 1)
        ro.globalenv["data_values"] = ro.FloatVector(data_coo.data)
        ro.globalenv["nrow"] = data_mat.shape[0]
        ro.globalenv["ncol"] = data_mat.shape[1]
        ro.globalenv["seed"] = seed

    else:
        raise ValueError("The input matrix is not sparse")

    # Run scDblFinder in R
    ro.r('''
        library(scDblFinder)
        library(SingleCellExperiment)
        library(Matrix)

        set.seed(seed)
        
        # Reconstruct the sparse matrix in R
        sparse_mat <- sparseMatrix(
            i = data_row,
            j = data_col,
            x = data_values,
            dims = c(nrow, ncol)
        )
        
        sce <- scDblFinder(
            SingleCellExperiment(
                list(counts = sparse_mat)
            )
        )
        doublet_score <- sce$scDblFinder.score
        doublet_class <- sce$scDblFinder.class
    ''')

    # Retrieve the results
    doublet_score = np.array(ro.globalenv["doublet_score"])
    doublet_class = np.array(ro.globalenv["doublet_class"])

    # Store the results in the AnnData object
    adata.obs["scDblFinder_score"] = doublet_score
    adata.obs["scDblFinder_class"] = doublet_class
    print(adata.obs["scDblFinder_class"].value_counts())
    
    adata_singlet = adata[adata.obs['scDblFinder_class'] == 'singlet'].copy()
    
    return adata_singlet
