import scanpy as sc
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import numpy as np
import scipy.sparse as sp

def scdbl(adata: sc.AnnData, seed: int = 123) -> sc.AnnData:
    """
    Run scDblFinder on the input AnnData object.

    Args:
        adata (AnnData): Input AnnData object with normalized counts.
        seed (int): Random seed for reproducibility. Default is 123.

    Returns:
        AnnData: AnnData object with singlet cells stored in the 'adata_singlet' object.
    
    Remark:
        Takes around 15 minutes to run Parks et al. (2020) (zenodo.org/records/5500511)
        thymic development dataset (~250k cells x ~30k genes).
        On an Apple M3 Pro 12 cores 36 GB RAM.
    """
    
    # Import R libraries
    scDblFinder = importr("scDblFinder")
    SingleCellExperiment = importr("SingleCellExperiment")
    Matrix = importr("Matrix")
    
    # Prepare the data matrix
    data_mat = adata.X.T  # Transpose the matrix
    
    # Check if the input matrix is sparse
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
        # Load R libraries
        library(scDblFinder)
        library(SingleCellExperiment)
        library(Matrix)

        # Set random seed for reproducibility
        set.seed(seed)
        
        # Reconstruct the sparse matrix in R
        sparse_mat <- sparseMatrix(
            i = data_row,
            j = data_col,
            x = data_values,
            dims = c(nrow, ncol)
        )
        
        # Run scDblFinder
        sce <- scDblFinder(
            SingleCellExperiment(
                list(counts = sparse_mat)
            )
        )
        
        # Extract the doublet scores and classes
        doublet_score <- sce$scDblFinder.score
        doublet_class <- sce$scDblFinder.class
    ''')

    # Retrieve the results
    doublet_score = np.array(ro.globalenv["doublet_score"])
    doublet_class = np.array(ro.globalenv["doublet_class"])

    # Store the results in the AnnData object
    adata.obs["scDblFinder_score"] = doublet_score
    adata.obs["scDblFinder_class"] = doublet_class
    
    # Print the counts of singlet and doublet cells
    print(adata.obs["scDblFinder_class"].value_counts())
    
    # Return the AnnData object with singlet cells
    adata_singlet = adata[adata.obs['scDblFinder_class'] == 'singlet'].copy()
    
    return adata_singlet
