import numpy as np 
from scipy.stats import median_abs_deviation as mad
import scanpy as sc
import rpy2

from numba import njit
from tqdm.auto import trange
from joblib import Parallel, delayed
import scipy.sparse as sp


def is_outlier(adata: sc.AnnData, metric: str, nmads: int) -> np.ndarray:
    """
    Checks if data points are outliers based on a given metric and number of
    median absolute deviations to consider as outliers.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.
    metric : str
        Name of the metric used to check for outliers.
    nmads : int
        Number of median absolute deviations to consider as outliers.

    Returns
    -------
    np.ndarray
        A boolean array indicating whether each data point is an outlier or not.
    """
    # Get the metric values
    M = adata.obs[metric]
    
    # Calculate the outlier cutoff based on the median and number of MADs
    outlier_cutoff = np.median(M) + nmads * mad(M)
    
    # Check if each data point is an outlier
    outlier = (M < np.median(M) - nmads * mad(M)) | (outlier_cutoff < M)
    
    return outlier

def reorder_labels_and_matrix(matrix, labels):
    
    import rpy2.robjects.numpy2ri
    from rpy2.robjects.packages import importr
        
    rpy2.robjects.numpy2ri.activate()

    # Import R package seriation
    seriation = importr('seriation')


    # Convert Python numpy matrix to R matrix
    r_matrix = rpy2.robjects.r.matrix(matrix, nrow=matrix.shape[0], ncol=matrix.shape[1])
    
    # Call seriate function from seriation package
    ordered_indices = seriation.seriate(r_matrix, method = 'PCA_angle')
    print(ordered_indices[1])
    # Convert R indices back to Python + making sure they have the right shape
    ordered_indices = (np.array(ordered_indices[1])-1).flatten() 

    # Convert ordered indices back to Python
    ordered_matrix = matrix[np.ix_(ordered_indices, ordered_indices)]
    ordered_labels = [labels[i] for i in ordered_indices]

    return ordered_matrix, ordered_labels


def shuffle_row(row):
    if sp.issparse(row):
        # Convert the sparse row to COO format for easy manipulation
        row_coo = row.tocoo()
        
        # Get non-zero indices and values
        idx, values = row_coo.col, row_coo.data
        
        # Shuffle the values
        np.random.shuffle(values)
        
        # Reconstruct the sparse row
        shuffled_row = sp.csr_matrix((values, (np.zeros_like(idx), idx)), shape=(1, row.shape[1]))
        
        return shuffled_row
    else:
        # Handle dense arrays, though not expected in a sparse matrix
        np.random.shuffle(row)
        return row

def shuffle_adata(adata: sc.AnnData) -> sc.AnnData:
    # Parallel processing of rows
    shuffled_rows = Parallel(n_jobs=-1)(
        delayed(shuffle_row)(adata.X[i]) for i in trange(adata.X.shape[0])
    )
    
    # Reconstruct the sparse matrix from the shuffled rows
    adata.X = sp.vstack(shuffled_rows)
    
    return adata