import numpy as np 
from scipy.stats import median_abs_deviation as mad
import scanpy as sc


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

