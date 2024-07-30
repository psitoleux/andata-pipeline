import numpy as np 
import scanpy as sc

from typing import Tuple
from scipy.stats import pearsonr, spearmanr



def correlations_along_vector(
    adata: sc.AnnData,
    v : np.ndarray
        ) -> Tuple[np.ndarray, tuple, tuple]:

    """
    
    """
    
    v = v / np.sqrt(np.sum(v**2))
    pearson_correlations = np.zeros(adata.obsm['X_pca'].shape[1])
    pearson_p_values = np.zeros(adata.obsm['X_pca'].shape[1])
    
    spearman_correlations = np.zeros(adata.obsm['X_pca'].shape[1])
    spearman_p_values = np.zeros(adata.obsm['X_pca'].shape[1])

    
    for i in range(adata.obsm['X_pca'].shape[1]):
        w = adata.obsm['X_pca'][:,i] / np.linalg.norm(adata.obsm['X_pca'][:,i])
        
        pearson_correlations[i], pearson_p_values[i] = pearsonr(v, w)
        spearman_correlations[i], spearman_p_values[i] = spearmanr(v, w)
        
    
    
    # Compute the most correlated axis to vector v 
    
    variation_axis = np.zeros(adata.obsm['X_pca'].shape[1])
    
    for i in range(adata.obsm['X_pca'].shape[1]):
        r = pearson_correlations[i]
        variation_axis[i] = np.dot(v, adata.obsm['X_pca'][:,i] / np.linalg.norm(adata.obsm['X_pca'][:,i])  )
        
    variation_axis /= np.sqrt(np.sum(variation_axis**2))
    
    return variation_axis, (pearson_correlations, pearson_p_values), (spearman_correlations, spearman_p_values)
