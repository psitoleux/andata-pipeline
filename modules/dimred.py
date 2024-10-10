import scanpy as sc
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

from scipy import sparse as sp
from sklearn.utils.sparsefuncs import mean_variance_axis 

def glmpca_(adata : sc.AnnData, n_comps : int = 50) -> sc.AnnData:
    # TODO apply GLMPCA + write it in the adata )
    
    
    data_mat = adata.X.T
    
    ro.globalenv['n_comps'] = n_comps
    
    if sp.issparse(data_mat):
        # Convert the sparse matrix to COO format
        data_coo = data_mat.tocoo()
        
        # Transfer data to the R environment
        ro.globalenv["data_row"] = ro.IntVector(data_coo.row + 1)  # R is 1-indexed
        ro.globalenv["data_col"] = ro.IntVector(data_coo.col + 1)
        ro.globalenv["data_values"] = ro.FloatVector(data_coo.data)
        ro.globalenv["nrow"] = data_mat.shape[0]
        ro.globalenv["ncol"] = data_mat.shape[1]
 

 
    
    ro.r('''
        library(fastglmpca)
        library(Matrix)
        
        
        sparse_mat <- sparseMatrix(
            i = data_row,
            j = data_col,
            x = data_values,
            dims = c(nrow, ncol)
        )

        
        fit0 <- init_glmpca_pois(sparse_mat, K = n_comps)
        
        print('fit0 created')
        
        fit <- fit_glmpca_pois(sparse_mat, fit0, verbose = TRUE)
        
        
    ''')
    
    
    U, V, d = ro.r('fit$U'), ro.r('fit$V'), ro.r('fit$d')
    
    andata.obsm['X_pca'] = np.array(V)
    andata.varm['PCs'] = np.array(U)
    andata.uns['pca']['variance'] = np.array(d)
    
    _, total_variance = mean_variance_axis(np.log(adata.X+1), axis = 0 )
    
    adata.uns['pca']['variance_ratio'] = np.array(d) / np.log( np.exp(total_variance) + 1 )
    
        
    
 
    return adata

def hpca(adata : sc.AnnData, n_comps : int = 50) -> sc.AnnData:
    
    X = adata[:, adata.var.highly_variable].X
    
    X_centered = X - X.mean(axis = 0)
    
    cov = X_centered.T @ X_centered
    
    w, 
    
    
    
    


def glmpca(adata : sc.AnnData, n_comps : int = 50) -> sc.AnnData:
    
    from glmpca.glmpca import glmpca as GLMPCA
    
    
    out = GLMPCA(adata.X.T, L = n_comps)
    
    U = out['factors']
    V = out['loadings']
    
    print(U.shape, V.shape)
    
    print(np.linalg.norm(U[0,:]))
    print(np.linalg.norm(U[:,0]))
    
    print('\n')
    
    
    print(np.linalg.norm(V[0,:]))
    print(np.linalg.norm(V[:,0]))
    
    return None