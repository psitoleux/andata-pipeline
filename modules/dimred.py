import scanpy as sc
from rpy2.robjects.packages import importr
import rpy2.robjects as ro


def glmpca(adata : sc.AnnData) -> sc.AnnData:
    # TODO apply GLMPCA + write it in the adata 
    
    fastglmpca = importr("fastglmpca")
    
    
    ro.r('''
        library(fastglmpca)
        
        
        
    ''')

 
    pass 

