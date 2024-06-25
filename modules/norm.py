from abmodule import Module
import scanpy as sc


def log1p(adata : sc.AnnData, L : int = 1e4) -> sc.AnnData:
    sc.pp.normalize_total(adata, target_sum=L)
    sc.pp.log1p(adata)
    
    return adata