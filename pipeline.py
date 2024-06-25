import numpy as np 
import scanpy as sc
from abc import ABC
from modules.abmodule import Module

config = {
    "outlier_keys" : ["n_counts", "n_genes"],
    
    "ambient" : None,
    "doublets" : None,
    
    "normalization" : True,
    "feature_selection" : True,
    "dim_reduction" : True,
    
    "batch_corr" : True
}
    
class Pipeline(Module):
    
    def __init__(self, adata: sc.AnnData, config: dict):
        self.raw = adata
        self.adata = adata
        self.config = config
        self.outlier_keys = config.get("outlier_keys")
        
        self.ambient = None if config.get("ambient") is None else ""
        
        self.doublets = self._get_doublets_method(config.get("doublets"))
        self.normalization = self._get_normalization_method(config.get("normalization"))
        self.feature_selection = self._get_feature_selection_method(config.get("feature_selection"))
        self.dim_reduction = self._get_dim_reduction_method(config.get("dim_reduction"), config)
        self.batch_corr = self._get_batch_corr_method(config.get("batch_corr"), config)
    
    def _get_doublets_method(self, doublets_config):
        if doublets_config is None:
            return None
        elif doublets_config == "scdbl":
            from modules.doublet import scdbl
            return scdbl
        else:
            return None

    def _get_normalization_method(self, normalization_config):
        if normalization_config is None:
            return None
        elif normalization_config == "log1p":
            from modules.norm import log1p
            return log1p
        elif normalization_config == "pearson_residuals":
            return sc.experimental.pp.normalize_pearson_residuals
        elif normalization_config == "sanity":
            from modules.norm import sanity_normalization
            return sanity_normalization
        else:
            return None
    
    def _get_feature_selection_method(self, feature_selection_config):
        if feature_selection_config is None:
            return None
        elif feature_selection_config == "hvg":
            return sc.pp.highly_variable_genes
        elif feature_selection_config == "deviance":
            from modules.featsel import deviance
            return deviance
        else:
            return None

    def _get_dim_reduction_method(self, dim_reduction_config, config):
        if dim_reduction_config is None:
            return None
        elif dim_reduction_config == "pca":
            return lambda adata: sc.pp.pca(adata, n_comps=config.get("pca_n_comps", 50))
        elif dim_reduction_config == "glmpca":
            from modules.glmpca import glmpca
            return lambda adata: glmpca(adata, n_comps=config.get("pca_n_comps", 50))
        else:
            return None
    
    def _get_batch_corr_method(self, batch_corr_config, config):
        if batch_corr_config is None:
            return None
        elif batch_corr_config == "bbknn":
            from modules.batchcorr import bbknn
            return lambda adata: bbknn(adata, batch_key=config.get("batch_key"))
        else:
            return None
        
        
    def tag(self, tags : dict) -> None:
        """
        A function to tag variables in the AnnData object based on the provided tags dictionary.
        
        Parameters:
            tags (dict): A dictionary where the keys are variables to tag, and the values are tuples 
                         containing the pattern to match and a flag (0 for startswith, 1 for contains).
        
        Returns:
            None
        """
        for k, v in tags.items():
            if v[1] == 0:
                self.adata.var[k] = adata.var_names.str.startswith(v[0])
            else: 
                self.adata.var[k] = adata.var_names.str.contains(v[0])
            
            
        # note for mitochondrial ribo and hb we have tags = {"mt" : ("MT-", 0), "ribo" : ("RPS|RPL", 0), "hb" : ("^HB[^(P)]", 1)}
    def delraw(self) -> None:
        del self.adata.raw
         
    def call(self):
        
        for step in [self.normalization, self.feature_selection, self.dim_reduction, self.batch_corr]:
            if step is not None:
                self.adata = step(self.adata)
                
        return self.adata
    
    
    