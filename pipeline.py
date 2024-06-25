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
    
    def __init__(self, adata : sc.AnnData, config : dict):
        
        self.raw = adata
        self.adata = adata

        self.config = config
        
        self.outlier_keys = config["outlier_keys"]
         
        if config["ambient"] is None:
            self.ambient = None
        elif config["ambient"] == "":
            
        
        if config["doublets"] is None:
            self.doublets = None
        elif config["doublets"] == "scdbl":
            from modules.doublet import scdbl
            self.doublets = scdbl
        
        if config["normalization"] is None:
            self.normalization = None
        elif config['dim_reduction'] == "glmpca":
            self.normalization = None
        elif config["normalization"] == "log1p": 
            from modules.norm import log1p
            self.normalization = log1p
        elif config["normalization"] == "pearson_residuals":
            self.normalization = sc.experimental.pp.normalize_pearson_residuals
        else:
            self.normalization = None
        
        if config["feature_selection"] is None:
            self.feature_selection = None
        elif config["feature_selection"] == "hvg":
            self.feature_selection = sc.pp.highly_variable_genes
        elif config["feature_selection"] == "deviance":
            from modules.featsel import deviance
            self.feature_selection = deviance
        
        if config["dim_reduction"] is None:
            self.dim_reduction = None
        elif config["dim_reduction"] == "pca":
            self.dim_reduction = lambda adata : sc.pp.pca(adata, n_comps = config["pca_n_comps"])
        elif config["dim_reduction"] == "glmpca":
            from modules.glmpca import glmpca
            self.dim_reduction = lambda adata : glmpca(adata, n_comps = config["pca_n_comps"])
        
        
        if config["batch_corr"] is None:
            self.batch_corr = None
        elif config["batch_corr"] == "bbknn":
            from modules.batchcorr import bbknn
            self.batch_corr = lambda adata : bbknn(adata, batch_key = config["batch_key"])
            
        
        
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
    
    
    