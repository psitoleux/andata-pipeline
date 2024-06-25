import numpy as np 
import scanpy as sc
from abc import ABC
from modules.abmodule import Module
    
class Pipeline(Module):
    
    def __init__(self, adata : sc.AnnData, config : dict):
        
        self.raw = adata
        self.config = config
        
        self.outlier_keys = config["outlier_keys"]
        
        
        self.ambient = config["ambient"]
        self.doublets = config["doublets"]
        
        self.normalization = config["normalization"]
        self.feature_selection = config["feature_selection"]
        self.dim_reduction = config["dim_reduction"]
        
        self.batch_corr = config["batch_corr"]
        
        self.adata = adata
        
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
        
        pass
    
    
    