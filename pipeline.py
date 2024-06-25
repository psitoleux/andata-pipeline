import numpy as np 
import scanpy as sc
from abc import ABC



class DoubletFinder(Module):
    def __init__(self, method : str, **kwargs):
        
        pass
    
    def call(self, adata : sc.AnnData) -> sc.AnnData:
        
        pass
    
class AmbientRemover(Module):
    def __init__(self, method : str, **kwargs):

        pass
    
    def call(self, adata : sc.AnnData) -> sc.AnnData:
        
        pass
    
    
class FeatureSelection(Module):
    
    def __init__(self, method : str, **kwargs):
        
        pass
    
    def call(self, adata : sc.AnnData) -> sc.AnnData:
        
        pass
    

class Normalization(Module):
    
    def __init__(self, method : str, **kwargs):
        
        pass
    
    def call(self, adata : sc.AnnData) -> sc.AnnData:
        
        pass
    
class DimReduction(Module):
    
    def __init__(self, method : str, **kwargs):
        
        pass
    
    def call(self, adata : sc.AnnData) -> sc.AnnData:
        
        pass
    
class BatchCorrection(Module):
    
    def __init__(self, method : str, **kwargs):
        
        pass
    
    def call(self, adata : sc.AnnData) -> sc.AnnData:
        
        pass

class Pipeline:
    
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
            
    def delraw(self) -> None:
        del self.adata.raw
         
    """
            
        # mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
        
    """
    def run(self):
        
        pass
    
    
    