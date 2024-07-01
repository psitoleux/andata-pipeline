import numpy as np 
import scanpy as sc
from parser import get_config

from utils import is_outlier

import os
from urllib.request import urlretrieve


config = {
    "outlier_keys" : ["n_counts", "n_genes"],
    
    "ambient" : None,
    "doublets" : None,
    
    "normalization" : True,
    "feature_selection" : True,
    "dim_reduction" : True,
    
    "batch_corr" : True
}
    
class Pipeline():
    
    def __init__(self, config: dict = None) -> None:
        
        self.config = get_config() if config is None else config
 
        self.input = config.get("input")
        self.output = config.get("output")
        self.online_link = config.get("online_link")
        self.figures_dir = './figures/'
        
        if self.input is None and self.online_link is not None:
            print("Downloading data...")
            self.input = self._download_file(self.online_link)
        
        if not os.path.isfile(self.input):
            if self.online_link is not None:
                print(f"File {self.input} not found locally. Downloading data...")
                self.input = self._download_file(self.online_link)
            else:
                raise FileNotFoundError(f"Input file {self.input} not found and no online link provided.")        
        
        print("Loading data...") 
        adata = sc.read(self.input, cache = True)
    
        print('Making gene names unique...')
        adata.var_names_make_unique()
        
        self.raw = adata.copy()
        self.adata = adata
        self._update_config()        
        
        del adata
    
    def get_new_config(self) -> None:
        
        self.config = get_config()
        self._update_config()
        self._reset_adata()
        
        return None
    
    def _update_config(self) -> None:
        """
        Update the configuration parameters of the pipeline.

        This function updates the configuration parameters of the pipeline
        based on the values specified in the self.config dictionary.

        Returns:
            None
        """

        # Update outlier keys
        self.outlier_keys = self.config.get("outlier_keys")

        # Update ambient method
        self.ambient = self._get_ambient_method(self.config.get("ambient"))

        # Update doublets method
        self.doublets = self._get_doublets_method(self.config.get("doublets"))

        # Update normalization method
        self.normalization = self._get_normalization_method(self.config.get("normalization"))

        # Update feature selection method
        self.feature_selection = self._get_feature_selection_method(self.config.get("feature_selection"))

        # Update dimensionality reduction method
        self.dim_reduction = self._get_dim_reduction_method(self.config.get("dim_reduction"))

        # Update batch correction method
        self.batch_corr = self._get_batch_corr_method(self.config.get("batch_corr"))

        # Update visualization method
        self.visualization = self._get_visualization_method(self.config.get("visualization"))

        # TODO correct ordering of steps ?
        self.steps = [self.ambient, self.doublets, self.normalization, self.feature_selection,
                    self.dim_reduction, self.batch_corr]

        return None
    
    

    def _reset_adata(self) -> None:
        self.adata = self.raw.copy()
        
        return None
         
         
    def _download_file(self, url: str) -> str:
        file_name = url.split('/')[-1]
        urllib.request.urlretrieve(url, file_name)
        return file_name
    
    def _get_ambient_method(self, ambient_config):
        if ambient_config is None:
            return None
        elif ambient_config == "soupx":
            from modules.ambient import soupx
            return soupx
        else:
            return None
    
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
        elif feature_selection_config == "pearson":
            return sc.experimental.pp.highly_variable_genes
        elif feature_selection_config == "deviance":
            from modules.featsel import deviance
            return deviance
        else:
            return None

    def _get_dim_reduction_method(self, dim_reduction_config):
        if dim_reduction_config is None:
            return None
        elif dim_reduction_config == "pca":
            return lambda adata: sc.pp.pca(adata, n_comps=self.config.get("pca_n_comps", 50))
        elif dim_reduction_config == "glmpca":
            from modules.glmpca import glmpca
            return lambda adata: glmpca(adata, n_comps=self.config.get("pca_n_comps", 50))
        else:
            return None
    
    def _get_batch_corr_method(self, batch_corr_config):
        if batch_corr_config is None:
            return None
        elif batch_corr_config == "bbknn":
            from scanpy.external.pp import bbknn 
            return lambda adata: bbknn(adata, batch_key=self.config.get("batch_key"))
        elif batch_corr_config == "bbknn_reg":
            from modules.batch_corr import bbknn_reg
            return lambda adata: bbknn_reg(adata, batch_key=self.config.get("batch_key"))
        else:
            return None
        
    def _get_visualization_method(self, visualization_config):
        if visualization_config is None:
            return None
        elif visualization_config == "umap":
            return lambda adata, n_comp=2: sc.tl.umap(adata, n_components=n_comp)
        # TODO PacMap 
        else:
            return None
        
    
    def outliers(self) -> None:
        
        print("Removing outliers...") 
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=self.outlier_keys, inplace=True, percent_top=[20], log1p=True)
        
        self.adata.obs['outlier'] =  (is_outlier(self.adata, "log1p_total_counts", 5)
    | is_outlier(self.adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(self.adata, "pct_counts_in_top_20_genes", 5))
        
        print("Removing MT outliers...") 
        if 'mt' in self.outlier_keys:
            self.adata.obs['mt_outlier'] =  is_outlier(self.adata, "pct_counts_mt", 3) | (self.adata.obs["pct_counts_mt"] > 8)
            
        self.adata = self.adata[(~self.adata.obs.outlier) & (~self.adata.obs.mt_outlier)].copy()
        
        print('Outliers removed') 
        
    def tag(self, tags : dict) -> None:
        """
        A function to tag variables in the AnnData object based on the provided tags dictionary.
        
        Parameters:
            tags (dict): A dictionary where the keys are variables to tag, and the values are tuples 
                         containing the pattern to match and a flag (0 for startswith, 1 for contains).
        
        Returns:
            None
            
        Example:
            tags = {"mt" : ("MT-", 0), "ribo" : (("RPS","RPL"), 0), "hb" : ("^HB[^(P)]", 1)} 

        """
        for k, v in tags.items():
            if v[1] == 0:
                self.adata.var[k] = self.adata.var_names.str.startswith(v[0])
            else: 
                self.adata.var[k] = self.adata.var_names.str.contains(v[0])
            
            
    def delraw(self) -> None:
        del self.adata.raw
         
    def preprocess(self) -> sc.AnnData:
        
        self.tag({"mt" : ("MT-", 0), "ribo" : (("RPS","RPL"), 0), "hb" : ("^HB[^(P)]", 1)})
        self.outliers()
        
        print('running normalization, feature selection, dim reduction, and batch correction...')

        if self.normalization is not None:
            self.normalization(self.adata)
        
        if self.ambient is not None:
            self.adata = self.ambient(self.adata)
            
        if self.doublets is not None:
            self.adata = self.doublets(self.adata)
            
        
        
        for step in [ self.feature_selection, self.dim_reduction, self.batch_corr]:
            if step is not None:
                step(self.adata)
                
        return self.adata
    
    def visualize(self, color_key = 'method') -> None:
       
        print('outputting pca plot...') 
        fig_pca = sc.pl.pca(self.adata, color = color_key, annotate_var_explained = True, return_fig = True) # plot 2D pca 
        
        
        print('creating 2d embedding plot...')
        self.visualization(adata=self.adata) # create umap-like plots
        sc.pl.umap(self.adata, color = color_key, )
        
        # TODO allow 3D plots with plotly
        from plots import pca_3d
        print('creating 3d interactive pca plot...')
        pca_3d(self.adata, color_key = color_key) 
        
        return None
    
    def analysis(self) -> None:
        
        return None
