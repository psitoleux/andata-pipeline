import numpy as np
import pandas as pd
import scanpy as sc
from parser import get_config

from utils import is_outlier

import os, json
from datetime import datetime
from urllib.request import urlretrieve
import pickle

import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go


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
    
    def __init__(self, config : dict = None) -> None:
        
        
        self.config = get_config() if config is None else config
 
        self.input_file = config.get("input")
        
        if self.input_file == 'thymus_raw':
            self.input_file = 'HTA07.A01.v02.entire_data_raw_count.h5ad'
        elif self.input_file == 'thymus_fig1': 
            self.input_file = 'HTA08.v01.A05.Science_human_fig1.h5ad'
            
        
            
        
        
        self.online_link = config.get("online_link")
        
        self._save = config.get("save")
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_dir = config.get("output") + '/' + self.timestamp + '/'

        self.figures_dir = self.output_dir + 'figures/'
     
        # Check if figures directory exists, if not, create it
        
        if self._save:
            os.makedirs(self.output_dir, exist_ok=True)
            os.makedirs(self.figures_dir, exist_ok=True)
    
        if self.input_file is None and self.online_link is not None:
            print("Downloading data...")
            self.input_file = self._download_file(self.online_link)
        
        if not os.path.isfile(self.input_file):
            if self.online_link is not None:
                print(f"File {self.input} not found locally. Downloading data...")
                self.input_file = self._download_file(self.online_link)
            else:
                raise FileNotFoundError(f"Input file {self.input} not found and no online link provided.")        
        

        
        print("Loading data...") 
        adata = sc.read(self.input_file, cache = True)
    
        print('Making gene names unique...')
        adata.var_names_make_unique()
        
        self.raw = adata.copy()
        self.adata = adata
        self._update_config()
        
        
        del adata
        
        if self._save:
            self.save_config()
    
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
        
        # Update timestamp (to have distsinct filenames)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S") 
        
        # Update outlier keys
        self.outlier_keys = self.config.get("outlier_keys")
        self.do_qc = self.config.get("qc")
        
        if not self.do_qc:
            print('! Warning: Outliers will not be excluded.')

        # Update ambient method
        self.ambient = self._get_ambient_method(self.config.get("ambient"))

        # Update doublets method
        self.doublets = self._get_doublets_method(self.config.get("doublets"))

        # Update normalization method
        self.normalization = self._get_normalization_method(self.config.get("normalization"))

        # Update feature selection method
        self.n_top_genes = 1_000 #self.config.get("n_top_genes")
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
        elif normalization_config == "plog1p":
            from modules.norm import pure_log1p
            return pure_log1p
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
            from modules.featsel import highly_variable_genes
            return lambda adata, n_top_genes=self.n_top_genes: highly_variable_genes(adata, n_top_genes=n_top_genes)
        elif feature_selection_config == "pearson":
            return sc.experimental.pp.highly_variable_genes
        elif feature_selection_config == "deviance":
            from modules.featsel import deviance
            return deviance
        elif feature_selection_config == "heg":
            from modules.featsel import highly_expressed_genes
            return lambda adata, n_top_genes=self.n_top_genes: highly_expressed_genes(adata, n_top_genes=n_top_genes)
        elif feature_selection_config == "moeg":
            from modules.featsel import most_often_expressed_genes
            return lambda adata, n_top_genes=self.n_top_genes: most_often_expressed_genes(adata, n_top_genes=n_top_genes)
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
        
    def _get_cycle_genes(self) -> list:
         
        with open('gene_lists/cycle_genes.pkl', 'rb') as f:
            self.cycle_genes = pickle.load(f)
        
        return self.cycle_genes
    
    def _get_yamanaka_factors(self) -> list:
        
        yamanaka_factors = [
            "POU5F1",  # Also known as OCT4
            "SOX2",
            "KLF4",
            "MYC"      # Also known as c-MYC
            ]
        
        return yamanaka_factors
    
    
    def _save_gene_lists(self, tosave : list, filename : str) -> None:
        
        with open("gene_lists/" + filename + ".pkl", 'wb') as f:
            pickle.dump(tosave, f)
            
        return None
    
    
            
    
    def outliers(self) -> None:
        
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=self.outlier_keys, inplace=True, percent_top=[20], log1p=True)
        
        self.adata.obs['outlier'] =  (is_outlier(self.adata, "log1p_total_counts", 5)
    | is_outlier(self.adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(self.adata, "pct_counts_in_top_20_genes", 5))
       
        if 'mt' in self.outlier_keys:
            self.adata.obs['mt_outlier'] =  is_outlier(self.adata, "pct_counts_mt", 3) | (self.adata.obs["pct_counts_mt"] > 8)
          
          
        
        if self.do_qc:
            print("Removing outliers...") 
  
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
            self.adata = self.normalization(self.adata)
        
        if self.ambient is not None and self.do_qc:
            self.adata = self.ambient(self.adata)
            
        if self.doublets is not None and self.do_qc:
            self.adata = self.doublets(self.adata)
            
        if self.feature_selection is not None:
            self.adata = self.feature_selection(self.adata) 
        
        for step in [self.dim_reduction, self.batch_corr]:
            if step is not None:
                step(self.adata)
                
                
        
        nunique = self.adata.obs.nunique()

        self.adata.obs = self.adata.obs.drop(nunique[nunique == 1].index, axis=1)
        
        print('preprocessing complete')
         
        return self.adata
    
    def test(self) -> None:
        
        print('testing...')
        print(self._get_cycle_genes())
        print('passed')
    
    def visualize(self, n_dim = 2) -> None:
        from plots import scatter3D, pca3D
        
        if not ('neighbors' in self.adata.uns):
            sc.pp.neighbors(self.adata)
         
        self.visualization(self.adata, n_dim)            
            
        print('outputting pca plot...') 
        fig_pca = sc.pl.pca(self.adata, color = 'method', annotate_var_explained = True, return_fig = True) # plot 2D pca 
        self.save_mpl(fig_pca)
         
        X_ = np.array([self.adata.obsm['X_umap'][:, 0], self.adata.obsm['X_umap'][:, 1], self.adata.obsm['X_pca'][:, 0]]).T
        fig_umap_pca = scatter3D(X_, colors = self.adata.obs, title='UMAP & PC1', labels = ['UMAP 1', 'UMAP 2', 'PC1'])
        
        self.save_plotly(fig_umap_pca, 'umap_pca')
        
        fig_pca_3D = pca3D(self.adata, idx = [0, 1, 2])
        
        self.save_plotly(fig_pca_3D, 'pca3d')    
        
        #if n_dim == 3:
            
        #    pca3D(self.adata)
        #    scatter3D(self.adata.obsm['X_umap'], colors = self.adata.obs, title='UMAP')

        return None
    
    def analysis(self) -> None:
        
        from plots import set_matplotlib_style
        
        print('setting matplotlib style...')
        set_matplotlib_style()
        
        print('excluding cycle genes...')
        cycle_genes = self._get_cycle_genes()
        high_var_no_cycle_idx = (self.adata.var['highly_variable'] != (self.adata.var['highly_variable'] 
                                                        & self.adata.var['GeneName'].isin(cycle_genes))) 
        
        
        relevant_genes = self.adata.var['GeneName'][high_var_no_cycle_idx]
        
        print('getting dense array...')
        A = self.adata[:, high_var_no_cycle_idx].X.toarray()
        
        print('computing correlations, covariances...')
        correlations = np.corrcoef(A, rowvar=False)
        covariances = np.cov(A, rowvar=False)
        
        print('inverting covariance to get coupling...')
        coupling_matrix = np.linalg.pinv(covariances)
        
        import rpy2.robjects.numpy2ri
        from rpy2.robjects.packages import importr
        
        rpy2.robjects.numpy2ri.activate()

        # Import R package seriation
        seriation = importr('seriation')
        
        def reorder_labels_and_matrix(matrix, labels):
            # Convert Python numpy matrix to R matrix
            r_matrix = rpy2.robjects.r.matrix(matrix, nrow=matrix.shape[0], ncol=matrix.shape[1])
            
            # Call seriate function from seriation package
            ordered_indices = seriation.seriate(r_matrix, method = 'PCA')
            
            # Convert R indices back to Python + making sure they have the right shape
            ordered_indices = (np.array(ordered_indices)-1).flatten() 

            # Convert ordered indices back to Python
            ordered_matrix = matrix[np.ix_(ordered_indices, ordered_indices)]
            ordered_labels = [labels[i] for i in ordered_indices]

            return ordered_matrix, ordered_labels
        
        if False: 
            from plots import heatmap_with_annotations
            
            corr_ordered, labels_corr  = reorder_labels_and_matrix(correlations, relevant_genes)
            
            fig = heatmap_with_annotations(corr_ordered, labels_corr)
            
            coupling_ordered, labels_coupling = reorder_labels_and_matrix(log1p(coupling_matrix), relevant_genes )
            
            
            log1p = lambda x: np.sign(x) * np.log(1+np.abs(x))
            
            fig = heatmap_with_annotations((coupling_ordered), labels_coupling)
            
            covariances_ordered, labels_cov = reorder_labels_and_matrix(covariances, relevant_genes)
            
            fig = heatmap_with_annotations(log1p(covariances_ordered), labels_cov)
        
        print('plotting variance ratio...')
        fig, ax = plt.subplots(figsize=(16, 9))
        
        var_ratio = self.adata.uns['pca']['variance_ratio']
        n_coms = len(var_ratio)
        plt.plot(range(1, n_coms + 1), var_ratio, 'o-')
        
        def top_k_off_diagonal_indices_symmetric(matrix, k):
            # Ensure the matrix is square
            assert matrix.shape[0] == matrix.shape[1], "The input matrix must be square."
            
            n = matrix.shape[0]
            
            # Create a mask for the upper triangular off-diagonal elements
            mask = np.triu(np.ones((n, n), dtype=bool), k=1)
            
            # Extract the upper triangular off-diagonal elements
            off_diagonal_elements = matrix[mask]
            
            # Get the indices of the sorted elements in descending order
            sorted_indices = np.argsort(off_diagonal_elements)[::-1]
            
            # Select the top k indices
            top_k_indices_flat = sorted_indices[:k]
            
            # Convert flat indices back to 2D indices
            top_k_indices = np.vstack(np.unravel_index(np.where(mask.flatten())[0][top_k_indices_flat], (n, n))).T
            
            return top_k_indices
        
        def plot_top_k_joints(data : np.ndarray, matrix : np.ndarray, genes
                              , title : str = ''
                              , k :int = 9):
            
            top_k_indices = top_k_off_diagonal_indices_symmetric(matrix, k)
            
            return plot_k_joints(data, matrix, genes, top_k_indices, title, k)
            
         
        def plot_k_joints(data : np.ndarray
                          , matrix : np.ndarray
                          , genes : pd.Series
                          , indices : np.ndarray
                          , title : str = ''
                          , k : int = 9):
                
            
            l = int(np.sqrt(k))
            

            fig, axs = plt.subplots(l, l, figsize=(l*4, l*4))

            for i,ax in zip(range(k), axs.flat):
                
                genex, geney = genes.iloc[indices[i, 0]], genes.iloc[indices[i, 1]]


                share_both_zeros = np.sum((data[:, indices[i, 0]] == 0) & (data[:, indices[i, 1]] == 0)) / data.shape[0]
                nb_non_zeros_both = np.sum((data[:, indices[i, 0]] != 0) & (data[:, indices[i, 1]] != 0)) 
                nb_non_zeros_x = np.sum(data[:, indices[i, 0]] != 0) 
                nb_non_zeros_y = np.sum(data[:, indices[i, 1]] != 0)
                
                alpha = 0.5
                if nb_non_zeros_both > 1000:
                    alpha = 0.1
                if nb_non_zeros_both > 10000:
                    alpha = 0.01
                
                
                ax.plot(data[:, indices[i, 0]]
                        ,data[:, indices[i, 1]], 'o', markersize=2, alpha = alpha,
                        label = genex + ' - ' + geney,
                        color = 'C' + str(i % 10)) 
 
                ax.set_title('share both zeros ' + '{0:.4f}'.format(share_both_zeros) 
                             + ' \n  nb both non-zeros ' + str(nb_non_zeros_both) 
                             + ' \n matrix value ' + '{0:.2f}'.format(matrix[indices[i, 0], indices[i, 1]]))
                ax.set_xlabel(genex + ' - ' + str(nb_non_zeros_x) + ' non-zeros' )
                ax.set_ylabel(geney + ' - ' + str(nb_non_zeros_y) + ' non zeros' )
                
                ax.set_aspect('equal')
                
            fig.suptitle(title) 
            plt.tight_layout()

            return fig
        
        print('plotting joints...')

        self.save_mpl(
            plot_top_k_joints(A, coupling_matrix, relevant_genes, title='coupling-matrix')
        )
        
        self.save_mpl(
            plot_top_k_joints(A, -coupling_matrix, relevant_genes, title = 'negative-coupling-matrix')
        )
        
        self.save_mpl(
            plot_top_k_joints(A, correlations, relevant_genes, title = 'correlations')
        )
        
        self.save_mpl(
            plot_top_k_joints(A, -correlations, relevant_genes, title = 'negative-correlations')
        )
        
        self.save_mpl(
            plot_top_k_joints(A, covariances, relevant_genes, title = 'covariances')
        )
        
        self.save_mpl(
            plot_top_k_joints(A, -covariances, relevant_genes, title = 'negative-covariances')
        )
        
        def random_offdiag_idx(N, k):
            
            indices = [(i,j) for i in range(N) for j in range(i+1, N)]
            iidx = np.random.choice(len(indices), size=k, replace=False).astype(int)
            
            return indices[iidx]
            
            
        def plot_k_random_joints(data : np.ndarray, matrix : np.ndarray, genes, k = 9):
                
            indices = random_offdiag_idx(data.shape[1], k)
                
            return plot_k_joints(data, matrix, genes, indices, title = 'random-pairs', k = k)

        
        self.save_mpl(
            plot_k_random_joints(A, coupling_matrix, relevant_genes)
        )

        
        return None
    
    def save_mpl(self, fig : plt.Figure, format : str = 'png', title : str = None) -> None:
        
        if title is None:
            title = fig.get_suptitle()
        
        fig.savefig(self.figures_dir + '/' + title + '.' + format
                    , bbox_inches='tight'
                    ,format = format
                    , dpi = 600
                    )
        del fig
        
        
    def save_plotly(self, fig : go.Figure, name : str) -> None:
        
        fig.write_html(self.figures_dir + '/' + name + '.html')
        
        del fig
    
    def save_config(self) -> None:
        # Save configuration
        config_filename = os.path.join(self.output_dir, f"config_{self.timestamp}.json")
        print(f"Saving configuration to {config_filename}...")
        with open(config_filename, 'w') as config_file:
            # Write the configuration to the file
            json.dump(self.config, config_file, indent=4)

    def save_adata(self) -> None:
        """
        Saves the AnnData object and the configuration to separate files in a subdirectory.
        The subdirectory is named "output_<timestamp>" and is created in the output directory specified in the configuration.
        """
        print(f"Saving AnnData to {adata_filename}...")
        self.adata.write(adata_filename)
        
        print("Save complete.")