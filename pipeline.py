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
import seaborn as sns
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
        self.output_dir = config.get("output")  + self.timestamp + '/'

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
        self.seed = self.config.get("seed")
        
        # Update outlier keys
        self.outlier_keys = self.config.get("outlier_keys")
        self.do_qc = self.config.get("qc")
        
        if not self.do_qc:
            print('! Warning: Outliers will not be excluded.')
            
        self.shuffle = self.config.get("shuffle")
        
        self.obs_mask = self.config.get("obs_mask")
        
        if bool(self.obs_mask):
            for key in self.obs_mask:
                if key not in self.adata.obs.columns:
                    raise KeyError(f"Key {key} not found in adata.obs.")
                
                self.adata = self.adata[self.adata.obs[key].isin(self.obs_mask[key])]
                
                
        if self.shuffle:
            print("Shuffling data...")
            np.shuffle(self.adata.X)
                
        
        # Update ambient method
        self.ambient = self._get_ambient_method(self.config.get("ambient"))

        # Update doublets method
        self.doublets = self._get_doublets_method(self.config.get("doublets"))
        self.filter_doublets = self.config.get("filter_dbl")

        # Update normalization method
        self.normalization = self._get_normalization_method(self.config.get("normalization"))

        # Update feature selection method
        self.n_top_genes = self.config.get("n_top_genes")
        self.feature_selection = self._get_feature_selection_method(self.config.get("feature_selection"))

        # Update dimensionality reduction method
        self.dim_reduction = self._get_dim_reduction_method(self.config.get("dim_red"))

        # Update batch correction method
        self.batch_corr = self._get_batch_corr_method(self.config.get("batch_corr"))

        # Update visualization method
        self.visualization = self._get_visualization_method(self.config.get("viz"))

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
            return lambda adata: scdbl(adata=adata, filter = self.filter_doublets, seed=self.seed)
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
        elif normalization_config == "pearson":
            def pearson_residuals(adata):
                print('applying pearson residuals recipe...')
                sc.experimental.pp.recipe_pearson_residuals(adata
                                                            , random_state=self.seed
                                                            , inplace=True
                                                            , n_top_genes=self.n_top_genes
                                                            , theta=self.config.get("pearson_overdis", 100))
                print('var columns: ', adata.var.columns)
                #print('doing pearson residual normalization...')
                from modules.norm import log1p
                
                #sc.experimental.pp.normalize_pearson_residuals(adata, inplace=True)
                

                return log1p(adata) 
            return pearson_residuals
        elif normalization_config == "sanity":
            from modules.norm import sanity_normalization
            return sanity_normalization
        else:
            return None
    
    def _get_feature_selection_method(self, feature_selection_config):
        if feature_selection_config is None:
            return None
        elif self.config.get("dim_red") == "glmpca":
            from modules.featsel import highly_variable_genes
            def glmpca_feature_selection(adata, n_top_genes=self.n_top_genes):
                
                adata.X = np.log1p(adata.X)
                adata = highly_variable_genes(adata, n_top_genes=n_top_genes)
                adata.X = np.expm1(adata.X)
                
                return adata 
            return glmpca_feature_selection
        elif feature_selection_config == "hvg" or feature_selection_config == "seurat":
            from modules.featsel import highly_variable_genes
            return lambda adata, n_top_genes=self.n_top_genes: highly_variable_genes(adata, n_top_genes=n_top_genes)
        elif feature_selection_config == "deviance":
            from modules.featsel import deviance
            return deviance
        elif feature_selection_config == "heg":
            from modules.featsel import highly_expressed_genes
            return lambda adata, n_top_genes=self.n_top_genes: highly_expressed_genes(adata, n_top_genes=n_top_genes)
        elif feature_selection_config == "moeg":
            from modules.featsel import most_often_expressed_genes
            return lambda adata, n_top_genes=self.n_top_genes: most_often_expressed_genes(adata, n_top_genes=n_top_genes)
        elif feature_selection_config == "seurat_v3" or feature_selection_config == "cell_ranger":
            from modules.featsel import highly_variable_genes
            return lambda adata, n_top_genes=self.n_top_genes: highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=feature_selection_config)
        else:
            return None

    def _get_dim_reduction_method(self, dim_reduction_config):
        if dim_reduction_config is None:
            return None
        elif dim_reduction_config == "pca":
            return lambda adata: sc.pp.pca(adata, n_comps=self.config.get("pca_n_comps", 50))
        elif dim_reduction_config == "glmpca":
            from modules.dimred import glmpca
            return lambda adata: glmpca(adata, n_comps=self.config.get("pca_n_comps", 50))
        else:
            return None
    
    def _get_batch_corr_method(self, batch_corr_config):
        if batch_corr_config is None:
            return None
        elif batch_corr_config == "bbknn":
            from modules.batch_corr import bbknn 
            return lambda adata: bbknn(adata, batch_key=self.config.get("batch_key"))
        elif batch_corr_config == "bbknn_reg":
            from modules.batch_corr import bbknn_reg
            return lambda adata: bbknn_reg(adata
                                           , batch_key=self.config.get("batch_key")
                                           , bbknn_key=self.config.get("bbknn_key")
                                           , confounder_key=self.config.get("confounder_key"))
        else:
            return None
        
    def _get_visualization_method(self, visualization_config):
        if visualization_config is None:
            return None
        elif visualization_config == "umap":
            return lambda adata, n_comp=2: sc.tl.umap(adata, n_components=n_comp, random_state=self.seed)
        # TODO PacMap ??
        elif visualization_config == "trimap":
            return lambda adata, n_comp=2: sc.external.tl.trimap(adata, n_components=n_comp, metric = "cosine")
        else:
            return None
        
    def _get_cycle_genes(self) -> list:
         
        with open('gene_lists/cycle_genes.pkl', 'rb') as f:
            self.cycle_genes = pickle.load(f)
        
        return self.cycle_genes
    
    def _get_pbmc_marker_genes(self) -> list:
         
        return ['IL7R',            # CD4 T cells
                'LYZ', 'CD14',     # CD14+ Monocytes
                'MS4A1',           # B cells
                'CD8A',            # CD8 T cells
                 'NKG7',    # NK cells
                 'CST3',]  # Dendritic Cells ]
    
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
            
            print('feature selection complete')
            
        
        if self.dim_reduction is not None:
            print('running dim reduction...')
            self.dim_reduction(self.adata)
            
        if self.batch_corr is not None:
            adata = self.batch_corr(self.adata)
                
                
        
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
        
        if not ('neighbors' in self.adata.uns.keys()) and (self.batch_corr is None):
            sc.pp.neighbors(self.adata, metric = 'cosine')
            
        n_genes_axis = self.directions_analysis()
         
        self.visualization(self.adata, n_dim)
        print('outputting pca plot...') 
        
        pca_color = 'method'
        
        if not pca_color in self.adata.obs.columns:
            pca_color = 'donor'
        
        
        viz_method = self.config.get('viz')

        for color in ['method', 'donor']:
            
            fig_pca = sc.pl.pca(self.adata, color = color, annotate_var_explained = True, return_fig = True) # plot 2D pca 
            self.save_mpl(fig_pca, title = 'pca-batch-' + color)
            
            
            fig = sc.pl.umap(self.adata, color = color, return_fig = True, show = False) # plot 2D umap
            self.save_mpl(fig, title = 'umap-batch-' + color)
        
        for gene in self._get_pbmc_marker_genes():
            fig_pca_marker = sc.pl.pca(self.adata, color = gene, show = False, return_fig = True)
            self.save_mpl(fig_pca_marker, title = gene + '-marker-pca')
            
            fig_marker = sc.pl.umap(self.adata, color = gene, show = False, return_fig = True)
            self.save_mpl(fig_marker, title = gene + '-marker-umap')
        
        
        if False:
            
            X_ = np.array([self.adata.obsm['X_' + viz_method][:, 0], self.adata.obsm['X_' + viz_method][:, 1], self.adata.obsm['X_pca'][:, 0]]).T
            fig_umap_pca = scatter3D(X_, colors = self.adata.obs, title= viz_method + ' & PC1', labels = [viz_method + ' 1', viz_method + ' 2',  'PC1'])

            self.save_plotly(fig_umap_pca, viz_method + '_pca')
            
            fig_pca_3D = pca3D(self.adata, idx = [0, 1, 2])
            
            scatter_trace = fig_pca_3D.data[0]
            
            x, y, z = scatter_trace.x, scatter_trace.y, scatter_trace.z
            
            scale_factor = np.max(np.mean(np.abs([x, y, z]), axis = 1))
            
            n_genes_axis *= scale_factor
            
            fig_pca_3D.add_trace(go.Scatter3d(
                x = [-n_genes_axis[0], n_genes_axis[0]],
                y = [-n_genes_axis[1], n_genes_axis[1]],
                z = [-n_genes_axis[2], n_genes_axis[2]],
                mode = 'lines',
                line = dict(color = 'black', width = 2),
                name = 'n_genes direction'
            ) )
            
            
            self.save_plotly(fig_pca_3D, 'pca3d')  
            
            fig_umap = scatter3D(self.adata.obsm['X_umap'], colors = self.adata.obs, title='UMAP', labels=['UMAP 1', 'UMAP 2', 'UMAP 3'])
            
            self.save_plotly(fig_umap, 'umap',)

        return None
    
    def analysis(self, PCA: bool = False
                 , plot_matrices = False
                 , genes : str = 'highly_variable'
                 , ) -> None:
        
        from plots import set_matplotlib_style
        
        print('setting matplotlib style...')
        set_matplotlib_style()
        
        print('excluding cycle genes...')
        
        cycle_genes = self._get_cycle_genes()
        marker_genes = self._get_pbmc_marker_genes()
        
        i = 0
        if True:
            
            list_hvg_idx = self.adata.var['highly_variable'] & self.adata.var['GeneName'].isin(marker_genes) 
            
            print('number:', len(list_hvg_idx))
            
            print(self.adata.var['GeneName'][list_hvg_idx])
        
        if genes == 'highly_variable':
            if 'highly_variable' in self.adata.var.columns:
                genes_idx = (self.adata.var['highly_variable'] != (self.adata.var['highly_variable'] 
                                                                & self.adata.var['GeneName'].isin(cycle_genes))) 
        elif genes == 'marker':
            genes_idx = self.adata.var['GeneName'].isin(marker_genes)
        else:
            high_var_no_cycle_idx = list(set(self.adata.var['GeneName']) - set(cycle_genes))
            
        relevant_genes = self.adata.var['GeneName'][genes_idx]
        
        
        print('getting dense array...')
        if PCA:
            A = self.adata.obsm['X_pca']@self.adata.varm['PCs'].T 
            A = A[:, genes_idx]
        else:
            A = self.adata[:, genes_idx].X.toarray()
        
        print('computing correlations...')
        correlations = np.corrcoef(A, rowvar=False)
        print('computing covariances...')
        covariances = np.cov(A, rowvar=False)

        # Setting regularization parameter
        alpha = self.adata.uns['pca']['variance'][-1] 
        print('inverting covariance to get coupling...')
        coupling_matrix = np.linalg.pinv(covariances )
                                        #+ alpha * np.eye(covariances.shape[0]), )
        
        if plot_matrices: 
            from utils import reorder_labels_and_matrix
            from plots import heatmap_with_annotations  
            corr_ordered, labels_corr  = reorder_labels_and_matrix(correlations, relevant_genes)
            
            fig = heatmap_with_annotations(corr_ordered, labels_corr)
            
            self.save_plotly(fig, 'correlation_matrix')
            
            log1p = lambda x: np.sign(x) * np.log(1+np.abs(x))

            coupling_ordered, labels_coupling = reorder_labels_and_matrix(log1p(coupling_matrix), relevant_genes )
            
            self.save_plotly(fig, 'coupling_matrix')
            
            fig = heatmap_with_annotations((coupling_ordered), labels_coupling)
            
            covariances_ordered, labels_cov = reorder_labels_and_matrix(covariances, relevant_genes)
            
            fig = heatmap_with_annotations(log1p(covariances_ordered), labels_cov)
            
            self.save_plotly(fig, 'covariance_matrix')
        
        print('plotting variance ratio...')
        fig, ax = plt.subplots(figsize=(16, 9))
        
        twin = ax.twinx()
        
        var_ratio = self.adata.uns['pca']['variance_ratio']
        n_coms = len(var_ratio)
        left_plot = ax.plot(range(1, n_coms + 1), var_ratio, 'o-', label='Variance Ratio', color = 'C2')
        
        right_plot = twin.plot(range(1, n_coms + 1), np.cumsum(var_ratio), 'o-', label='Cumulative Variance Ratio', color = 'C5')
        
        ax.set_ylabel('Variance Ratio')
        twin.set_ylabel('Cumulative Variance Ratio')
        
        ax.set_xlabel('Number of Components')
        
        ax.yaxis.label.set_color(left_plot[0].get_color())
        twin.yaxis.label.set_color(right_plot[0].get_color())
        
        ax.tick_params(axis='y', colors=left_plot[0].get_color())
        twin.tick_params(axis='y', colors=right_plot[0].get_color())
        
        ax.spines["right"].set_edgecolor(right_plot[0].get_color())
        
        ax.grid(color = left_plot[0].get_color(), alpha = 0.2)
        twin.grid(color = right_plot[0].get_color(), alpha = 0.2)
        
        self.save_mpl(fig, title = 'variance_ratio')
        
        from plots import plot_top_k_joints
        print('plotting joints...')
        
        count_non_zeros = not PCA
        aspect_equal = not PCA 
        
        args = {
            'count_zeros' : count_non_zeros,
            'aspect_equal' : aspect_equal,
            'genes' : relevant_genes 
        }

        self.save_mpl(
            plot_top_k_joints(A, coupling_matrix, relevant_genes, title='coupling-matrix', count_zeros = count_non_zeros, aspect_equal = aspect_equal)
        )
        
        self.save_mpl(
            plot_top_k_joints(A, -coupling_matrix, relevant_genes, title = 'negative-coupling-matrix', count_zeros = count_non_zeros, aspect_equal = aspect_equal)
        )
        
        self.save_mpl(
            plot_top_k_joints(A, correlations, relevant_genes, title = 'correlations',  count_zeros = count_non_zeros, aspect_equal = aspect_equal)
        )
        
        self.save_mpl(
            plot_top_k_joints(A, -correlations, relevant_genes, title = 'negative-correlations', count_zeros = count_non_zeros, aspect_equal = aspect_equal)
        )
        
        """
        self.save_mpl(
            plot_top_k_joints(A, covariances, relevant_genes, title = 'covariances',    count_zeros = count_non_zeros, aspect_equal = aspect_equal)
        )
        
        self.save_mpl(
            plot_top_k_joints(A, -covariances, relevant_genes, title = 'negative-covariances', count_zeros = count_non_zeros, aspect_equal = aspect_equal)
        )
        
        from plots import plot_k_random_joints
        
        self.save_mpl(
            plot_k_random_joints(A, coupling_matrix, relevant_genes)
        )
        """
        from plots import joint_distribution
        
        for i in range(len(marker_genes)):
            for j in range(i+1, len(marker_genes)):
                
                self.save_mpl(
                    joint_distribution( marker_genes[i], marker_genes[j], adata=self.adata,show_nonzero_counts= not PCA),
                    title = 'joint-distribution' + marker_genes[i] + '-' + marker_genes[j], 
                )
        
        return None
    
    def grid_direction_analysis(self, key = 'donor'):
        from plots import heatmap_with_annotations 
        values = self.adata.obs[key].unique()
        k = len(values)
        A = np.zeros((k, k))
        
        for i in range(k):
            for j in range(i+1, k):
                subadata = self.adata[self.adata.obs[key].isin([values[i], values[j]])]
                _, _, pca_share = self.directions_analysis(key, subadata)
                A[i, j] = A[j, i] = pca_share
                
                
        fig = heatmap_with_annotations(A, values)
        self.save_plotly(fig, 'grid_analysis-' + key)
        
        
        return A
        
    
    def directions_analysis(self, v = 'n_genes', adata = None):
        from analysis import correlations_along_vector
        from plots import set_matplotlib_style
        
        if adata is None:
            adata = self.adata
        
        set_matplotlib_style()
        
        if v == 'donor' or v == 'method':
            direction_of_interest = pd.factorize(adata.obs[v])[0]
        elif v in adata.obs.columns:
            direction_of_interest = adata.obs[v]
        
        axis, pearson, spearman = correlations_along_vector(adata, direction_of_interest)

        fig, ax = plt.subplots(figsize=(16, 9))
        
        pearson_corr, pearson_pval = pearson
        
        k = len(pearson_corr)
        
        twin = ax.twinx()
        
        col_left = (235 / 255, 106 / 255, 96/ 255)
        col_right = (0 / 255, 79 / 255, 115 / 255)
        
        left_plot = ax.plot(range(k), np.abs(pearson_corr), 'o-', label = 'pearson', color = col_left)
        right_plot = twin.plot(range(k), np.sqrt(np.cumsum(pearson_corr**2)), 'o-', color = col_right)
        
        ax.set_ylabel('correlation')
        twin.set_ylabel('cumulative correlation')
        
        ax.set_xlabel('k')
        ax.legend(facecolor = 'w')
        
        ax.yaxis.label.set_color(left_plot[0].get_color())
        twin.yaxis.label.set_color(right_plot[0].get_color())
        
        ax.spines["right"].set_edgecolor(right_plot[0].get_color())
        
        ax.grid(color = col_left, alpha = 0.2)
        twin.grid(color = col_right, alpha = 0.2)
        
        
        plt.legend(facecolor='white')
        
        self.save_mpl(fig, title = v + '-directions')
        
        fig, ax = plt.subplots(figsize=(16, 9))
        
        sns.distplot(pearson_pval, label = 'pearson', color = col_left)
        
        plt.yscale('log') 
        
        plt.legend(facecolor = 'w')
        
        self.save_mpl(fig, title = v + '-directions-pvals')
        
        pca_share = np.sum(pearson_corr**2 * adata.uns['pca']['variance_ratio'])
        
        return axis, pearson, pca_share
        
        
    
    def save_mpl(self, fig : plt.Figure, format : str = 'png', title : str = '') -> None:
        
        if title == '':
            title = 'plot' + datetime.now().strftime("%Y%m%d-%H%M%S")
        
        title = title + fig.get_suptitle()
        
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