import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go


def joint_distribution(gene_i: str, gene_j: str, adata: sc.AnnData) -> matplotlib.figure.Figure:
    """
    Plot the joint distribution of two genes.

    Parameters
    ----------
    gene_i : str
        The name of the first gene.
    gene_j : str
        The name of the second gene.
    adata : AnnData
        The annotated data object.

    Returns
    -------
    matplotlib.figure.Figure
        The figure with the joint distribution plot.
    """
    # Get the indices of the genes in adata.var_names
    gidx1, gidx2 = adata.var_names.get_loc(gene_i), adata.var_names.get_loc(gene_j)
    
    # Create a figure with a 16x16 inch size
    fig, ax = plt.subplots(figsize=(16, 16))
    
    # Extract the expression values of the genes as numpy arrays
    exp1, exp2 = adata.X[:, gidx1].A1, adata.X[:, gidx2].A1
    
    # Plot the scatter plot of the expression values
    plt.scatter(exp1, exp2, c='k', s=1, alpha=0.5)
    
    # Add labels to the x and y axes
    plt.xlabel(gene_i + ' (non-zero counts ' + str(np.sum(exp1 != 0) / adata.X.shape[0]) + ')')
    plt.ylabel(gene_j + ' (non-zero counts ' + str(np.sum(exp2 != 0) / adata.X.shape[0]) + ')')
    
    # Add a title to the plot
    title_str = '(' + str((exp1 > 0 & exp2 > 0).sum() / adata.X.shape[0]) + ')'
    plt.title(gene_i + ' vs ' + gene_j + ' (both non-zero counts ' + title_str + ')')
    
    # Add a text label to the plot
    plt.text(-1., -1., '0-0 ' + str(np.sum((exp1 + exp2) == 0) / adata.X.shape[0]), fontsize=15)

    # Return the figure
    return fig
    
def pca_3d(adata: sc.AnnData, color_key : str = 'method') -> go.Figure:
        
    x,y,z = adata.obsm['X_pca'][:,0], adata.obsm['X_pca'][:,1], adata.obsm['X_pca'][:,2]
    variance_ratio = adata.uns['pca']['variance_ratio']
     
    labels = ['PC' + str(i) + "{:10.4f}".format(variance_ratio[i]) for i in range(1, 4)]
    
    fig = px.scatter_3d(x=x, y=y, z=z,
              color=adata.obs[color_key])
              
    fig.update_layout(xaxis_title=labels[0] , yaxis_title=labels[1], zaxis_title=labels[2])
    
    fig.update_traces(marker_size = 1)
    fig.show()
    
    
    return fig 