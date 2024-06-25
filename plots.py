import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
import matplotlib.pyplot as plt

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
    
    # Extract the expression values of the genes
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
    
    