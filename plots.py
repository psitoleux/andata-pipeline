import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

from plotutils import *

def set_matplotlib_style(name = "ggplot"):
    plt.style.use(name)
    
    colors = get_tol()

    matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(
        color= colors)

def get_tol():
    
    return ['#332288', '#117733', '#44AA99', '#88CCEE'
              , '#DDCC77', '#CC6677', '#AA4499', '#882255']

def joint_distribution(gene_i: str, gene_j: str
                       , adata: sc.AnnData
                       , color : str = 'k'
                       , PCA : bool = False) -> matplotlib.figure.Figure:
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
    show_nonzero_counts : bool
        Whether to show the number of non-zero counts in the labels (default: True)

    Returns
    -------
    matplotlib.figure.Figure
        The figure with the joint distribution plot.
    """
    # Get the indices of the genes in adata.var_names
    gidx1, gidx2 = adata.var_names.get_loc(gene_i), adata.var_names.get_loc(gene_j)
    
    # Create a figure with a 16x16 inch size
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    
    # Extract the expression values of the genes as numpy arrays
    if PCA:
        idx_i = adata.var['highly_variable'][adata.var['highly_variable'] == True].index.get_loc(gene_i)
        idx_j = adata.var['highly_variable'][adata.var['highly_variable'] == True].index.get_loc(gene_j)
        
        PC_projection = adata.varm['PCs'][[idx_i, idx_j], :]
        
        recon_expr = adata.obsm['X_pca']@PC_projection.T
        
        expr1, expr2 = recon_expr[:, 0], recon_expr[:, 1]
    else:
        expr1, expr2 = adata.X[:, gidx1].todense(), adata.X[:, gidx2].todense()
    
    
    show_nonzero_counts = False
    
    ax.plot(expr1,expr2, 'o', markersize=2, alpha = 0.01,
                label = gene_i + ' - ' + gene_j,
                color = color) 
 
    # Add labels to the x and y axes
    if show_nonzero_counts:
        plt.xlabel(gene_i + ' (non-zero counts ' + str(np.sum(expr1 != 0) / adata.X.shape[0]) + ')')
        plt.ylabel(gene_j + ' (non-zero counts ' + str(np.sum(expr2 != 0) / adata.X.shape[0]) + ')')
    else:
        plt.xlabel(gene_i)
        plt.ylabel(gene_j)
    
    # Add a title to the plot
    if show_nonzero_counts:
        title_str = '(' + str((expr1 > 0 & expr2 > 0).sum() / adata.X.shape[0]) + ')'
        plt.title(gene_i + ' vs ' + gene_j + ' (both non-zero counts ' + title_str + ')')
    else:
        plt.title(gene_i + ' vs ' + gene_j)
    
    # Add a text label to the plot
    if show_nonzero_counts:
        plt.text(-1., -1., '0-0 ' + str(np.sum((expr1 + expr2) == 0) / adata.X.shape[0]), fontsize=15)

    # Return the figure
    return fig
    
def pca3D(adata: sc.AnnData, idx: list = [0, 1, 2]) -> go.Figure:
    """
    Creates a 3D PCA scatter plot using Plotly.

    Parameters:
    adata (sc.AnnData): Annotated data containing PCA information.
    color_key (str): Key in adata.obs to use for coloring.
    idx (list): Indices of PCA components to plot (default: [0, 1, 2])

    Returns:
    go.Figure: Plotly 3D scatter plot figure.
    """

    # Extract PCA coordinates
    x, y, z = adata.obsm['X_pca'][:, idx[0]], adata.obsm['X_pca'][:, idx[1]], adata.obsm['X_pca'][:, idx[2]]
    
    # Extract variance ratios for axes labels
    variance_ratio = adata.uns['pca']['variance_ratio']
    labels = [f'PC{i+1} ({variance_ratio[i]:.4f})' for i in idx]

    # Call scatter3D function to generate the plot
    fig = scatter3D(x=np.array([x, y, z]).T, colors=adata.obs, title='PCA 3D Plot', labels=labels)

    return fig


def scatter2D(x: np.ndarray
           , colors: pd.DataFrame
           , title: str
           , idx: np.ndarray = [0, 1]
           , labels = ['X', 'Y']) -> go.Figure:
    """
    This function creates a 2D scatter plot based on the input data.

    Parameters:
    x (np.ndarray): Input data for the plot
    idx (np.ndarray): Indices for selecting the columns
    labels (list): Labels for x, y axes
    colors (pd.DataFrame): Color data for the plot

    Returns:
    go.Figure: 2D scatter plot figure
    """
    
    categorical_cols, continuous_cols = determine_variable_types(colors)

    # Define color options dictionary
    color_options = create_color_options(colors, categorical_cols, continuous_cols)

    # Generate color palette for categorical variable
    category_colors = get_tol()

    # Extract x, y coordinates
    df = extract_coordinates(x, colors, idx, dim = 2)
    

    # Create a new 2D scatter plot
    fig = go.Figure()

    fig = add_traces(fig, df, color_options, category_colors, dim=2)
    
    fig.data[0].visible = True
    
    # Create menu buttons
    
    menu_buttons = create_menu_buttons(fig, color_option, categorical_cols, df)
    
    fig = update_layout_menus(fig, menu_buttons)
    
    # Set title
    fig.update_layout(title=title, xaxis_title=labels[0], yaxis_title=labels[1])
    
    # Axes labels
    
    fig.update_layout(xaxis_title=labels[0], yaxis_title=labels[1])
    
    fig.update_trace(marker=dict(size=1))
    
    fig.update_layout(title=title)
    
    fig.show()
    
    return fig
    

def scatter3D(x: np.ndarray
           , colors: pd.DataFrame
           , title: str
           , idx: np.ndarray = [0, 1, 2]
           , labels = ['X', 'Y', 'Z']) -> go.Figure:
    """
    This function creates a 3D scatter plot based on the input data.

    Parameters:
    x (np.ndarray): Input data for the plot
    idx (np.ndarray): Indices for selecting the columns
    labels (list): Labels for x, y, and z axes
    colors (pd.DataFrame): Color data for the plot

    Returns:
    go.Figure: 3D scatter plot figure
    """
    
    categorical_cols, continuous_cols = determine_variable_types(colors)


    # Define color options dictionary
    color_options = create_color_options(colors, categorical_cols, continuous_cols)

    # Generate color palette for categorical variable
    category_colors = get_tol()

    # Extract x, y, z coordinates
    
    df = extract_coordinates(x, colors, idx, dim = 3)
    

    # Create a new 3D scatter plot
    fig = go.Figure()

    # Add traces for each color option
    fig = add_traces(fig, df, color_options, category_colors, dim = 3)
    
    # Make the first trace visible
    fig.data[0].visible = True

    # Create the menu buttons
    menu_buttons = create_menu_buttons(fig, color_options, categorical_cols, df)

    # Add dropdown menu to the plot
    
    fig = update_layout_menu(fig, menu_buttons)
    
    # Axes labels
    
    fig.update_layout(scene=dict(xaxis_title=labels[0], yaxis_title=labels[1], zaxis_title=labels[2]))
    
     # Small markers 
    
    fig.update_traces(marker=dict(size=1))

    # Set plot title
    fig.update_layout(title=title)
     
    # Display the plot
    fig.show()
    
    
    return fig 




def heatmap_with_annotations(matrix: np.ndarray, labels: list):


    # Create a list of gene name annotations
    annotations = [
        [dict(text=[i], x=i, y=j, xref='x1', yref='y1', showarrow=False) 
        for i in range(len(labels))] 
        for j in range(len(labels))
    ]

    # Create the heatmap figure
    fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=labels,
            y=labels,
            colorscale='Viridis',
            colorbar=dict(title='Correlation'),
            hoverongaps=False,
            hovertemplate='Gene A: %{y}<br>Gene B: %{x}<br>Correlation: %{z}<extra></extra>'
        ))

    # Update layout
    fig.update_layout(
        title='Gene Expression Correlation Matrix',
        xaxis_title='Genes',
        yaxis_title='Genes',
        width=1600,
        height=1600,
        xaxis=dict(side='top')
    )

    # Show the plot
    fig.show()
    
    return fig

def plot_top_k_joints(data : np.ndarray, matrix : np.ndarray, genes
                        , title : str = ''
                        , k :int = 9
                        , count_zeros : bool = False
                        , aspect_equal : bool = False
                        , include_value : bool = False
                        , matrix_name : str = ''):
    
    top_k_indices = top_k_off_diagonal_indices_symmetric(matrix, k)
    
    return plot_k_joints(data, matrix, genes, top_k_indices, title, k, count_zeros, aspect_equal, include_value)
    

def plot_2_joints(data1 : np.ndarray
                  , data2 : np.ndarray
                  , matrix1 : np.ndarray
                  , matrix2 : np.ndarray
                  , genes : pd.Series
                  , idx : np.array
                  , title : str = ''
                  , matrix_name1 : str = ''
                  , matrix_name2 : str = ''
                  , include_value : bool = True
                  , color : str = 'C0'):
    
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    
    genex, geney = genes.iloc[idx[0]], genes.iloc[idx[1]]

    
    for i, data, matrix, name,ax in zip(range(2), [data1, data2], [matrix1, matrix2] ,[matrix_name1, matrix_name2], ax.flat):
        
        sub_title = name
        
        if include_value:
            sub_title =  '{0:.2f}'.format(matrix[idx[0], idx[1]])
        xlabel, ylabel = genex, geney
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(sub_title)
        
        if i == 0:
            alpha = 0.2
        elif i == 1:
            alpha = 0.01
        
        ax.plot(data[:, idx[0]], data[:, idx[1]], '.', color = color, alpha = alpha, markersize = 1)
    
    plt.tight_layout()
    fig.suptitle(title)
    
    return fig

def plot_k_joints(data : np.ndarray
                    , matrix : np.ndarray
                    , genes : pd.Series
                    , indices : np.ndarray
                    , title : str = ''
                    , k : int = 9
                    , count_zeros : bool = False
                    , aspect_equal : bool = False
                    , include_value : bool = False
                    , matrix_name : str = ''):
        
    
    l = int(np.sqrt(k))
    

    fig, axs = plt.subplots(l, l, figsize=(l*4, l*4))

    for i,ax in zip(range(k), axs.flat):
        
        genex, geney = genes.iloc[indices[i, 0]], genes.iloc[indices[i, 1]]

        sub_title = matrix_name
        if include_value:
            sub_title =  '\n value ' + '{0:.2f}'.format(matrix[indices[i, 0], indices[i, 1]])
        xlabel, ylabel = genex, geney
        
        alpha = 0.01
        if count_zeros:
            share_both_zeros = np.sum((data[:, indices[i, 0]] == 0) & (data[:, indices[i, 1]] == 0)) / data.shape[0]
            nb_non_zeros_both = np.sum((data[:, indices[i, 0]] != 0) & (data[:, indices[i, 1]] != 0)) 
            nb_non_zeros_x = np.sum(data[:, indices[i, 0]] != 0) 
            nb_non_zeros_y = np.sum(data[:, indices[i, 1]] != 0)
            
            sub_title =  ' \n share both zeros ' + '{0:.4f}'.format(share_both_zeros) + \
                ' \n nb both non-zeros ' + str(nb_non_zeros_both) + sub_title
            xlabel = xlabel + ' - ' + str(nb_non_zeros_x) + ' non-zeros'
            ylabel = ylabel + ' - ' + str(nb_non_zeros_y) + ' non zeros'
            
            alpha = 0.5
            if nb_non_zeros_both > 1000:
                alpha = 0.1
            if nb_non_zeros_both > 10000:
                alpha = 0.01
            
             
        
        ax.plot(data[:, indices[i, 0]]
                ,data[:, indices[i, 1]], 'o', markersize=2, alpha = alpha,
                label = genex + ' - ' + geney,
                color = 'C' + str(i % 10)) 

        ax.set_title(sub_title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        if aspect_equal:
            ax.set_aspect('equal')
        
    fig.suptitle(title) 
    plt.tight_layout()

    return fig


def plot_k_random_joints(data : np.ndarray, matrix : np.ndarray, genes, k = 9):
        
    indices = random_offdiag_idx(data.shape[1], k)
        
    return plot_k_joints(data, matrix, genes, indices, title = 'random-pairs', k = k)


def plot_compare_symmetric_matrices(data1 : np.ndarray
                                    , data2 : np.ndarray
                                    , name1 : str
                                    , name2 : str
                                    , mtx_kind : str
                                    , scale = 'linear'
                                    ):
    
    N = data1.shape[1]
    idx = np.triu_indices(N, 1)
    
    
    fig, ax = plt.subplots(figsize = (9, 9))
    
    
    if scale == 'linear':
        ax.plot(data1[idx], data2[idx], 'o', alpha = 0.01)
    elif scale == 'log':
        pass
    
    ax.set_xlabel(name1)
    ax.set_ylabel(name2)
    
    ax.set_title(mtx_kind)
    
    return fig
    
    
    
    
    
    