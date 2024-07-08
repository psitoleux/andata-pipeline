import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

from plotutils import *

def set_matplotlib_style():
    plt.style.use("ggplot")
    
    colors = get_tol()

    matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(
        color= colors)

def get_tol():
    
    return ['#332288', '#117733', '#44AA99', '#88CCEE'
              , '#DDCC77', '#CC6677', '#AA4499', '#882255']

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
