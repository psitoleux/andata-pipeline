import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

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
    # Extract x, y, z coordinates
    
    df = extract_coordinates(x, colors, idx, dim = 3)
    
    categorical_cols, continuous_cols = determine_variable_types(df)


    # Define color options dictionary
    color_options = create_color_options(df, categorical_cols, continuous_cols)

    # Generate color palette for categorical variable
    category_colors = px.colors.qualitative.Set1

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

def determine_variable_types(df: pd.DataFrame):
    categorical_cols = []
    continuous_cols = []
    
    for col in df.columns:
        ctype = df[col].dtype
        if isinstance(ctype, pd.CategoricalDtype) or isinstance(ctype, pd.StringDtype):
            categorical_cols.append(col)
        else:
            continuous_cols.append(col)
    
    return categorical_cols, continuous_cols

def create_color_options(df: pd.DataFrame, categorical_cols: list, continuous_cols: list):
    color_options = {}
    
    for col in categorical_cols:
        color_options[f'{col} (Categorical)'] = col
    
    for col in continuous_cols:
        color_options[f'{col} (Continuous)'] = col
    
    return color_options


def create_menu_buttons(fig: go.Figure, color_options: dict, categorical_cols: list, df : pd.DataFrame):
    menu_buttons = []
    index = 0
    for label, col in color_options.items():
        visibility = [False] * len(fig.data)
        if col in categorical_cols:  # Categorical
            for _ in df[col].unique():
                visibility[index] = True
                index += 1
        else:  # Continuous
            visibility[index] = True
            index += 1
        
        menu_buttons.append(
            dict(
                args=[{"visible": visibility}],
                label=label,
                method="update"
            )
        )       

    return menu_buttons

def extract_coordinates(x: np.ndarray, df : pd.DataFrame, idx: list, dim: int):
    df['x'], df['y'] = x[:, idx[0]], x[:, idx[1]]
    if dim == 3:
        df['z'] = x[:, idx[2]]
    return df

def add_traces(fig: go.Figure, df: pd.DataFrame, color_options: dict, category_colors: list, dim: int):
    for label, col in color_options.items():
        if col in df.columns:
            if isinstance(df[col].dtype, (pd.CategoricalDtype, pd.StringDtype)):
                for i, cat in enumerate(df[col].unique()):
                    if dim == 2:
                        fig.add_trace(
                            go.Scatter(
                                x=df[df[col] == cat]['x'],
                                y=df[df[col] == cat]['y'],
                                mode='markers',
                                marker=dict(color=category_colors[i % len(category_colors)]),
                                name=f"{label}: {cat}",
                                visible=False
                            )
                        )
                    elif dim == 3:
                        fig.add_trace(
                            go.Scatter3d(
                                x=df[df[col] == cat]['x'],
                                y=df[df[col] == cat]['y'],
                                z=df[df[col] == cat]['z'],
                                mode='markers',
                                marker=dict(color=category_colors[i % len(category_colors)]),
                                name=f"{label}: {cat}",
                                visible=False
                            )
                        )
            else:
                if dim == 2:
                    fig.add_trace(
                        go.Scatter(
                            x=df['x'],
                            y=df['y'],
                            mode='markers',
                            marker=dict(
                                color=df[col],
                                colorbar=dict(title=label),
                                colorscale='Viridis'
                            ),
                            name=label,
                            visible=False
                        )
                    )
                elif dim == 3:
                    fig.add_trace(
                        go.Scatter3d(
                            x=df['x'],
                            y=df['y'],
                            z=df['z'],
                            mode='markers',
                            marker=dict(
                                color=df[col],
                                colorbar=dict(title=label),
                                colorscale='Viridis'
                            ),
                            name=label,
                            visible=False
                        )
                    )
    return fig


def update_layout_menu(fig: go.Figure, menu_buttons: list):
   fig.update_layout(
        updatemenus=[
            dict(
                active=0,
                buttons=menu_buttons,
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.1,
                xanchor="left",
                y=1.15,
                yanchor="top"
            )], 
    )
   
   return fig
