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

    fig.update_layout( scene = dict(xaxis_title=labels[0] 
                                    , yaxis_title=labels[1]
                                    , zaxis_title=labels[2]))
                      

    fig.update_traces(marker_size = 1)
    fig.show()
    
    scene = dict(
                    xaxis_title='X AXIS TITLE',
                    yaxis_title='Y AXIS TITLE',
                    zaxis_title='Z AXIS TITLE'),

    return fig 

def plot3D(x: np.ndarray, colors: pd.DataFrame, idx: np.ndarray = [0, 1, 2], labels: list = ['X', 'Y', 'Z']) -> go.Figure:
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
    x, y, z = x[:, idx[0]], x[:, idx[1]], x[:, idx[2]]

    # Determine categorical and continuous variables
    categorical_cols = []
    continuous_cols = []

    # Identify categorical and continuous columns
    for col in df.columns:
        ctype = df[col].dtype
        if isinstance(ctype, pd.CategoricalDtype) or isinstance(ctype, pd.StringDtype):
            categorical_cols.append(col)
        elif pd.api.types.is_numeric_dtype(df[col]):
            continuous_cols.append(col)

    # Define color options dictionary
    color_options = {}

    # Add categorical color options
    for col in categorical_cols:
        color_options[f'{col} (Categorical)'] = col

    # Add continuous color options
    for col in continuous_cols:
        color_options[f'{col} (Continuous)'] = col

    # Generate color palette for categorical variable
    category_colors = px.colors.qualitative.Set1

    # Create a new 3D scatter plot
    fig = go.Figure()

    # Add traces for each color option
    for label, col in color_options.items():
        if col in categorical_cols:  # Categorical
            for i, cat in enumerate(df[col].unique()):
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
        else:  # Continuous
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

    # Make the first trace visible
    fig.data[0].visible = True

    # Create the menu buttons
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

    # Add dropdown menu to the plot
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
            scene=dict(xaxis_title=labels[0], yaxis_title=labels[1], zaxis_title=labels[2])
    )

    # Display the plot
    fig.show()
    
    
    return fig 
