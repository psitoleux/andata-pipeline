import pandas as pd
import plotly.graph_objects as go
import numpy as np
import plotly.express as px
import plotly.io as pio





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


def random_offdiag_idx(N, k):
    
    indices = [(i,j) for i in range(N) for j in range(i+1, N)]
    indices = np.array(indices)
    iidx = np.random.choice(len(indices), size=k, replace=False).astype(int)
    
    return indices[iidx]
    

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
