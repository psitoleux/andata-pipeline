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

