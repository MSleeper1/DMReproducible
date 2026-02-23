import pandas as pd
import sys
import plotly.express as px
from scipy.cluster.hierarchy import linkage, leaves_list
import plotly.graph_objects as go
import numpy as np


# Reads in csv to df with cell type as index
def prep_uxm_df(uxm_file):
    # Read file to df
    df = pd.read_csv(uxm_file)
    # Set CellType as index
    df.set_index("CellType", inplace=True)
    return df


#--- Clustered heatmap ---#

def plotly_hm_clustered(df, plot_name):

    # Transpose to have samples on the x-axis and cell types on y-axis
    df_t = df.T
    # Perform hierarchical clustering on rows (samples)
    row_linkage = linkage(df_t, method='average', metric='euclidean')
    row_order = leaves_list(row_linkage)

    # Perform hierarchical clustering on columns (cell types)
    col_linkage = linkage(df_t.T, method='average', metric='euclidean')
    col_order = leaves_list(col_linkage)

    # Reorder rows and columns based on clustering
    df_clustered = df_t.iloc[row_order, col_order]

    # Plot interactive heatmap with clustered order
    fig = px.imshow(
        df_clustered,
        labels=dict(x="Cell Type", y="Sample", color="Contribution"),
        x=df_clustered.columns,
        y=df_clustered.index,
        color_continuous_scale='Viridis',
        aspect="auto",
    )

    fig.update_layout(
        title="Clustered Cellular Contributions Heatmap",
        xaxis_tickangle=45,
        width=1200,
        height=800,
    )

    # Save to HTML
    fig.write_html(plot_name)

def heatmap_complex(df, plot_name):
    # Transpose: samples as rows, cell types as columns
    df_t = df.T

    # Optional: filter cell types by max contribution to reduce options
    threshold = 0.01
    filtered_cols = df_t.columns[df_t.max(axis=0) > threshold]
    df_t_filtered = df_t[filtered_cols]

    # Cluster rows and columns
    row_linkage = linkage(df_t_filtered, method='average', metric='euclidean')
    row_order = leaves_list(row_linkage)
    col_linkage = linkage(df_t_filtered.T, method='average', metric='euclidean')
    col_order = leaves_list(col_linkage)

    # Reorder dataframe
    df_clustered = df_t_filtered.iloc[row_order, col_order]

    # Prepare list of cell types for dropdown (show all individually + an option for all)
    cell_types = list(df_clustered.columns)
    dropdown_options = [{'label': 'All Cell Types', 'value': 'all'}] + \
                    [{'label': ct, 'value': ct} for ct in cell_types]

    # Initialize figure with all cell types
    fig = go.Figure()

    heatmap = go.Heatmap(
        z=df_clustered.values,
        x=df_clustered.columns,
        y=df_clustered.index,
        colorscale='Viridis',
        colorbar=dict(title="Contribution")
    )
    fig.add_trace(heatmap)

    # Update layout with dropdown to filter columns (cell types)
    fig.update_layout(
        title="Clustered Cellular Contributions Heatmap with Cell Type Filter",
        xaxis=dict(tickangle=45),
        width=1200,
        height=800,
        updatemenus=[
            dict(
                active=0,
                buttons=[
                    dict(
                        label=option['label'],
                        method="update",
                        args=[{"z": [df_clustered[([option['value']] if option['value'] != 'all' else cell_types)].values],
                            "x": [df_clustered[([option['value']] if option['value'] != 'all' else cell_types)].columns]},
                            {"xaxis": {"tickangle": 45}}]
                    )
                    for option in dropdown_options
                ],
                direction="down",
                showactive=True,
                x=0,
                y=1.15,
                xanchor='left',
                yanchor='top'
            )
        ]
    )

    fig.write_html(plot_name)

### MAIN SCRIPT BEGINS ###

if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            uxm_csv = str(sys.argv[1])
        except ValueError:
            print("Argument is not a valid srting.")
        
        # create df for plotting
        df = prep_uxm_df(uxm_csv)

        # # Sort cell types (rows) by average contribution, descending
        # df_sorted = df.loc[df.mean(axis=1).sort_values(ascending=False).index]
        # # Filter out rows where the sum across all columns (samples) is zero
        # df_bar_filtered = df_sorted.loc[df_sorted.sum(axis=1) > 0]
        # Remove rows and columns with zero variance for normalization
        df_clean = df.loc[df.var(axis=1) > 0, df.var(axis=0) > 0]

        # plotting interactive clustered heatmap
        clustered_heatmap_name = uxm_csv.replace(".csv", ".heatmap_clustered.html")
        plotly_hm_clustered(df_clean, clustered_heatmap_name)

        complex_heatmap_name = uxm_csv.replace(".csv", ".heatmap_complex.html")
        heatmap_complex(df, complex_heatmap_name)
        
    else:
        print("No UXM file provided.")

