import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
import sys
import numpy as np
from scipy.spatial.distance import pdist
from PIL import Image 

#-----------------------#
#   Functions Defined   #
#-----------------------#


#--- Reading, checking, and cleaning data ---#

# Reads in csv to df with cell type as index
def prep_uxm_df(uxm_file):
    # Read file to df
    df = pd.read_csv(uxm_file)
    # Set CellType as index
    df.set_index("CellType", inplace=True)
    return df

# Checks data to see if issues exist before normalizing data
def check_data(df):
        # Check for NaNs or infinite values:
        print("Any NaNs?", df.isnull().values.any())
        print("Any infinite values?", np.isinf(df.values).any())
        print("Rows with zero variance:", (df.var(axis=1) == 0).sum())
        print("Columns with zero variance:", (df.var(axis=0) == 0).sum())
        print("Data shape:", df.shape)
        print(df.dtypes)
        # Compute pairwise distances manually
        distances = pdist(df_clean.values, metric="euclidean")
        print("Distances:", distances)
        print("Any non-finite in distances?", np.isfinite(distances).all())

#--- Stacked Barcharts ---#

# Plots cellular contribution as stacked bar plot
def plot_uxm_stacked(df, plot_name):
    # plot
    fig, ax = plt.subplots(figsize=(14, 8))
    # Transpose so samples are on x-axis
    df.T.plot(
        kind="bar",
        stacked=True,
        figsize=(14, 8),
        ax=ax,
        colormap="tab20"  # nice categorical colormap
    )
    # Labels & Formatting
    ax.set_ylabel("Estimated Cellular Contribution")
    ax.set_xlabel("Sample")
    ax.set_title("Estimated Cellular Composition of Tissue Samples", fontsize=16)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Cell Type")
    plt.tight_layout()
    # Save
    plt.savefig(plot_name, dpi=300)


#--- Heatmaps ---#

# Plots cellular contribution as heatmap
def plot_heatmap(df, plot_name):
    # Create the heatmap
    plt.figure(figsize=(14, 10))
    sns.heatmap(
        df,
        cmap="viridis",      # color palette (can try "mako", "plasma", "coolwarm")
        cbar_kws={'label': 'Estimated Cellular Contribution'},
        linewidths=0.3
    )
    # Labels
    plt.title("Cell Type Contribution Heatmap", fontsize=16)
    plt.xlabel("Sample")
    plt.ylabel("Cell Type")
    plt.tight_layout()
    # Save
    plt.savefig(plot_name, dpi=300)

# Plots cellular contribution as clustered heatmap
def plot_clustered_heatmap(df, plot_name):
    # Create the clustered heatmap
    sns.clustermap(
        df,
        cmap="viridis",       # try "mako", "coolwarm", "plasma" for different looks
        figsize=(14, 12),
        cbar_kws={'label': 'Estimated Cellular Contribution'},
        linewidths=0.3,
        metric="euclidean",   # distance metric for clustering
        method="average"      # linkage method (try "ward", "complete", "single")
    )
    plt.suptitle("Clustered Heatmap of Cell Type Contributions", y=1.02, fontsize=16)
    # Save
    plt.savefig(plot_name, dpi=300, bbox_inches='tight')


#--- Paneled images with multiple plots ---#

# Plots heatmap and bar chart (top 10 cell types)
def plots_joined_reduced(df, plot_name):

    TOP_N = 10  # Number of top cell types to keep in the bar chart legend
    
    # Keep original column order (matches stacked bar chart)
    sample_order = df.columns.tolist()

    # Create "Other" category for stacked bar chart
    top_celltypes = df.index[:TOP_N]
    df_bar = df.copy()
    
    # Sum all contributions not in top N into "Other"
    df_bar.loc["Other"] = df_bar.loc[~df_bar.index.isin(top_celltypes)].sum()
    
    # Keep only top N + Other for stacked bar chart
    df_bar = df_bar.loc[list(top_celltypes) + ["Other"]]
    
    # Create combined figure
    fig = plt.figure(figsize=(22, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1])
    
    # Left Panel: Stacked Bar Chart
    ax1 = plt.subplot(gs[0])
    df_bar.T.plot(
        kind="bar",
        stacked=True,
        ax=ax1,
        colormap="tab20"
    )
    ax1.set_ylabel("Estimated Cellular Contribution")
    ax1.set_xlabel("Sample")
    ax1.set_title(f"Cellular Composition by Sample (Top {TOP_N} + Other)", fontsize=16)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Cell Type")
    
    # Right Panel: Full Heatmap
    ax2 = plt.subplot(gs[1])
    sns.heatmap(
        df[sample_order],  # full set of cell types for heatmap
        cmap="viridis",
        cbar_kws={'label': 'Estimated Cellular Contribution'},
        linewidths=0.3,
        ax=ax2
    )
    ax2.set_title("Heatmap of Cell Type Contributions", fontsize=16)
    ax2.set_xlabel("Sample")
    ax2.set_ylabel("Cell Type")
    # Save
    plt.tight_layout()
    plt.savefig(plot_name, dpi=300, bbox_inches='tight')


# plots normalized and clustered heatmap
def clustered_hm_normalized(df, plot_name):
    # Normalize rows to 0–1
    df_normalized = df.sub(df.min(axis=1), axis=0)
    df_normalized = df_normalized.div(df_normalized.max(axis=1), axis=0)
    sns.clustermap(
        df_normalized,
        cmap="viridis",
        metric="euclidean",
        method="average",
        linewidths=0.3,
        figsize=(8, 10),
        cbar_kws={'label': 'Normalized Contribution (0–1)'}
    )
    plt.suptitle("Clustered Heatmap of Cell Type Contributions (Normalized)", y=1.02, fontsize=16)
    # Save
    plt.savefig(plot_name, dpi=300, bbox_inches='tight')


# plots heatmaps of normalized and non normalized data side by side
def paneled_hm_normalized(df, plot_name):
    # Normalize
    df_normalized = df.sub(df.min(axis=1), axis=0)
    df_normalized = df_normalized.div(df_normalized.max(axis=1), axis=0)
    # Create and save original cluster heatmap
    cluster_orig = sns.clustermap(
        df,
        cmap="viridis",
        metric="euclidean",
        method="average",
        linewidths=0.3,
        figsize=(8, 10),
        cbar_kws={'label': 'Contribution'}
    )
    cluster_orig.savefig("cluster_orig.png")
    plt.close(cluster_orig.figure)

    # Create and save normalized cluster heatmap
    cluster_norm = sns.clustermap(
        df_normalized,
        cmap="viridis",
        metric="euclidean",
        method="average",
        linewidths=0.3,
        figsize=(8, 10),
        cbar_kws={'label': 'Normalized Contribution (0–1)'}
    )
    cluster_norm.savefig("cluster_norm.png")
    plt.close(cluster_norm.figure)

    # Load images
    img_orig = Image.open("cluster_orig.png")
    img_norm = Image.open("cluster_norm.png")

    # Plot side-by-side
    fig, axs = plt.subplots(1, 2, figsize=(24, 12))
    axs[0].imshow(img_orig)
    axs[0].axis('off')
    axs[0].set_title("Original Contributions", fontsize=16)

    axs[1].imshow(img_norm)
    axs[1].axis('off')
    axs[1].set_title("Normalized Contributions", fontsize=16)

    plt.tight_layout()
    plt.savefig(plot_name, dpi=300)



### MAIN SCRIPT BEGINS ###

if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            uxm_csv = str(sys.argv[1])
        except ValueError:
            print("Argument is not a valid srting.")
        
        # create df for plotting
        df = prep_uxm_df(uxm_csv)

        # Sort cell types (rows) by average contribution, descending
        df_sorted = df.loc[df.mean(axis=1).sort_values(ascending=False).index]
        # Filter out rows where the sum across all columns (samples) is zero
        df_bar_filtered = df_sorted.loc[df_sorted.sum(axis=1) > 0]
        # Remove rows and columns with zero variance for normalization
        df_clean = df_sorted.loc[df.var(axis=1) > 0, df.var(axis=0) > 0]

        # plotting stacked bar
        bar_plot_name = uxm_csv.replace(".csv", ".stacked_bar.png")
        plot_uxm_stacked(df_bar_filtered, bar_plot_name)

        # plotting heatmap
        heatmap_plot_name = uxm_csv.replace(".csv", ".heatmap.png")
        plot_heatmap(df_sorted, heatmap_plot_name)

        # plotting clustered heatmap
        clustered_plot_name = uxm_csv.replace(".csv", ".clustered_heatmap.png")
        plot_clustered_heatmap(df_sorted, clustered_plot_name)

        # plotting multiple plots with reduced legend to show top 10
        panel_plot_reduced_name = uxm_csv.replace(".csv", ".paneled_heatmap_barplot_reduced_legend.png")
        plots_joined_reduced(df_sorted, panel_plot_reduced_name)
        
        # heatmap with normalized values 0-1
        hm_normalized_name = uxm_csv.replace(".csv", ".heatmap_normalized.png")
        clustered_hm_normalized(df_clean, hm_normalized_name)

        # paneled image with normalized heatmap and not normalized heatmap side by side
        hm_normalized_paneled_name = uxm_csv.replace(".csv", ".paneled_heatmap_normalized.png")
        paneled_hm_normalized(df_clean, hm_normalized_paneled_name)

    else:
        print("No UXM file provided.")
