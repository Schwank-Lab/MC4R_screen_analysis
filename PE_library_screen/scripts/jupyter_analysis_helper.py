# Helper functions for editing summarization notebook

import pandas as pd
from os import listdir
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from scipy.stats import gaussian_kde
import seaborn as sn

# Set the style for the plots
lsize = 14
titlesize = 16
axesthickness = 1.5
plt.rcParams['axes.linewidth'] = axesthickness
pdffonttpye = 42
psfonttype = 42
fontfamily = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = pdffonttpye
mpl.rcParams['ps.fonttype'] = psfonttype
mpl.rcParams['font.family'] = fontfamily


def list_files(directory, pattern='combined'):
    """List all files in a directory containing a specific pattern."""
    return [f for f in listdir(directory) if pattern in f]


def create_scatter_plot(df, x_col, y_col, title, xlabel, ylabel, xlim=100, ylim=100, savepath=None):
    """Create a scatter plot with regression line and correlation values."""
    # Reset to matplotlib defaults first
    plt.rcdefaults()

    # Then reapply your custom styling using the global variables
    plt.rcParams['axes.linewidth'] = axesthickness
    mpl.rcParams['pdf.fonttype'] = pdffonttpye
    mpl.rcParams['ps.fonttype'] = psfonttype
    mpl.rcParams['font.family'] = fontfamily

    print("Full length of df:", len(df))
    df = df[[x_col, y_col]].dropna()
    print("Plotted number after dropping NaNs:", len(df))
    if len(df) == 0:
        raise ValueError("No data to plot after dropping NaNs.")
    x = np.array(df[x_col])
    y = np.array(df[y_col])

    # Calculate density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x_sorted, y_sorted, z_sorted = x[idx], y[idx], z[idx]

    # Calculate correlations
    spearman_corr = stats.spearmanr(x, y)[0]
    pearson_corr = stats.pearsonr(x, y)[0]
    print(f'Spearman correlation: {spearman_corr}')
    print(f'Pearson correlation: {pearson_corr}')

    # Create a new figure for each condition-column pair
    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)

    # Create scatter plot with density coloring
    ax.scatter(x_sorted, y_sorted, c=z_sorted, s=0.3)

    # Calculate and plot regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    x_line = np.array([min(x), max(x)])
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, linestyle=':', color='gray', linewidth=1)

    # Set title and labels
    ax.set_title(title, fontsize=titlesize)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ylim)

    # Set tick parameters
    ax.tick_params(axis='x', labelsize=16, width=1.5)
    ax.tick_params(axis='y', labelsize=16, width=1.5)

    # Add correlation text in the top-left corner
    ax.text(0.05, 0.95, f'R = {spearman_corr:.2f}\nr = {pearson_corr:.2f}',
            transform=ax.transAxes, fontsize=16, verticalalignment='top')

    # Remove top and right spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.tight_layout()
    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.show()


def correlation_heatmap(df, cols, title, savepath=None, x_labels=None, y_labels=None, mask_upper=True, vmin=0.5):
    """Create a correlation heatmap for the given columns."""
    # Select only the relevant columns and drop rows with NaN values
    dataframe = df[cols].copy().dropna()

    if len(dataframe) == 0:
        raise ValueError("No data to plot after dropping NaNs.")

    # Calculate the correlation matrix
    corr = dataframe.corr('pearson')

    # Create the mask for the upper triangle if requested
    if mask_upper:
        mask = np.triu(np.ones_like(corr, dtype=bool))
        np.fill_diagonal(mask, False)
    else:
        mask = None

    # Create the figure
    fig_dims = (9, 4)
    fig, ax = plt.subplots(figsize=fig_dims)

    ax.set_title(f'{title}\nCorrelations at {len(dataframe)} positions (Pearson)', fontsize=24)

    # Create the heatmap
    cmap = sn.color_palette("Blues", as_cmap=True)
    with sn.plotting_context(font_scale=1.8):
        hm = sn.heatmap(
            corr,
            annot=True,
            fmt='.2g',
            vmin=0.5,
            vmax=1,
            linewidths=1.5,
            linecolor='white',
            cbar=False,
            cmap=cmap,
            mask=mask,
            annot_kws={"size": 21},
            xticklabels=x_labels,
            yticklabels=y_labels,
            ax=ax
        )

        hm.tick_params(axis='x', labelsize=18)
        hm.tick_params(axis='y', labelsize=18)

    # Adjust appearance
    ax.tick_params(left=False, bottom=False)
    ax.set_facecolor("white")
    plt.yticks(rotation=0)

    # Save the figure if a path is provided
    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')

    plt.show()


def focused_correlation_heatmap(corr_matrix, title, x_labels=None, y_labels=None,
                                savepath=None, vmin=0,
                                fig_width=None, fig_height=4):
    """
    Create a heatmap for a focused portion of a correlation matrix.

    Parameters:
    -----------
    corr_matrix : pandas.DataFrame
        Pre-calculated correlation matrix to display
    title : str
        Title for the plot
    x_labels, y_labels : list, optional
        Custom labels for x and y axes. If None, uses matrix indices.
    savepath : str, optional
        Path to save the figure
    vmin : float
        Minimum value for the colormap
    fig_width, fig_height : float
        Width and height of the figure
    """
    # Set default figure width based on number of columns if not provided
    if fig_width is None:
        fig_width = max(9, corr_matrix.shape[1] * 1.5)

    # Use matrix indices as labels if not specified
    if x_labels is None:
        x_labels = corr_matrix.columns
    if y_labels is None:
        y_labels = corr_matrix.index

    # Create the figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Set title with larger font size
    ax.set_title(title, fontsize=24)

    # Create the heatmap using context manager
    cmap = sn.color_palette("Blues", as_cmap=True)
    with sn.plotting_context(font_scale=1.8):
        hm = sn.heatmap(
            corr_matrix,
            annot=True,
            fmt='.2g',
            vmin=vmin,
            vmax=1,
            linewidths=1.5,
            linecolor='white',
            cbar=True,  # Include color bar for focused heatmaps
            cmap=cmap,
            annot_kws={"size": 21},
            xticklabels=x_labels,
            yticklabels=y_labels,
            ax=ax
        )

        # Increase axis labels font size
        hm.tick_params(axis='x', labelsize=18, rotation=90)  # Angled for typically longer x labels
        hm.tick_params(axis='y', labelsize=18, rotation=0)

    # Adjust appearance
    ax.tick_params(left=False, bottom=False)
    ax.set_facecolor("white")

    # Save the figure if a path is provided
    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')

    plt.tight_layout()
    plt.show()


def create_histogram(df, column, title, xlabel='Editing efficiency (%)',
                     ylabel='Number of barcodes', range_limits=(0, 100), bins=50, savepath=None):
    """
    Create a histogram for editing efficiency data with mean line.

    Parameters:
    -----------
    df : pandas DataFrame
        DataFrame containing the data to plot
    column : str
        Column name in df to create histogram for
    title : str
        Title for the plot
    xlabel : str, optional
        X-axis label
    ylabel : str, optional
        Y-axis label
    range_limits : tuple, optional
        Range for the x-axis (min, max)
    bins : int, optional
        Number of histogram bins
    savepath : str, optional
        If provided, save the plot to this path

    Returns:
    --------
    None : Displays or saves the plot
    """
    # Reset to matplotlib defaults first
    plt.rcdefaults()

    # Then reapply your custom styling using the global variables
    plt.rcParams['axes.linewidth'] = axesthickness
    mpl.rcParams['pdf.fonttype'] = pdffonttpye
    mpl.rcParams['ps.fonttype'] = psfonttype
    mpl.rcParams['font.family'] = fontfamily

    # Print data statistics
    print(f"Full length of df: {len(df)}")
    data = df[column].dropna()
    print(f"Plotted number after dropping NaNs: {len(data)}")

    if len(data) == 0:
        raise ValueError(f"No data to plot for column {column} after dropping NaNs.")

    # Calculate statistics
    mean_val = data.mean()
    median_val = data.median()
    std_val = data.std()

    print(f"Mean: {mean_val:.2f}%")
    print(f"Median: {median_val:.2f}%")
    print(f"Std Dev: {std_val:.2f}%")

    # Create the plot
    fig, ax = plt.subplots(figsize=(5, 3), dpi=150)

    # Plot histogram
    ax.hist(data, bins=bins, range=range_limits,
            color='lightgrey', edgecolor='black', linewidth=1)

    # Add a vertical dotted line at the mean
    ax.axvline(mean_val, color='black', linestyle='dotted')

    # Add text with mean value
    ax.text(0.65, 0.85, f'Mean: {mean_val:.1f}%',
            transform=ax.transAxes, fontsize=12)

    # Set title and labels
    ax.set_title(title, fontsize=titlesize)
    ax.set_xlabel(xlabel, fontsize=lsize)
    ax.set_ylabel(ylabel, fontsize=lsize)
    ax.set_xlim(range_limits)

    # Set tick parameters
    ax.tick_params(axis='x', labelsize=lsize, width=1.5)
    ax.tick_params(axis='y', labelsize=lsize, width=1.5)

    # Remove top and right spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.tight_layout()
    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.show()
