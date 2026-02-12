import pandas as pd
import matplotlib.pyplot as plt
from ICIW_Plots import cm2inch, make_square_subplots
import numpy as np

def create_ray_independence_plots_SI():
    """
    Generates and saves a set of 2x2 subplots visualizing the dependence of raytracing simulation results on the number of rays used. The function compares the standard deviation of irradiance and computation time for two wavelengths (365 nm and 530 nm), each with and without reflectors, as a function of the number of rays. Data is read from an Excel file, and the resulting figure is saved as an SVG file.

    The function is intended to help assess the statistical independence and convergence of raytracing results with increasing ray count, and to visualize the computational cost associated with higher ray numbers.

    :Parameters:
        None

    :Returns:
        None

    :Raises:
        FileNotFoundError: If the Excel data file is not found at the specified path.
        KeyError: If the expected sheet or data ranges are missing in the Excel file.
        ValueError: If the data cannot be converted to float as expected.

    :Side Effects:
        - Reads data from 'data/raytracing/ray_independence_test.xlsx'.
        - Saves a figure to 'figures/ray_independence_plots.svg'.
        - Displays the generated plots using matplotlib's interactive window.

    :Notes:
        - Assumes the Excel file contains a sheet named 'Planilha1' with data in specific cell ranges.
        - Requires the custom plotting style 'ICIWstyle' and 'visualization/publication_style.mplstyle' to be available.
        - Uses custom utility functions `cm2inch` and `make_square_subplots` from the `ICIW_Plots` module.
        - The function does not accept any parameters and is intended to be run as a script or called from another script.

    :Example:

        >>> create_ray_independence_plots_SI()
        # Generates and displays the ray independence plots, and saves them as an SVG file.

    """
    
    # Set plot style
    plt.style.use("ICIWstyle")
    plt.style.use("visualization/publication_style.mplstyle")

    # read results/raytracing/ray_indepence_test.xlsx
    path = 'results\\raytracing\\ray_independence_test.xlsx'
    df = pd.read_excel(path, sheet_name=None)

    data = {}

    # extract data for 365 nm with and without reflectors
    data['a) 365 nm with reflectors'] = df['Planilha1'].iloc[2:9, 0:3].to_numpy().astype(float)
    data['b) 365 nm without reflectors'] = df['Planilha1'].iloc[2:9, 3:6].to_numpy().astype(float)

    # extract data for 530 nm with and without reflectors
    data['c) 530 nm with reflectors'] = df['Planilha1'].iloc[14:22, 0:3].to_numpy().astype(float)
    data['d) 530 nm without reflectors'] = df['Planilha1'].iloc[14:22, 3:6].to_numpy().astype(float)
    # create 2x2 subplots with data[label[1]] on left and data[label[2]] on right
    fig = plt.figure(figsize=(16*cm2inch, 18*cm2inch))
    axes = make_square_subplots(
        fig,
        ax_layout=(2, 2),
        ax_width=6*cm2inch,
        v_sep=2.8*cm2inch,
        h_sep=3.5*cm2inch,
        xlabel='number of rays / 1',
        ylabel='standard deviation / W m$^{-2}$',
    )

    labels = list(data.keys())
    legend_labels, legend_handles = [], []
    for i, ax in enumerate(axes.flat):
        x = data[labels[i]][:, 0]
        y = data[labels[i]][:, 2]
        ax.plot(x, y, marker='o', linestyle='-', label="standard deviation / W m$^{-2}$", color='C0')
        # plot vertical line at x=1.6e7
        ax.axvline(x=1.6e7, color='k', linestyle=':')
        # plot time on twin axis
        ax2 = ax.twinx()
        time = data[labels[i]][:, 1]
        ax2.plot(x, time, marker='d', linestyle='--', color='C3', label="time / s")
        ax2.set_ylabel('time / s')
        # ax2.tick_params(axis='y', labelcolor='C3')
        # set log scale for x axis
        ax.set_xscale('log')
        ax.set_title(labels[i], y=-0.37, fontsize=10)

        if ax == axes[0, 0]:
            legend_labels = [ax.get_legend_handles_labels()[1][0], ax2.get_legend_handles_labels()[1][0]]
            legend_handles = [ax.lines[0], ax2.lines[0]]

    fig.legend(legend_handles, legend_labels, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.))
    plt.savefig('figures/ray_independence_plots.svg', dpi=300, bbox_inches='tight')
    plt.show()


