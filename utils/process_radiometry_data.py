import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.colors as mcolors
import os
from pathlib import Path
from utils.led import LST1_01F06_GRN1_00, LST1_01G01_UV01_00

def read_radiometry_file(filepath: str) -> pd.DataFrame:
    """
    Loads radiometry data from a CSV file and returns it as a pandas DataFrame.

    The function attempts to read the file using two different decimal and separator conventions to accommodate variations in CSV formatting. The first attempt assumes European-style CSVs (comma decimal, semicolon separator). If the resulting DataFrame is empty, it retries with standard CSV formatting (dot decimal, comma separator). Column and index names are converted to floats for consistency.

    :param filepath: Path to the CSV file containing radiometry data.
    :type filepath: str

    :returns: DataFrame containing irradiance data, with both index and columns as floats representing spatial positions (typically in cm).
    :rtype: pandas.DataFrame

    :raises FileNotFoundError: If the file does not exist at the given path.
    :raises pd.errors.ParserError: If the file cannot be parsed as a CSV.

    :side effects: None.

    :notes:
        - Assumes the first row and column of the CSV are headers and indices, respectively.
        - Handles both European and standard CSV formats.
        - All positions are assumed to be in centimeters unless otherwise specified.
    """

    # Read the CSV file
    irr = pd.read_csv(filepath, decimal=',', sep=';', header=0, index_col=0)
    # ensure that column names are numeric
    irr.columns = irr.columns.str.replace(',', '.').astype(float)

    # if dataframe empty read in differently
    if irr.empty:
        irr = pd.read_csv(filepath, decimal='.', sep=',', header=0, index_col=0)
        irr.columns = irr.columns.astype(float)
        irr.index = irr.index.astype(float)
    return irr

def process_radiometry_file(filepath: str, correction_factor: float = 1.0, height: float = 33, width: float = 34) -> tuple[pd.DataFrame, float, float]:
    """
    Processes radiometry data from a CSV file, applies a correction factor, centers the data, and computes summary statistics.

    The function reads the radiometry data, centers the detection plane around zero, applies a correction factor to the irradiance values, and calculates the mean irradiance and homogeneity of the distribution. The region of interest is defined by the specified height and width (in cm).

    :param filepath: Path to the CSV file containing radiometry data.
    :type filepath: str
    :param correction_factor: Multiplicative factor to correct the irradiance data (default: 1.0).
    :type correction_factor: float, optional
    :param height: Height of the detection plane in centimeters (default: 33).
    :type height: float, optional
    :param width: Width of the detection plane in centimeters (default: 34).
    :type width: float, optional

    :returns: Tuple containing:
        - Processed irradiance data as a DataFrame (index/columns: centered positions in cm).
        - Mean irradiance in W/m² (float).
        - Homogeneity of the irradiance distribution (float, 0 to 1).
    :rtype: tuple[pandas.DataFrame, float, float]

    :raises FileNotFoundError: If the file does not exist.
    :raises pd.errors.ParserError: If the file cannot be parsed.

    :side effects: None.

    :notes:
        - The homogeneity is defined as 1 - (std/mean) of the irradiance values.
        - The region is centered based on the mean of the column positions.
        - Assumes spatial units are centimeters.
    """

    # Read the CSV file
    df = read_radiometry_file(filepath)
    center = df.columns.to_numpy().mean()

    half_width = width / 2
    half_height = height / 2
    df = df.loc[(df.index >= center - half_height) & (df.index <= center + half_height), (df.columns >= center - half_width) & (df.columns <= center + half_width)]
    # reindex and recolumns to be centered around 0
    df.index = df.index - center
    df.columns = df.columns - center

    # Apply correction factor to calculate the corrected irradiance
    irr = df * correction_factor

    # calculate mean irradiance and homogeneity
    mean_irradiance = irr.to_numpy().mean()
    homogeneity = 1 - irr.to_numpy().std() / mean_irradiance

    return irr, mean_irradiance, homogeneity

def process_all_cand_files(folder, correction_factor=1.0, height=33, width=34):
    """
    Processes all radiometry candidate CSV files in a folder and summarizes their irradiance statistics.

    The function iterates over all CSV files in the specified folder whose filenames contain 'UV' or 'GRN', processes each file, and collects the irradiance DataFrame, mean irradiance, and homogeneity in a results dictionary.

    :param folder: Path to the folder containing radiometry candidate CSV files.
    :type folder: str
    :param correction_factor: Multiplicative factor to correct the irradiance data (default: 1.0).
    :type correction_factor: float, optional
    :param height: Height of the detection plane in centimeters (default: 33).
    :type height: float, optional
    :param width: Width of the detection plane in centimeters (default: 34).
    :type width: float, optional

    :returns: Dictionary mapping filenames to a dictionary with keys 'irradiance' (DataFrame), 'mean_irradiance' (float), and 'homogeneity' (float).
    :rtype: dict[str, dict[str, Any]]

    :raises FileNotFoundError: If the folder does not exist.
    :raises pd.errors.ParserError: If a file cannot be parsed.

    :side effects: Reads all matching CSV files in the specified folder.

    :notes:
        - Only files ending with '.csv' and containing 'UV' or 'GRN' in the filename are processed.
        - The irradiance DataFrame is centered and corrected as in :func:`process_radiometry_file`.
    """
    results = {}
    for file in os.listdir(folder):
        if file.endswith('.csv') and ('UV' in file or 'GRN' in file):
            filepath = os.path.join(folder, file)
            irr, mean_irradiance, homogeneity = process_radiometry_file(filepath, correction_factor, height, width)
            results[file] = {
                'irradiance': irr,
                'mean_irradiance': mean_irradiance,
                'homogeneity': homogeneity
            }
    return results

def calculate_correction_factor(reference_filepath: str, expected_value: float) -> float:
    """
    Computes a correction factor for radiometry data using a reference measurement and an expected total irradiance value.

    The function integrates the measured irradiance over the detection area (using pixel size inferred from the DataFrame's index/columns), compares it to the expected value, and returns the scaling factor needed to match the expected irradiance.

    :param reference_filepath: Path to the reference CSV file containing measured irradiance data.
    :type reference_filepath: str
    :param expected_value: Expected total irradiance value in W/m².
    :type expected_value: float

    :returns: Correction factor to apply to radiometry data (float).
    :rtype: float

    :raises FileNotFoundError: If the reference file does not exist.
    :raises pd.errors.ParserError: If the file cannot be parsed.
    :raises IndexError: If the DataFrame does not contain enough data to infer pixel size.

    :side effects: None.

    :notes:
        - The units of the original data (e.g., umol m⁻² s⁻¹ or W m⁻²) do not affect the correction factor calculation.
        - Assumes uniform pixel size and that index/columns are in centimeters.
    """

    df = read_radiometry_file(reference_filepath)
    dx = (df.columns.to_numpy()[1] - df.columns.to_numpy()[0]) * 0.01 # width of one pixel in meters
    dy = (df.index.to_numpy()[1] - df.index.to_numpy()[0]) * 0.01 # height of one pixel in meters
    measured_value = df.to_numpy().sum() * (dx * dy)  # calculate total irradiance
    correction_factor = expected_value / measured_value # scale value to expected value. Note: If the original data is in umol m⁻² s⁻¹ or W m⁻², does not matter for the correction factor calculation, as the conversion factor between both units cancels out.
    return correction_factor

def plot_irradiance_radiometry(ax: Axes, irr: pd.DataFrame, title: str = "", y_title: float = -0.5, v_max: float = -1., extent: tuple = (-2, -1, -2, -1), shrink_cbar: float = 1.0) -> Axes:
        """
        Visualizes the irradiance distribution as a color-mapped image with a colorbar and labeled axes.

        The function displays the irradiance DataFrame as an image on the provided matplotlib Axes, with optional title, color scaling, and extent. The colorbar is labeled with irradiance units. The plot is set to have equal aspect ratio and axis labels in centimeters.

        :param ax: Matplotlib Axes object on which to plot the irradiance distribution.
        :type ax: matplotlib.axes.Axes
        :param irr: DataFrame containing irradiance data to plot (index/columns: positions in cm).
        :type irr: pandas.DataFrame
        :param title: Title of the plot (default: '').
        :type title: str, optional
        :param y_title: Y position of the title (default: -0.5).
        :type y_title: float, optional
        :param v_max: Maximum value for the color scale. If -1, auto-scales to nearest multiple of 5 above max value (default: -1).
        :type v_max: float, optional
        :param extent: Plot extent as (xmin, xmax, ymin, ymax). If default, uses DataFrame min/max (default: (-2, -1, -2, -1)).
        :type extent: tuple, optional
        :param shrink_cbar: Factor to shrink the colorbar (default: 1.0).
        :type shrink_cbar: float, optional

        :returns: The matplotlib Axes object with the plot.
        :rtype: matplotlib.axes.Axes

        :raises ValueError: If the DataFrame is empty or extent is invalid.

        :side effects: Modifies the provided Axes object and creates a colorbar.

        :notes:
            - The color scale is set to start at 0 and end at v_max.
            - The function auto-determines contour levels based on v_max.
            - The colorbar is labeled with '$E$ / W m⁻²'.
        """

        # Set color bar limits
        if v_max != -1.:
            vmax = v_max
        else:
            vmax = (int(irr.to_numpy().max() // 5) + 1) * 5 # round up to next multiple of 5
        norm = mcolors.Normalize(vmin=0, vmax=vmax, clip=False)

        # define contour levels based on vmax
        if vmax <= 30:
            step = 5
        if vmax <= 50:
            step = 10
        elif vmax <= 100:
            step = 20
        else:
            step = 50 
        levels = np.arange(0, vmax + step, step)

        # Create colormesh plot
        # Ensure extent values are floats
        if extent == (-2, -1, -2, -1):  
            x_min = float(irr.columns.min())
            x_max = float(irr.columns.max())
            y_min = float(irr.index.min())
            y_max = float(irr.index.max())
            extent = (x_min, x_max, y_min, y_max)
        tc = ax.imshow(
            irr,
            extent=extent,
            origin='lower',
            cmap='viridis',
            norm=norm,
            aspect='equal'
        )

        # set subtitle
        if title:
            ax.set_title(title, y=y_title, fontsize=10)

        # Create colorbar with reduced size to match plot
        cbar = plt.colorbar(tc, ax=ax, ticks=levels, shrink=shrink_cbar)
        cbar.set_label('$E$ / W m⁻²')

        ax.set_xlabel('$X$ / cm')
        ax.set_ylabel('$Y$ / cm')
        ax.set_aspect('equal', 'box')

        # Set x and y limits
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])

        return ax