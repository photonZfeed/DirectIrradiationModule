import matplotlib.pyplot as plt
from math import nan
from ICIW_Plots import make_rect_subplots, cm2inch
import os
from utils.process_radiometry_data import calculate_correction_factor, process_all_cand_files
from utils.process_raytracing_data import process_raytracing_file, calc_metrics
from systematic_sampler import Sampler
from utils.grid import Grid
from geometric_model import GeometricModel
from utils.led import LST1_01F06_GRN1_00, LST1_01G01_UV01_00
import numpy as np


def create_bar_plots(keys: list, geometric: dict, raytracing: dict, radiometric: dict, export_path: str = "", filename: str = "bar_plots_selected_configurations.svg", export_type: str = "svg"):
    """
    Generate comparative grouped bar plots for mean irradiance and homogeneity across multiple evaluation methods.

    This function visualizes and compares the results of three different evaluation methods—geometric modeling, raytracing simulation, and radiometric measurement—by creating grouped bar plots for each candidate configuration. For each method, the mean irradiance (``E_mean``) and homogeneity (``H``) are plotted side by side, allowing for direct visual comparison. The function supports both interactive display and export to file.

    Parameters
    ----------
    keys : list of str
        List of configuration names (labels) to be plotted. Each entry corresponds to a candidate configuration, e.g., 'UV1', 'UV1 - Reflector'.
    geometric : dict
        Dictionary mapping configuration names to geometric model results. Each value should be a dict with keys:
        - 'E_mean' (float): Mean irradiance (W/m²).
        - 'H' (float): Homogeneity (dimensionless, 0–1).
    raytracing : dict
        Dictionary mapping configuration names to raytracing simulation results. Each value should be a dict with keys:
        - 'E_mean' (float): Mean irradiance (W/m²).
        - 'H' (float): Homogeneity (dimensionless, 0–1).
    radiometric : dict
        Dictionary mapping configuration names to radiometric measurement results. Each value should be a dict with keys:
        - 'E_mean' (float): Mean irradiance (W/m²).
        - 'H' (float): Homogeneity (dimensionless, 0–1).
    export_path : str, optional
        Directory path where the plot will be saved. If empty (default), the plot is shown interactively instead of being saved.
    filename : str, optional
        Name of the output file (default: 'bar_plots_selected_configurations.svg'). Only used if ``export_path`` is provided.
    export_type : str, optional
        File format for export (default: 'svg'). Passed to ``matplotlib.pyplot.savefig``.

    Returns
    -------
    None
        This function does not return any value. It either displays the plot interactively or saves it to disk.

    Raises
    ------
    KeyError
        If expected keys are missing in the input dictionaries for any configuration.
    Exception
        Propagates exceptions from matplotlib or file I/O if saving fails.

    Side Effects
    ------------
    - Displays a matplotlib figure interactively (if ``export_path`` is not provided).
    - Saves a figure to disk (if ``export_path`` is provided).

    Notes
    -----
    - Assumes that all input dictionaries use the same set of configuration keys as provided in ``keys``.
    - Missing data for a method/configuration is represented as NaN and will appear as empty bars.
    - Uses a custom subplot layout and styling for publication-quality figures.
    - The x-axis labels are cleaned to remove ' - Reflector' for clarity.
    - The function expects ``make_rect_subplots`` and ``cm2inch`` utilities to be available in the environment.

    Example
    -------
    .. code-block:: python

        keys = ['UV1', 'UV2', 'UV1 - Reflector']
        geometric = {'UV1': {'E_mean': 20.1, 'H': 0.95}, ...}
        raytracing = {'UV1': {'E_mean': 19.8, 'H': 0.93}, ...}
        radiometric = {'UV1': {'E_mean': 18.5, 'H': 0.91}, ...}
        create_bar_plots(keys, geometric, raytracing, radiometric, export_path='./figures', filename='comparison', export_type='png')

    """
    types = ["Geometric Model", "Raytracing", "Radiometry"]

    # Prepare data
    E_means = {t: [] for t in types}
    H_values = {t: [] for t in types}
    labels = []

    for key in keys:
        labels.append(key)
        # Geometric
        if key in geometric:
            E_means["Geometric Model"].append(geometric[key]["E_mean"])
            H_values["Geometric Model"].append(geometric[key]["H"])
        else:
            E_means["Geometric Model"].append(nan)
            H_values["Geometric Model"].append(nan)
        # Raytracing
        if key in raytracing:
            E_means["Raytracing"].append(raytracing[key]["E_mean"])
            H_values["Raytracing"].append(raytracing[key]["H"])
        else:
            E_means["Raytracing"].append(nan)
            H_values["Raytracing"].append(nan)
        # Radiometric
        if key in radiometric:
            E_means["Radiometry"].append(radiometric[key]["E_mean"])
            H_values["Radiometry"].append(radiometric[key]["H"])
        else:
            E_means["Radiometry"].append(nan)
            H_values["Radiometry"].append(nan)

    x = range(len(labels))
    width = 0.25
    offsets = [-width, 0, width]

    fig = plt.figure(figsize=(18*cm2inch, 14*cm2inch))
    axs = make_rect_subplots(
        fig,
        ax_layout=(2, 1),
        ax_width=17*cm2inch,
        ax_height=5*cm2inch,
        v_sep=0.25,
        xlabel="",
        ylabel=[r"$E_\mathrm{mean}$ / W m⁻²", r"$H$ / 1"],
        sharelabel=True,
        sharex=True
    ).flatten()

    alpha_values = [1.0, 0.67, 0.33]
    for i, t in enumerate(types):
        axs[0].bar([xi + offsets[i] for xi in x], E_means[t], width=width, label=t, color="C0", alpha=alpha_values[i])
        axs[1].bar([xi + offsets[i] for xi in x], H_values[t], width=width, label=t, color="C0", alpha=alpha_values[i])

    labels = [label.replace(" - Reflector", "") for label in labels]
    axs[0].set_xticks(x)
    axs[0].set_xticklabels(labels, ha='right')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(labels, ha='right')
    # axs[0].set_ylim(0, 45)
    axs[1].set_ylim(0, 1.0)

    # legend above the plots
    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3)

    if export_path:
        plt.savefig(os.path.join(export_path, filename+"."+export_type), format=export_type, bbox_inches='tight')
    else:
        plt.show()

def create_bar_plot_screened_cands():

    """
    Generates comparative bar plots for screened candidate LED module configurations, visualizing results from three different evaluation methods: geometric modeling, raytracing simulation, and radiometric measurement.

    This function processes data for two types of LEDs (UV and green), automatically loading and mapping results from radiometric CSV files, geometric model simulations (from JSON configuration files), and raytracing text files. For each LED type, it creates grouped bar plots showing the mean irradiance (E_mean, in W/m²) and homogeneity (H, dimensionless, 0–1) for each candidate configuration, with and without reflectors. The plots are styled for publication and can be displayed interactively.

    Parameters
    ----------
    None
        This function does not accept any parameters. All paths and settings are hardcoded for the current project structure.

    Returns
    -------
    None
        The function does not return any value. It displays the generated plots using matplotlib's interactive window.

    Raises
    ------
    FileNotFoundError
        If required data files (CSV, JSON, or TXT) are missing in the expected directories.
    KeyError
        If expected keys are missing in the loaded data dictionaries.
    Exception
        Propagates exceptions from underlying data processing functions if data is malformed or incompatible.

    Side Effects
    ------------
    - Reads multiple files from disk (radiometry, raytracing, and configuration directories).
    - Displays matplotlib figures interactively (unless modified to export).

    Notes
    -----
    - Assumes the existence of specific directory structures and file naming conventions for radiometry, raytracing, and configuration data.
    - Uses a custom matplotlib style ('ICIWstyle' and 'visualization/publication_style.mplstyle').
    - The number of LEDs and their properties are hardcoded for each LED type.
    - The function is intended for exploratory analysis and figure generation in the context of LED module optimization.

    Usage Example
    -------------
    >>> create_bar_plot_screened_cands()
    # Displays comparative bar plots for all screened candidate configurations.

    """

    # Set plot style
    plt.style.use("ICIWstyle")
    plt.style.use("visualization/publication_style.mplstyle")

    # define height and width of detection plane
    height = 33  # cm
    width = 34  # cm

    for led_name in ["LST1-01G01-UV01-00", "LST1-01F06-GRN1-00"]:


        if led_name == "LST1-01G01-UV01-00":
            led = LST1_01G01_UV01_00  # 365 nm LED
            keys = ["UV1", "UV2", "UV3", "UV1 - Reflector", "UV2 - Reflector", "UV3 - Reflector"]
        else:
            led = LST1_01F06_GRN1_00 # 530 nm LED
            keys = ["GRN1", "GRN2", "GRN3", "GRN1 - Reflector", "GRN2 - Reflector", "GRN3 - Reflector"]

        # get radiometric data
        folder = f"./data/radiometry/{led.name}" # folder with all candidate files
        if led_name == "LST1-01G01-UV01-00":
            reference_file = f"./data/radiometry/{led.name}/UV3_reflector.csv" # reference file for correction factor
        else:
            reference_file = f"./data/radiometry/{led.name}/GRN3_reflector.csv" # reference file for correction factor
        N_led = 16 # number of LEDs in the module
        expected_value = N_led*led.total_power # expected irradiance value in W m^-2 according to LED datasheet
        correction_factor = calculate_correction_factor(reference_file, expected_value) # calculate correction factor based on reference file
        results_radiometry = process_all_cand_files(folder, correction_factor, height, width) # process all candidate files in the folder

        # map results to key names defined above, e.g UV1_no_reflector.csv -> UV1, UV1_reflector.csv -> UV1 - Reflector
        radiometric = {}
        for key in keys:
            if " - Reflector" in key:
                file_name = key.replace(" - Reflector", "_reflector.csv")
            else:
                file_name = key + "_no_reflector.csv"
            if file_name in results_radiometry:
                radiometric[key] = {}
                radiometric[key]["E_mean"] = results_radiometry[file_name]["mean_irradiance"]
                radiometric[key]["H"] = results_radiometry[file_name]["homogeneity"]

        # get geometric model data
        G = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
        Sampler1 = Sampler(G=G) 
        model = GeometricModel(G=G, led=led, resolution_xy=0.5)
        geometric = {}

        # read all json files from the folder with near optimal configurations
        for file in os.listdir(f"./sampled_configs/near_optimal_configs/{led.name}"):
            if file.endswith('.json'):
                config, height = Sampler1.read_configuration_from_json(os.path.join(f"./sampled_configs/near_optimal_configs/{led.name}", file))
                irr = model.simulate(config, height)
                geometric[file.replace(".json","")] = {
                    "E_mean": np.mean(irr),
                    "H": 1 - np.std(irr) / np.mean(irr),
                }

        # raytracing data
        raytracing = {}
        for key in keys:
            if " - Reflector" in key:
                file_name = key.replace(" - Reflector", "_reflector.txt")
            else:
                file_name = key + "_no_reflector.txt"
            file_path = os.path.join(f"./data/raytracing/{led.name}", file_name)
            if os.path.exists(file_path):
                raytracing_data = process_raytracing_file(file_name, os.path.join(f"./data/raytracing/{led.name}"))
                raytracing_result = calc_metrics(raytracing_data["Irradiance"])
                raytracing[key] = {
                    "E_mean": raytracing_result[0],
                    "H": raytracing_result[1],
                }
            
        create_bar_plots(keys, geometric, raytracing, radiometric) #, export_path="./figures", filename=f"bar_plots_screened_candidates_{led_name}", export_type="svg")