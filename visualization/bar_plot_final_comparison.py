import matplotlib.pyplot as plt
from utils.led import LST1_01G01_UV01_00, LST1_01F06_GRN1_00
from math import nan
from ICIW_Plots import make_square_subplots, cm2inch
import ICIW_Plots.colors as ICIWcolors
import os
from utils.process_radiometry_data import process_radiometry_file, calculate_correction_factor
from utils.process_radiometry_indirect import process_radiometry_data as process_radiometry_indirect_data

def create_bar_plot_final_comparison():
    r"""
    Generates and displays a comparative bar plot of mean irradiance and homogeneity index for two LED types (365 nm and 530 nm) under three different module configurations: indirect irradiation, manual optimization, and model-based optimization. The function reads radiometry data from CSV files, processes and normalizes the results, and visualizes the comparison in a two-panel figure. The resulting plot is saved as an SVG file and shown interactively.

    **High-Level Description**
        - Reads and processes radiometry data for two LED types and three configurations per type.
        - Calculates mean irradiance (W m\ :sup:`-2`) and homogeneity index (unitless) for each case.
        - Plots grouped bar charts for both metrics, annotates values, and saves the figure.

    :param None: This function takes no arguments.

    :returns: None
        This function does not return any value. It produces a side effect by saving and displaying a plot.

    :raises FileNotFoundError: If any required CSV data file is missing.
    :raises ValueError: If the data in the CSV files is malformed or cannot be processed.

    :side effects:
        - Saves a figure to ``figures/bar_plot_final_comparison.svg`` in the current working directory.
        - Displays a matplotlib window with the generated plot.

    **Notes**
        - Assumes the existence of radiometry data files in the ``data/radiometry/{led.name}/`` directory structure.
        - Requires custom matplotlib stylesheets: ``ICIWstyle`` and ``visualization/publication_style.mplstyle``.
        - Uses external utility functions for data processing and correction factor calculation.
        - The number of LEDs per configuration is hardcoded (8 for manual, 16 for model-based).

    **Usage Example**

    .. code-block:: python

        from visualization.bar_plot_final_comparison import create_bar_plot_final_comparison
        create_bar_plot_final_comparison()

    This will generate and display the comparison plot, saving it as an SVG file.
    """
    # Set plot style
    plt.style.use("ICIWstyle")
    plt.style.use("visualization/publication_style.mplstyle")

    # read data for final comparison bar plot
    bar_data = {}
    for led in [LST1_01G01_UV01_00, LST1_01F06_GRN1_00]:
        # Indirect Irradiation
        filepath_indirect = f"data/radiometry/{led.name}/indirect_irradiation.csv"
        irr_indirect, I_mean_indirect, H_indirect = process_radiometry_indirect_data(filepath_indirect)
        
        # Manually Optimized Configuration
        N_LEDs = 8 # number of in the module
        expected_value = N_LEDs * led.total_power  # W m⁻²
        filepath_manual = f"data/radiometry/{led.name}/best_manual_reflector.csv"
        filepath_manual_ref = f"data/radiometry/{led.name}/best_manual_no_reflector.csv"
        correction_factor = calculate_correction_factor(filepath_manual_ref, expected_value)
        irr_manual, I_mean_manual, H_manual = process_radiometry_file(filepath_manual, correction_factor)
        
        # Optimized Model-Based Configuration
        N_LEDs = 16 # number of in the module
        expected_value = N_LEDs * led.total_power  # W m⁻²
        if led == LST1_01G01_UV01_00:
            filepath_model = f"data/radiometry/{led.name}/UV3_reflector.csv"
            filepath_model_ref = f"data/radiometry/{led.name}/UV3_reflector.csv"
        else:
            filepath_model = f"data/radiometry/{led.name}/GRN3_reflector.csv"
            filepath_model_ref = f"data/radiometry/{led.name}/GRN3_reflector.csv"
        correction_factor = calculate_correction_factor(filepath_model_ref, expected_value)
        irr_model, I_mean_model, H_model = process_radiometry_file(filepath_model, correction_factor)
        
        bar_data[led.name] = {
            "indirect": (I_mean_indirect, H_indirect),
            "manual": (I_mean_manual, H_manual),
            "model": (I_mean_model, H_model),
        }


    # Combine all data into groups for each wavelength and case
    labels = ['365 nm LEDs', '530 nm LEDs']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(19 * cm2inch, 10 * cm2inch))

    x = range(len(labels))
    width = 0.25

    # Bar positions for each group
    positions = [
        [p - width for p in x],  # Geometric Model
        [p for p in x],  # Raytracing
        [p + width for p in x],  # Radiometry
    ]

    # Mean Irradiance Bar Plot
    ax1.bar(positions[0], [bar_data["LST1-01G01-UV01-00"]["indirect"][0], bar_data["LST1-01F06-GRN1-00"]["indirect"][0]], width=width, label='Indirect Irradiation', color="C0", alpha=1)
    ax1.bar(positions[1], [bar_data["LST1-01G01-UV01-00"]["manual"][0], bar_data["LST1-01F06-GRN1-00"]["manual"][0]], width=width, label='Manual Optimization', color="C0", alpha=.67)
    ax1.bar(positions[2], [bar_data["LST1-01G01-UV01-00"]["model"][0], bar_data["LST1-01F06-GRN1-00"]["model"][0]], width=width, label='Model-Based Optimization', color="C0", alpha=.33)
    ax1.set_ylabel('Mean Irradiance / W m$^{-2}$')
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_ylim(0, 130)

    # Homogeneity Index Bar Plot
    ax2.bar(positions[0], [bar_data["LST1-01G01-UV01-00"]["indirect"][1], bar_data["LST1-01F06-GRN1-00"]["indirect"][1]], width=width, label='Indirect Irradiation', color="C0", alpha=1)
    ax2.bar(positions[1], [bar_data["LST1-01G01-UV01-00"]["manual"][1], bar_data["LST1-01F06-GRN1-00"]["manual"][1]], width=width, label='Manual Optimization', color="C0", alpha=0.67)
    ax2.bar(positions[2], [bar_data["LST1-01G01-UV01-00"]["model"][1], bar_data["LST1-01F06-GRN1-00"]["model"][1]], width=width, label='Model-Based Optimization', color="C0", alpha=0.33)
    ax2.set_ylabel('Homogeneity / 1')
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_ylim(0, 1.)

    # show values on top of bars
    for rects in ax1.patches:
        height = rects.get_height()
        ax1.annotate(f'{height:.1f}',
                    xy=(rects.get_x() + rects.get_width() / 2, height),
                    xytext=(0, 2),  # 2 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
        
    for rects in ax2.patches:
        height = rects.get_height()
        ax2.annotate(f'{height:.2f}',
                    xy=(rects.get_x() + rects.get_width() / 2, height),
                    xytext=(0, 2),  # 2 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.08))

    plt.tight_layout()
    fig.savefig(os.path.join("figures", "bar_plot_final_comparison.svg"), dpi=300, bbox_inches='tight')
    plt.show()
