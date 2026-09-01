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
        - Assumes the existence of radiometry data files in the ``results/radiometry/{led.name}/`` directory structure.
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

    # number of LEDs of the three compared irradiation concepts. The concepts use different numbers
    # of LEDs, which is why the mean irradiance is additionally reported per LED below.
    n_leds = {"indirect": 6, "manual": 8, "model": 16}

    # read data for final comparison bar plot
    bar_data = {}
    for led in [LST1_01G01_UV01_00, LST1_01F06_GRN1_00]:
        # Indirect Irradiation
        filepath_indirect = f"results/radiometry/{led.name}/indirect_irradiation.csv"
        irr_indirect, I_mean_indirect, H_indirect = process_radiometry_indirect_data(filepath_indirect)

        # Manually Optimized Configuration
        N_LEDs = n_leds["manual"] # number of LEDs in the module
        expected_value = N_LEDs * led.total_power  # W m⁻²
        filepath_manual = f"results/radiometry/{led.name}/best_manual_reflector.csv"
        filepath_manual_ref = f"results/radiometry/{led.name}/best_manual_no_reflector.csv"
        correction_factor = calculate_correction_factor(filepath_manual_ref, expected_value)
        irr_manual, I_mean_manual, H_manual = process_radiometry_file(filepath_manual, correction_factor)
        
        # Optimized Model-Based Configuration
        N_LEDs = n_leds["model"] # number of LEDs in the module
        expected_value = N_LEDs * led.total_power  # W m⁻²
        if led == LST1_01G01_UV01_00:
            filepath_model = f"results/radiometry/{led.name}/UV3_reflector.csv"
            filepath_model_ref = f"results/radiometry/{led.name}/UV3_reflector.csv"
        else:
            filepath_model = f"results/radiometry/{led.name}/GRN3_reflector.csv"
            filepath_model_ref = f"results/radiometry/{led.name}/GRN3_reflector.csv"
        correction_factor = calculate_correction_factor(filepath_model_ref, expected_value)
        irr_model, I_mean_model, H_model = process_radiometry_file(filepath_model, correction_factor)
        
        bar_data[led.name] = {
            "indirect": (I_mean_indirect, H_indirect, I_mean_indirect / n_leds["indirect"]),
            "manual": (I_mean_manual, H_manual, I_mean_manual / n_leds["manual"]),
            "model": (I_mean_model, H_model, I_mean_model / n_leds["model"]),
        }


    # Combine all data into groups for each wavelength and case. Short tick labels are used because
    # the figure now holds three panels.
    labels = ['365 nm', '530 nm']
    concepts = [
        ("indirect", f"Indirect ({n_leds['indirect']} LEDs)", 1.0),
        ("manual", f"Manual Optimization ({n_leds['manual']} LEDs)", 0.67),
        ("model", f"Model-Based Optimization ({n_leds['model']} LEDs)", 0.33),
    ]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(19 * cm2inch, 8 * cm2inch))

    x = range(len(labels))
    # The bars are drawn edge to edge, so their width also sets the distance between the value
    # labels above them. 0.29 gives the labels enough room to not touch each other, which matters
    # in the homogeneity panel where all bars are of similar height.
    width = 0.29

    # Bar positions for each group
    positions = [
        [p - width for p in x],  # Indirect Irradiation
        [p for p in x],  # Manual Optimization
        [p + width for p in x],  # Model-Based Optimization
    ]

    def values(concept: str, index: int) -> list:
        """Collect one metric of one irradiation concept for both LED types."""
        return [bar_data["LST1-01G01-UV01-00"][concept][index],
                bar_data["LST1-01F06-GRN1-00"][concept][index]]

    def annotate_bars(ax, fmt: str):
        """Write the value of each bar above the bar."""
        for rects in ax.patches:
            height = rects.get_height()
            ax.annotate(format(height, fmt),
                        xy=(rects.get_x() + rects.get_width() / 2, height),
                        xytext=(0, 2),  # 2 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom',
                        fontsize=8)  # smaller than the axis labels so that the values fit

    # (a) Mean irradiance
    # (b) Mean irradiance per LED. The three concepts use 6, 8 and 16 LEDs operated at the same
    #     driving current, so normalizing by the number of LEDs allows a comparison that is
    #     independent of the number of installed light sources.
    # (c) Homogeneity
    for axis, index, ylabel, ylim, fmt in [
        (ax1, 0, 'Mean Irradiance / W m$^{-2}$', (0, 130), '.1f'),
        (ax2, 2, 'Mean Irradiance per LED / W m$^{-2}$', (0, 9), '.2f'),
        (ax3, 1, 'Homogeneity / 1', (0, 1.), '.2f'),
    ]:
        for position, (concept, label, alpha) in zip(positions, concepts):
            axis.bar(position, values(concept, index), width=width, label=label, color="C0",
                     alpha=alpha)
        axis.set_ylabel(ylabel)
        axis.set_xticks(x)
        axis.set_xticklabels(labels)
        axis.set_ylim(*ylim)
        annotate_bars(axis, fmt)

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.08))

    plt.tight_layout()
    fig.savefig(os.path.join("figures", "bar_plot_final_comparison.svg"), dpi=300, bbox_inches='tight')
    plt.show()

    # print the values including the normalization per LED
    for led_name, entries in bar_data.items():
        for concept, (mean, homogeneity, mean_per_led) in entries.items():
            print(f"{led_name} | {concept:<9s} | {n_leds[concept]:2d} LEDs | "
                  f"E_mean {mean:6.2f} W m-2 | E_mean/LED {mean_per_led:5.2f} W m-2 | "
                  f"H {homogeneity:.3f}")
