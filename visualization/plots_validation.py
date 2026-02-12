import matplotlib.pyplot as plt
from geometric_model import GeometricModel
from optimizer import Optimizer
from utils.led import LST1_01G01_UV01_00, LST1_01F06_GRN1_00
from ICIW_Plots import make_square_subplots, cm2inch
import os
from utils.process_raytracing_data import process_raytracing_file, calc_metrics, plot_irradiance_raytracing
from utils.process_radiometry_data import process_radiometry_file, calculate_correction_factor, plot_irradiance_radiometry
from utils.grid import Grid
import numpy as np

def create_plots_validation():
    r"""
    Generates and displays a comparative plots of irradiance, mean irradiance, and homogeneity index for two LED types (365 nm and 530 nm) for the best manual configuration obtained by three different methods: geometric model, raytracing, and radiometry. The plot is saved as an SVG file in the 'figures' directory.

    **High-Level Description**
        - Calculates the geometric model results.
        - Processes raytracing data for the best manual configuration.
        - Processes radiometry data for the best manual configuration, applying a correction factor based on reference measurements.
        - Calculates mean irradiance (W m\ :sup:`-2`) and homogeneity index (unitless) for each case.
        - Create three irradiance plots for geometric model, raytracing, and radiometry for each LED.
        - Plots grouped bar charts for both metrics, annotates values, and saves the figure.

    :param None: This function takes no arguments.

    :returns: None
        This function does not return any value. It produces a side effect by saving and displaying a plot.

    :raises FileNotFoundError: If any required CSV data file is missing.
    :raises ValueError: If the data in the CSV files is malformed or cannot be processed.

    :side effects:
        - Saves figures to ``figures/`` in the current working directory called ``validation_best_manual_LST1_01G01_UV01_00``, ``validation_best_manual_LST1_01F06_GRN1_00``, and `bar_plots_validation.svg`.
        - Displays a matplotlib window with the generated plot.

    **Notes**
        - Assumes the existence of rdata files in the ``results/raytracing/{led.name}/`` and ``results/radiometry/{led.name}/`` directory.
        - Requires custom matplotlib stylesheets: ``ICIWstyle`` and ``visualization/publication_style.mplstyle``.
        - Uses external utility functions for data processing and correction factor calculation.
        - The number of LEDs per configuration is hardcoded (8 for manual, 16 for model-based).

    **Usage Example**

    .. code-block:: python

        from visualization.bar_plots_validation import create_bar_plot_validation
        create_bar_plot_validation()

    This will generate and display the comparison plot, saving it as an SVG file.
    """
    # Set plot style
    plt.style.use("ICIWstyle")
    plt.style.use("visualization/publication_style.mplstyle")

    # read data for final validation bar plot
    data = {}
    for led in [LST1_01G01_UV01_00, LST1_01F06_GRN1_00]:
        
        # Geometric model
        grid = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
        model = GeometricModel(grid, led, resolution_xy=0.5)
        optimizer = Optimizer(irradiance_bound=150.0)
        best_manual_config = [(-7.5,0), (7.5,0), (-10, -10), (10, 10), (-10, 10), (10, -10), (0,12.5), (0,-12.5)]
        irr_geom = model.simulate(best_manual_config, 13.0)
        E_mean_geom = np.mean(irr_geom)
        std_geom = np.std(irr_geom)
        H_geom = 1 - std_geom / E_mean_geom

        # Raytracing
        path = f"results/raytracing/{led.name}/"
        filename = "best_manual_no_reflector.txt"
        irr_ray = process_raytracing_file(path=path, filename=filename)
        E_mean_ray, H_ray = calc_metrics(irr_ray["Irradiance"])

        # Radiometry
        N_LEDs = 8 # number of in the module
        expected_value = N_LEDs * led.total_power  # W m⁻²
        filepath_rad = f"results/radiometry/{led.name}/best_manual_no_reflector.csv"
        filepath_rad_ref = f"results/radiometry/{led.name}/best_manual_no_reflector.csv"
        correction_factor = calculate_correction_factor(filepath_rad_ref, expected_value)
        irr_rad, E_mean_rad, H_rad = process_radiometry_file(filepath_rad, correction_factor)
        
        # Store results in data dictionary for bar plot
        data[led.name] = {
            "geometric": (E_mean_geom, H_geom),
            "raytracing": (E_mean_ray, H_ray),
            "radiometry": (E_mean_rad, H_rad),
        }

        # Create 3 subplots for irradiance distribution for geometric model, raytracing, and radiometry for each LED
        fig, axs = plt.subplots(1, 3, figsize=(19 * cm2inch, 10 * cm2inch))
        if led == LST1_01G01_UV01_00:
            v_max = 70
        else:
            v_max = 20
        y_title = - 0.6
        shrink_cbar = 0.32
        extent = (
            irr_ray['x'].min(), irr_ray['x'].max(),
            irr_ray['y'].min(), irr_ray['y'].max()
        )
        model.plot_irradiance(axs[0], irr_geom, v_max=v_max, shrink_cbar=shrink_cbar, extent=extent)
        plot_irradiance_raytracing(axs[1], irr_ray, v_max=v_max, y_title=y_title, shrink_cbar=shrink_cbar, extent=extent)
        plot_irradiance_radiometry(axs[2], irr_rad, v_max=v_max, y_title=y_title, shrink_cbar=shrink_cbar, extent=extent)
        plt.tight_layout()
        # increase space between subplots
        plt.subplots_adjust(wspace=0.8)
        fig.savefig(os.path.join("figures", f"validation_best_manual_{led.name}.svg"), dpi=300, bbox_inches='tight')
        plt.show()



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
    ax1.bar(positions[0], [data["LST1-01G01-UV01-00"]["geometric"][0], data["LST1-01F06-GRN1-00"]["geometric"][0]], width=width, label='Geometric Model', color="C0", alpha=1)
    ax1.bar(positions[1], [data["LST1-01G01-UV01-00"]["raytracing"][0], data["LST1-01F06-GRN1-00"]["raytracing"][0]], width=width, label='Raytracing', color="C0", alpha=.67)
    ax1.bar(positions[2], [data["LST1-01G01-UV01-00"]["radiometry"][0], data["LST1-01F06-GRN1-00"]["radiometry"][0]], width=width, label='Radiometry', color="C0", alpha=.33)
    ax1.set_ylabel('Mean Irradiance / W m$^{-2}$')
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_ylim(0, 50)

    # Homogeneity Index Bar Plot
    ax2.bar(positions[0], [data["LST1-01G01-UV01-00"]["geometric"][1], data["LST1-01F06-GRN1-00"]["geometric"][1]], width=width, label='Geometric Model', color="C0", alpha=1)
    ax2.bar(positions[1], [data["LST1-01G01-UV01-00"]["raytracing"][1], data["LST1-01F06-GRN1-00"]["raytracing"][1]], width=width, label='Raytracing', color="C0", alpha=0.67)
    ax2.bar(positions[2], [data["LST1-01G01-UV01-00"]["radiometry"][1], data["LST1-01F06-GRN1-00"]["radiometry"][1]], width=width, label='Radiometry', color="C0", alpha=0.33)
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
    fig.savefig(os.path.join("figures", "bar_plot_validation.svg"), dpi=300, bbox_inches='tight')
    plt.show()