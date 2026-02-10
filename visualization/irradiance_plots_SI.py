import os
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from utils.process_raytracing_data import process_raytracing_file, calc_metrics, plot_irradiance_raytracing
from utils.process_radiometry_data import process_radiometry_file, calculate_correction_factor, plot_irradiance_radiometry
from utils.led import LST1_01F06_GRN1_00, LST1_01G01_UV01_00
from ICIW_Plots import cm2inch, make_square_subplots
from geometric_model import GeometricModel
from systematic_sampler import Sampler
from utils.grid import Grid

def create_irradiance_plots_SI():
    """
    Generates and saves irradiance comparison plots for multiple LED configurations using both radiometry and raytracing data in SI units (W/m²).

    For each LED configuration (UV1-3, GRN1-3, and best manual for each type), this function creates a figure with six subplots arranged in a 2x3 grid:

        - (0, 0): Geometric model (no reflector)
        - (0, 1): Raytracing (no reflector)
        - (0, 2): Radiometry (no reflector)
        - (1, 1): Raytracing (with reflector)
        - (1, 2): Radiometry (with reflector)
        - (1, 0): Empty (for layout)

    Each subplot visualizes the irradiance distribution on the target plane. All subplots share a unified colorbar for direct comparison.

    The function reads input data from specific directories, applies correction factors to radiometry data, and simulates the geometric model for each configuration. Figures are saved as SVG files in the 'figures' directory.

    :raises FileNotFoundError: If required data files or configuration JSONs are missing.
    :raises ValueError: If data files are malformed or contain unexpected values.
    :raises Exception: For errors during plotting or file saving.

    :side effects:
        - Creates and saves SVG figures in the 'figures' directory.
        - May create the 'figures' directory if it does not exist.
        - Displays figures using matplotlib's interactive window.

    :notes:
        - Assumes the presence of radiometry and raytracing data in the 'data/radiometry' and 'data/raytracing' directories, respectively.
        - Assumes LED configuration JSONs are present in 'sampled_configs/near_optimal_configs/<led_name>/<label>.json'.
        - Uses external utility functions for data processing and plotting.
        - The function does not return any value; it is intended for side effects (plotting and saving figures).

    :example:

        >>> create_irradiance_plots_SI()
        # Generates and saves irradiance comparison plots for all defined LED configurations.

    :Parameters:
        None

    :Returns:
        None

    """
    import numpy as np
    # Define directories
    radiometry_dir = os.path.join('data', 'radiometry')
    raytracing_dir = os.path.join('data', 'raytracing')

    led_configs = [
       {
           "label": "best_manual",
           "ref_file_radiometry": "best_manual_no_reflector.csv",          
           "led_obj": LST1_01G01_UV01_00,
       },
       {
           "label": "UV1",
           "ref_file_radiometry": "UV3_reflector.csv",
           "led_obj": LST1_01G01_UV01_00,
       },
        {
           "label": "UV2",
           "ref_file_radiometry": "UV3_reflector.csv",
           "led_obj": LST1_01G01_UV01_00,
       },
        {
           "label": "UV3",
           "ref_file_radiometry": "UV3_reflector.csv",
           "led_obj": LST1_01G01_UV01_00,
       },
        {
           "label": "best_manual",
           "ref_file_radiometry": "best_manual_no_reflector.csv",
           "led_obj": LST1_01F06_GRN1_00,
       },
       {
           "label": "GRN1",
           "ref_file_radiometry": "GRN3_reflector.csv",
           "led_obj": LST1_01F06_GRN1_00,
       },
       {
           "label": "GRN2",
           "ref_file_radiometry": "GRN3_reflector.csv",
           "led_obj": LST1_01F06_GRN1_00,
       },
       {
           "label": "GRN3",
           "ref_file_radiometry": "GRN3_reflector.csv",
           "led_obj": LST1_01F06_GRN1_00,
       },

    ]

    for config in led_configs:

        # if not "UV1" in config['label']:
        #     continue

        fig = plt.figure(figsize=(17*cm2inch, 10*cm2inch))
        axs = make_square_subplots(fig=fig,
                                   ax_layout=(2, 3),
                                   ax_width=4.5*cm2inch,
                                   h_sep=0.05*cm2inch,
                                   v_sep=0.05*cm2inch,
                                   xlabel='X (cm)',
                                   ylabel='Y (cm)',
                                   sharex=False,
                                   sharey=False
                                   )

        # filepaths
        radiometry_led_dir = os.path.join(radiometry_dir, config['led_obj'].name)
        raytracing_led_dir = os.path.join(raytracing_dir, config['led_obj'].name)

        # Raytracing No Reflector
        irr_ray_no_ref = process_raytracing_file(config['label']+ "_no_reflector.txt", raytracing_led_dir)
        # Radiometry No Reflector
        rad_no_ref_path = os.path.join(radiometry_led_dir, config['label'] + "_no_reflector.csv")
        if "best_manual" in config['label']:
            N_leds = 8
        else:
            N_leds = 16
        expected_value = config['led_obj'].total_power * N_leds
        correction_factor_no_ref = calculate_correction_factor(rad_no_ref_path, expected_value)
        irr_rad_no_ref, _, _ = process_radiometry_file(rad_no_ref_path, correction_factor=correction_factor_no_ref)
        # Raytracing Reflector
        irr_ray_ref = process_raytracing_file(config['label'] + "_reflector.txt", raytracing_led_dir)
        # Radiometry Reflector
        rad_ref_path = os.path.join(radiometry_led_dir, config['label'] + "_reflector.csv")
        correction_factor_ref = calculate_correction_factor(rad_ref_path, expected_value)
        irr_rad_ref, _, _ = process_radiometry_file(rad_ref_path, correction_factor=correction_factor_ref)

        # Set v_max for all plots
        v_max = max(
            irr_rad_no_ref.values.max(),
            irr_ray_no_ref['Irradiance'].max(),
            irr_rad_ref.values.max(),
            irr_ray_ref['Irradiance'].max()
        )
        v_max = round(v_max / 10) * 10

        extent = (
            irr_ray_no_ref['x'].min(), irr_ray_no_ref['x'].max(),
            irr_ray_no_ref['y'].min(), irr_ray_no_ref['y'].max()
        )

        # read config from sampled_configs
        G = Grid(width=extent[1] - extent[0], height=extent[3] - extent[2], step=2.5, side_space=1.5, top_bottom_space=2.)
        sampler = Sampler(G=G)
        if config['label'] == "best_manual":
            conf = [(-7.5,0), (7.5,0), (-10, -10), (10, 10), (-10, 10), (10, -10), (0,12.5), (0,-12.5)]
            height = 13
        else:
            conf, height = sampler.read_configuration_from_json(os.path.join('sampled_configs', 'near_optimal_configs', f'{config["led_obj"].name}', f'{config["label"]}.json'))
        model = GeometricModel(G=G, led=config['led_obj'])
        irr = model.simulate(config=conf, z=height)

        shrink_cbar = 1 # shrink colorbar for all plots
        y_title = 0

        # a) Geometric Model No Reflector
        model.plot_irradiance(axs[0, 0], irr, v_max=v_max, extent=extent, shrink_cbar=shrink_cbar, title="", y_title=y_title)
        # b) Raytracing No Reflector
        plot_irradiance_raytracing(axs[0, 1], irr_ray_no_ref, v_max=v_max, extent=extent, shrink_cbar=shrink_cbar, title="", y_title=y_title)
        # c) Radiometry No Reflector
        plot_irradiance_radiometry(axs[0, 2], irr_rad_no_ref, v_max=v_max, extent=extent, shrink_cbar=shrink_cbar, title="", y_title=y_title)
        # d) Raytracing Reflector
        plot_irradiance_raytracing(axs[1, 1], irr_ray_ref, v_max=v_max, extent=extent, shrink_cbar=shrink_cbar, title="", y_title=y_title)
        # e) Radiometry Reflector
        plot_irradiance_radiometry(axs[1, 2], irr_rad_ref, v_max=v_max, extent=extent, shrink_cbar=shrink_cbar, title="", y_title=y_title)

        # remove all colorbars and create unified colorbar
        colorbar_axes = []
        for ax in fig.axes:
            # Check if axis is a colorbar by its class name or if it has a colorbar label
            if 'colorbar' in ax.__class__.__name__.lower():
                colorbar_axes.append(ax)
            elif hasattr(ax, 'get_label') and ax.get_label() == '<colorbar>':
                colorbar_axes.append(ax)
            # Also check for colorbar axes by their position (usually on the right side)
            elif hasattr(ax, 'get_position'):
                pos = ax.get_position()
                if pos.x0 > 0.8:  # Heuristic: colorbars are often placed right of main axes
                    colorbar_axes.append(ax)
        for cax in colorbar_axes:
            fig.delaxes(cax)

        # Add a single colorbar for all subplots
        norm = Normalize(vmin=0, vmax=v_max)
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        # Place the colorbar to the right of all subplots
        cbar = fig.colorbar(sm, ax=axs, orientation='vertical', fraction=0.032, pad=0.04)
        cbar.set_label('Irradiance / W m$^{-2}$')

        for ax in [axs[0,1], axs[1,1], axs[0,2], axs[1,2]]:
            ax.set_ylabel('')

        for ax in [axs[0,1], axs[0,2]]:
            ax.set_xlabel('')
            # remove tick labels
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        axs[1,2].set_yticklabels([])
        axs[1,1].set_ylabel('$Y$ / cm')

        for ax in axs.flatten():
            ax.tick_params(direction='in', which='both')
        
                
                
        # axs[1,0] empty
        axs[1,0].axis('off')

            

        # plt.tight_layout()
        if not os.path.exists('figures'):
            os.makedirs('figures')
        if "best_manual" in config['label']:
            plt.savefig(os.path.join('figures', f'irradiance_plots_SI_best_manual_{config["led_obj"].name}.svg'), format='svg')
        else:
            plt.savefig(os.path.join('figures', f'irradiance_plots_SI_{config["label"]}.svg'), format='svg')
        plt.show()
