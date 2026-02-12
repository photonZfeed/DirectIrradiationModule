import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from geometric_model import GeometricModel
from systematic_sampler import Sampler
from utils.grid import Grid
from utils.led import LST1_01G01_UV01_00, LST1_01F06_GRN1_00
from utils.process_radiometry_indirect import process_radiometry_data
from ICIW_Plots import make_square_subplots, cm2inch
from optimizer import Optimizer, read_optimization_results, save_optimization_results

def process_led(ax: Axes, led, optimizer, grid, t_values, results_path, sim_filename, title, load_previous: bool = False, find_raytracing_candidates: bool = False, plot_near_optimal_configs: bool = False):
    """
    Generates and visualizes the Pareto front for a specified LED configuration, including optional export and plotting of near-optimal configurations, and overlays additional reference points for comparison.

    This function reads or computes optimization results for a given LED and simulation setup, plots the Pareto front of mean irradiance ratio versus homogeneity index, and optionally exports near-optimal configurations for further raytracing analysis. It also supports plotting insets of selected configurations and overlays reference points such as indirect irradiation and manually optimized configurations.

    :param ax: Matplotlib axes object on which to plot the Pareto front and configuration insets.
        Type: :class:`matplotlib.axes.Axes`
    :param led: LED object to be optimized. Must have a ``name`` attribute for file path construction.
        Type: Any
    :param optimizer: Optimizer instance used for Pareto front calculation and metric evaluation.
        Type: :class:`Optimizer`
    :param grid: Grid object defining the spatial simulation domain.
        Type: :class:`Grid`
    :param t_values: Array of t values (unitless, shape (N,)) for Pareto front calculation.
        Type: :class:`numpy.ndarray`
    :param results_path: Path to the directory containing simulation and result files.
        Type: str
    :param sim_filename: Filename of the simulation results (e.g., ``simulation_results_...npz``).
        Type: str
    :param title: Title for the Pareto plot.
        Type: str
    :param load_previous: If True, loads previous optimization results from disk if available; otherwise, computes and saves new results. Default is False.
        Type: bool, optional
    :param find_raytracing_candidates: If True, exports near-optimal configurations (around the best found) as JSON files for raytracing. Default is False.
        Type: bool, optional
    :param plot_near_optimal_configs: If True, plots all near-optimal configurations from exported JSON files as insets; otherwise, plots only start, end, and best configurations. Default is False.
        Type: bool, optional

    :returns: The axes object with the Pareto plot and configuration insets.
    :rtype: :class:`matplotlib.axes.Axes`

    :raises FileNotFoundError: If required simulation or radiometry data files are missing.
    :raises ValueError: If the LED object does not have a ``name`` attribute or if configuration files are malformed.

    :side effects:
        - May write optimization and configuration JSON files to disk if ``find_raytracing_candidates`` is True or if new optimization results are computed.
        - May print status messages to stdout.

    .. note::
        - Assumes the existence of a compatible directory structure under ``results_path`` for reading/writing results and configuration files.
        - The function creates insets on the provided axes for visualizing selected configurations.
        - The optimizer is re-instantiated within the function, which may override the passed-in optimizer.

    **Usage Example**::

        import matplotlib.pyplot as plt
        from utils.led import LST1_01G01_UV01_00
        from optimizer import Optimizer
        from utils.grid import Grid
        fig, ax = plt.subplots()
        grid = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
        optimizer = Optimizer(irradiance_bound=150.0)
        t_values = np.linspace(0, 1, 100)
        process_led(ax, LST1_01G01_UV01_00, optimizer, grid, t_values, 'results/', 'simulation_results.npz', '365 nm LEDs')

    """
    # Create geometric model instance and optimizer
    model = GeometricModel(grid, led, resolution_xy=0.5)
    optimizer = Optimizer(irradiance_bound=150.0)

    # Prepare file paths
    opt_resultsfile = sim_filename.replace("simulation_results", "optimization_results")
    opt_resultpath = os.path.join(results_path, f"geometric/{led.name}/{opt_resultsfile}")
    sim_resultpath = os.path.join(results_path, f"geometric/{led.name}/{sim_filename}")

    # Helper to read or compute optimization results
    def get_optimization_results() -> tuple:
        """
        Read or compute optimization results for the given simulation.

        Returns
        -------
        tuple
            Tuple containing Pareto configs, z, ratios, homogeneities, scores, best config, best z, best ratio, best homogeneity, best score.
        """
        if load_previous and os.path.exists(opt_resultpath):
            return read_optimization_results(opt_resultpath)
        configs, z, means, maxs, stds = optimizer.read_simulations(sim_resultpath)
        print(f"Read {len(configs)} configurations for LED {led.name}")
        pareto_configs, pareto_z, pareto_ratios, pareto_homogeneities, pareto_scores = optimizer.create_pareto_front(configs, z, means, maxs, stds, t_values)
        best_config, best_z, best_ratio, best_homogeneity, best_score = optimizer.find_best_configuration(configs, z, means, maxs, stds)
        if not load_previous:
            save_optimization_results(
                opt_resultpath,
                pareto_configs, pareto_z, pareto_ratios, pareto_homogeneities, pareto_scores,
                best_config, best_z, best_ratio, best_homogeneity, best_score
            )
        return (pareto_configs, pareto_z, pareto_ratios, pareto_homogeneities, pareto_scores,
                best_config, best_z, best_ratio, best_homogeneity, best_score)

    (
        pareto_configs, pareto_z, pareto_ratios, pareto_homogeneities, pareto_scores,
        best_config, best_z, best_ratio, best_homogeneity, best_score
    ) = get_optimization_results()

    # Export near-optimal configs for raytracing if requested
    if find_raytracing_candidates:
        best_index_in_pareto = next((i for i, (pcfg, pz) in enumerate(zip(pareto_configs, pareto_z))
                                    if np.array_equal(pcfg, best_config) and pz == best_z), None)
        if best_index_in_pareto is not None:
            export_dir = os.path.join(results_path, "sampled_configs/near_optimal_configs", led.name)
            os.makedirs(export_dir, exist_ok=True)
            sampler = Sampler(grid)
            for idx in range(max(0, best_index_in_pareto-2), min(len(pareto_configs), best_index_in_pareto+1)):
                configuration = pareto_configs[idx]
                height = pareto_z[idx]
                new_idx = idx - max(0, best_index_in_pareto-2) + 1  # 1,2,3
                suffix = "UV" if led == LST1_01G01_UV01_00 else "GRN"
                filepath = os.path.join(export_dir, f"{suffix}{new_idx}.json")
                sampler.save_configuration_to_json(
                    configuration=configuration,
                    height=height,
                    filename=filepath,
                    save_fig_type='svg'
                )

    # Prepare insets for configuration plots
    Sampler1 = Sampler(grid, 8)
    inset_positions = [
        (0.05, 0.07, 0.25, 0.25),  # upper left
        (0.35, 0.4, 0.25, 0.25),   # center
        (0.6, 0.07, 0.25, 0.25),   # lower right
    ]
    inset_axes_list = [ax.inset_axes(pos) for pos in inset_positions]

    if plot_near_optimal_configs:
        # Plot all near-optimal configs from JSON
        configs_to_plot, heights_to_plot = [], []
        near_optimal_configs_path = os.path.join(results_path, "sampled_configs/near_optimal_configs", led.name)
        for file in sorted(os.listdir(near_optimal_configs_path)):
            if file.endswith(".json"):
                config, height = Sampler1.read_configuration_from_json(os.path.join(near_optimal_configs_path, file))
                configs_to_plot.append(config)
                heights_to_plot.append(height)
        for iax, config, height in zip(inset_axes_list, configs_to_plot, heights_to_plot):
            Sampler1.plot_configuration_small(iax, config)
            iax.set_title(f"{len(config)} LEDs, {int(height)} cm", fontsize=10, y=-0.35)
    else:
        # Plot only start, end, and best configuration
        configs_to_plot = [pareto_configs[0], pareto_configs[-1], best_config]
        heights_to_plot = [pareto_z[0], pareto_z[-1], best_z]
        for iax, config, height in zip(inset_axes_list, configs_to_plot, heights_to_plot):
            Sampler1.plot_configuration_small(iax, config)

    # Calculate metrics for manual config
    best_manual_config = [(-7.5,0), (7.5,0), (-10, -10), (10, 10), (-10, 10), (10, -10), (0,12.5), (0,-12.5)]
    best_manual_irr = model.simulate(best_manual_config, 13.0)
    best_manual_mean = np.mean(best_manual_irr)
    best_manual_std = np.std(best_manual_irr)
    best_manual_ratio, best_manual_homogeneity = optimizer.calc_metrics(
        np.array([best_manual_mean]), np.array([np.max(best_manual_irr)]), np.array([best_manual_std])
    )

    # Read indirect irradiation data
    _, E_mean_indirect, H_indirect = process_radiometry_data(
        os.path.join(results_path, f"radiometry/{led.name}/indirect_irradiation.csv"))

    # Prepare other points for Pareto plot
    other_points = {
        "Indirect Irradiation": ("d", E_mean_indirect/optimizer.irradiance_bound, H_indirect),
        "Manually Optimized Configuration": ("s", best_manual_ratio, best_manual_homogeneity),
        "Best Configuration": ("X", best_ratio, best_homogeneity),
    }

    ax = optimizer.plot_pareto_front(ax, pareto_ratios, pareto_homogeneities, other_points=other_points, title=title)
    return ax

def create_pareto_1_16_leds_5_15cm():
    """
    Generates and visualizes the Pareto front for LED arrays ranging from 1 to 16 LEDs positioned at heights between 5 cm and 15 cm.

    This function sets up the simulation grid, optimizer, and LED types, then loads or computes the Pareto front for each LED type using precomputed simulation results. The resulting Pareto plots display the trade-off between mean irradiance ratio and homogeneity index for each configuration. Insets show selected LED arrangements. The function creates a single figure with two subplots (one for each LED type) and displays the result.

    :param None:
        This function does not take any arguments.

    :returns: None
        The function creates and displays a matplotlib figure but does not return any value.

    :raises FileNotFoundError:
        If required simulation or result files are missing.

    :side effects:
        - Displays a matplotlib figure window with the Pareto plots and configuration insets.
        - May print status messages to stdout.
        - Optionally saves the figure to disk if the corresponding line is uncommented.

    .. note::
        - Assumes the existence of simulation and result files in the specified directory structure.
        - Uses precomputed simulation data for efficiency.
        - The function is intended for interactive use or for generating publication-quality figures.

    **Parameters used internally:**

    - ``grid`` (:class:`Grid`): Defines the spatial simulation domain (width=33 cm, height=34 cm, step=2.5 cm, side_space=1.5 cm, top_bottom_space=2.0 cm).
    - ``optimizer`` (:class:`Optimizer`): Optimizer instance with irradiance bound set to 150.0 mW/cm².
    - ``leds`` (list): List of tuples with LED objects and subplot titles.
    - ``t_values`` (:class:`numpy.ndarray`): Array of t values (unitless, shape (100,)) for Pareto front calculation.
    - ``results_path`` (str): Path to the directory containing simulation and result files (default: "results/").
    - ``sim_filename`` (str): Filename of the simulation results (default: "simulation_results_MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz").
    - ``figure_path`` (str): Path to the directory for saving figures (default: "figures/").

    **Usage Example:**

    .. code-block:: python

        from pareto_plots_H_vs_E_ratio import create_pareto_1_16_leds_5_15cm
        create_pareto_1_16_leds_5_15cm()

    """
    grid = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
    optimizer = Optimizer(irradiance_bound=150.0)
    leds  = [
        (LST1_01G01_UV01_00, "a) 365 nm LEDs"),
        (LST1_01F06_GRN1_00, "b) 530 nm LEDs"),
    ]
    t_values = np.linspace(0, 1, 100)
    results_path = "results/"
    sim_filename = "simulation_results_MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz"
    figure_path = "figures/"
    fig = plt.figure(figsize=(19*cm2inch, 12*cm2inch))
    axs = make_square_subplots(
        fig=fig,
        ax_width=8 * cm2inch,
        ax_layout=(1, 2),
        h_sep=2 * cm2inch,
        v_sep=0 * cm2inch,
        sharex=False,
        sharey=True,
        sharelabel=False,
        xlabel="$E_\\mathrm{mean}\\,/\\,E_\\mathrm{bound}$ / 1",
        ylabel="Homogeneity Index / 1"
    )
    for i, (led, title) in enumerate(leds):
        ax = axs[0, i]
        process_led(ax, led, optimizer, grid, t_values, results_path, sim_filename, title, load_previous=True, plot_near_optimal_configs=True)
    # Single legend for all subplots
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 0.94))
    # plt.savefig(os.path.join(figure_path, f"pareto_plots_{results_path.split('/')[-2]}.svg"))
    plt.show()

def create_pareto_8_leds_13cm():
    """
    Generates and visualizes the Pareto front for an array of 8 LEDs positioned at a height of 13 cm.

    This function sets up the simulation grid, optimizer, and LED types, then loads or computes the Pareto front for each LED type using precomputed simulation results. The resulting Pareto plots display the trade-off between mean irradiance ratio and homogeneity index for each configuration. Insets show selected LED arrangements. The function creates a single figure with two subplots (one for each LED type) and displays the result.

    :param None:
        This function does not take any arguments.

    :returns: None
        The function creates and displays a matplotlib figure but does not return any value.

    :raises FileNotFoundError:
        If required simulation or result files are missing.

    :side effects:
        - Displays a matplotlib figure window with the Pareto plots and configuration insets.
        - May print status messages to stdout.
        - Optionally saves the figure to disk if the corresponding line is uncommented.

    .. note::
        - Assumes the existence of simulation and result files in the specified directory structure.
        - Uses precomputed simulation data for efficiency.
        - The function is intended for interactive use or for generating publication-quality figures.

    **Parameters used internally:**

    - ``grid`` (:class:`Grid`): Defines the spatial simulation domain (width=33 cm, height=34 cm, step=2.5 cm, side_space=1.5 cm, top_bottom_space=2.0 cm).
    - ``optimizer`` (:class:`Optimizer`): Optimizer instance with irradiance bound set to 150.0 mW/cm².
    - ``leds`` (list): List of tuples with LED objects and subplot titles.
    - ``t_values`` (:class:`numpy.ndarray`): Array of t values (unitless, shape (100,)) for Pareto front calculation.
    - ``results_path`` (str): Path to the directory containing simulation and result files (default: "results/").
    - ``sim_filename`` (str): Filename of the simulation results (default: "simulation_results_systematic_8_leds_13_cm_482931_samples.npz").
    - ``figure_path`` (str): Path to the directory for saving figures (default: "figures/").

    **Usage Example:**

    .. code-block:: python

        from pareto_plots_H_vs_E_ratio import create_pareto_8_leds_13cm
        create_pareto_8_leds_13cm()

    """
    grid: Grid = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
    optimizer: Optimizer = Optimizer(irradiance_bound=150.0)
    leds: list = [
        (LST1_01G01_UV01_00, "a) 365 nm LEDs"),
        (LST1_01F06_GRN1_00, "b) 530 nm LEDs"),
    ]
    t_values: np.ndarray = np.linspace(0, 1, 100)
    results_path: str = "results/"
    sim_filename: str = "simulation_results_systematic_8_leds_13_cm_482931_samples.npz"
    figure_path: str = "figures/"
    fig = plt.figure(figsize=(19*cm2inch, 12*cm2inch))
    axs = make_square_subplots(
        fig=fig,
        ax_width=8 * cm2inch,
        ax_layout=(1, 2),
        h_sep=2 * cm2inch,
        v_sep=0 * cm2inch,
        sharex=False,
        sharey=True,
        sharelabel=False,
        xlabel="$E_\\mathrm{mean}\\,/\\,E_\\mathrm{bound}$ / 1",
        ylabel="Homogeneity Index / 1"
    )
    for i, (led, title) in enumerate(leds):
        ax = axs[0, i]
        process_led(ax, led, optimizer, grid, t_values, results_path, sim_filename, title, load_previous=True, plot_near_optimal_configs=False)
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 0.94))
    plt.savefig(os.path.join(figure_path, f"pareto_plots_{results_path.split('/')[-2]}.svg"))
    plt.show()