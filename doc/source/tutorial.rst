Getting Started: Full Optimization Workflow
============================================

This tutorial walks through the complete workflow for designing and optimizing a direct irradiation module, from defining the physical geometry to selecting the best LED configuration according to the defined objective function. The code shown here is collected in ``usage_example.py``, which runs all steps one after the other. All paths are relative to the root of the repository; run the code from there.


Setting Up the Grid and LED
---------------------------

The :class:`~utils.grid.Grid` object defines the detection plane (here 33 cm × 34 cm) and the possible LED mounting positions on it: a square grid with a pitch of ``step`` = 2.5 cm, kept ``side_space`` and ``top_bottom_space`` away from the edges, which gives the 13 × 13 = 169 mounting positions of the module. The resolution at which the irradiance is evaluated is set later by the geometric model (``resolution_xy``, 0.5 cm by default). The :class:`~utils.led.LED` object encodes the physical properties of a single LED, including its optical output power and its relative angular intensity distribution loaded from a CSV file.

.. code-block:: python

    from utils.grid import Grid
    from utils.led import LED

    G = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
    led = LED(
        name="LST1-01G01-UV01-00",
        optical_power=0.875,  # @ 500 mA in W
        intensity_data_file="utils/angle_intensity_data/LST1-01G01-UV01-00.csv"
    )


Generating Configurations: Systematic Sampler
----------------------------------------------

The :class:`~utils.systematic_sampler.Sampler` generates all symmetrically valid arrangements of a fixed number of LEDs on the grid at a single height. This exhaustive enumeration is practical when the LED count is small, because the symmetry constraints keep the search space manageable. The resulting configurations and their metadata are saved to a compressed ``.npz`` file for later use.

.. code-block:: python

    from utils.systematic_sampler import Sampler
    from utils.systematic_sampler import save_configurations as save_sys_configurations

    sampler = Sampler(G=G, led_count=8, height=13)
    samples_sys = sampler.generate_all_configurations()
    save_sys_configurations(samples_sys, sampler,
        "results/examples/example_systematic_8_leds_13_cm.npz")


Generating Configurations: Monte Carlo Sampler
-----------------------------------------------

When the design space is too large for exhaustive enumeration, the :class:`~utils.monte_carlo_sampler.MonteCarloSampler` draws random unique configurations from the full space. The ``seed`` parameter ensures reproducibility: :meth:`~utils.monte_carlo_sampler.MonteCarloSampler.generate_unique_samples` draws the configurations one after the other (``use_parallel`` only takes effect together with ``target_samples_per_led``), so the result depends on the seed only. With ``seed=42`` and 5 million samples, it reproduces the sample of the publication. The file path is constructed to encode the key sampling parameters so that it is self-documenting.

.. code-block:: python

    from utils.monte_carlo_sampler import MonteCarloSampler
    from utils.monte_carlo_sampler import save_configurations as save_mc_configurations

    mc_sampler = MonteCarloSampler(
        G=G, led_min=1, led_max=16,
        height_min=5, height_max=15, height_step=1,
        seed=42, sampling_mode="uniform", verbose=False
    )
    num_samples = 10000
    samples_mc = mc_sampler.generate_unique_samples(num_samples=num_samples, use_parallel=True)
    save_mc_configurations(samples_mc, mc_sampler,
        f"results/examples/example_MonteCarlo_{mc_sampler.led_min}-{mc_sampler.led_max}_leds_"
        f"{mc_sampler.height_min}-{mc_sampler.height_max}_cm_{num_samples}_samples.npz")


Simulating Irradiance for a Single Configuration
-------------------------------------------------

The :class:`~geometric_model.GeometricModel` takes the grid and LED objects and computes the irradiance distribution at each grid point. The result is a 2-D array of irradiance values (W m\ :sup:`-2`) that can be visualized immediately.

.. code-block:: python

    from geometric_model import GeometricModel
    import matplotlib.pyplot as plt

    model = GeometricModel(G=G, led=led)
    sample_config, z = sampler.read_configuration_from_json(
        "results/sampled_configs/near_optimal_configs/LST1-01G01-UV01-00/UV3.json")
    irr = model.simulate(config=sample_config, z=z)

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle("Irradiance Distribution of UV3 (Geometric Model)", fontsize=14)
    model.plot_irradiance(ax=ax, irr=irr)
    plt.show()
    plt.close(fig)


Batch Simulating a Configuration Set
-------------------------------------

Simulating thousands of configurations one at a time would be slow. :meth:`~geometric_model.GeometricModel.simulate_batch` processes configurations in parallel using a ``ProcessPoolExecutor`` and accepts a ``batch_size`` argument to control memory usage. The worker processes are started with ``spawn`` on every operating system, so code that calls it from a script must be guarded by ``if __name__ == "__main__":``, as in ``usage_example.py``. It reads the configurations saved in the previous Monte Carlo step and returns a structured result object that is then written to disk with :func:`~geometric_model.save_results`.

.. code-block:: python

    from geometric_model import save_results

    configs, z = model.read_configurations(
        f"results/examples/example_MonteCarlo_{mc_sampler.led_min}-{mc_sampler.led_max}_leds_"
        f"{mc_sampler.height_min}-{mc_sampler.height_max}_cm_{num_samples}_samples.npz")
    sim_results = model.simulate_batch(configurations=configs, z=z, batch_size=500)
    save_results(sim_results, model, "results/examples/example_MonteCarlo_simulation_results.npz")


Finding the Best Configuration
--------------------------------

The :class:`~optimizer.Optimizer` is initialized with an irradiance bound — the maximum permissible irradiance in W m\ :sup:`-2`. Any configuration whose peak irradiance exceeds this bound is excluded from consideration. :meth:`~optimizer.Optimizer.find_best_configuration` selects the configuration that maximizes the product of the mean-irradiance-to-bound ratio and the homogeneity index, providing a single-figure-of-merit ranking without requiring manual weight tuning.

.. code-block:: python

    from optimizer import Optimizer

    irradiance_bound = 150  # upper limit of the local irradiance in W m⁻²
    optimizer = Optimizer(irradiance_bound=irradiance_bound)
    configs, z, means, maxs, stds = optimizer.read_simulations(
        "results/examples/example_MonteCarlo_simulation_results.npz")
    best_config, best_z, best_ratio, best_homogeneity, best_score = \
        optimizer.find_best_configuration(configs, z, means, maxs, stds)
    print(f"Best configuration: {best_config}")
    print(f"Best height: {best_z} cm")


Multi-Objective Analysis: Pareto Front
----------------------------------------

In practice, mean irradiance and homogeneity are competing objectives: increasing one often degrades the other. :meth:`~optimizer.Optimizer.create_pareto_front` sweeps a weighting factor ``t`` from 0 (homogeneity only) to 1 (mean irradiance only), collecting the Pareto-optimal configurations at each trade-off point. The resulting front is visualized with :meth:`~optimizer.Optimizer.plot_pareto_front`, and the full results are saved for later retrieval.

.. code-block:: python

    from optimizer import save_optimization_results
    import numpy as np

    t_values = np.linspace(0, 1, 100)
    pareto_configs, pareto_z, pareto_ratios, pareto_homogeneities, pareto_scores = \
        optimizer.create_pareto_front(t_values=t_values, configurations=configs, z=z,
                                      means=means, max=maxs, stds=stds)
    fig, ax = plt.subplots(figsize=(8, 6))
    optimizer.plot_pareto_front(ax=ax, pareto_ratios=pareto_ratios,
                                pareto_homogeneities=pareto_homogeneities)
    plt.show()

    save_optimization_results(
        "results/examples/example_optimization_results.npz",
        pareto_configs, pareto_z, pareto_ratios, pareto_homogeneities, pareto_scores,
        best_config, best_z, best_ratio, best_homogeneity, best_score
    )


Analyzing Radiometry and Ray-Tracing Data
-----------------------------------------

Once a physical prototype is built, measured radiometry data and independent ray-tracing simulations can be compared with the predictions of the geometric model. The radiometric scans are not calibrated in absolute units: :func:`~analysis.process_radiometry_data.calculate_correction_factor` scales a reference scan (here UV3 with reflectors) so that its integral equals the radiant flux emitted by all LEDs according to the datasheet, and the same factor is applied to all scans of the group. This assumes that the entire emitted flux reaches the scanned area; :func:`analysis.uncertainty.radiometry_corrections` estimates the share that actually reaches it (Supporting Information Section S10). Ray-tracing data is processed analogously with :func:`~analysis.process_raytracing_data.process_raytracing_file`. Both datasets are visualized with the corresponding plot functions.

.. code-block:: python

    import os

    from analysis.process_radiometry_data import (
        calculate_correction_factor, process_radiometry_file, plot_irradiance_radiometry
    )
    from analysis.process_raytracing_data import process_raytracing_file, plot_irradiance_raytracing

    N_led = 16  # number of LEDs in the module
    radiometry_file = f"./results/radiometry/{led.name}/UV3_no_reflector.csv"
    reference_file  = f"./results/radiometry/{led.name}/UV3_reflector.csv"
    expected_value  = N_led * led.total_power  # emitted radiant flux of all LEDs in W

    correction_factor = calculate_correction_factor(reference_file, expected_value)
    irr, mean_irradiance, homogeneity = process_radiometry_file(radiometry_file, correction_factor)
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle(f"Radiometry Data Analysis for {os.path.basename(radiometry_file)}", fontsize=14)
    plot_irradiance_radiometry(ax=ax, irr=irr)
    plt.show()
    plt.close(fig)

    filename = "UV3_no_reflector.txt"
    raytracing_path = f"./results/raytracing/{led.name}/"
    data = process_raytracing_file(filename=filename, path=raytracing_path)
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle(f"Raytracing Data Analysis for {filename}", fontsize=14)
    plot_irradiance_raytracing(ax=ax, data=data)
    plt.show()
    plt.close(fig)


Comparing the Methods and Reproducing the Publication
-----------------------------------------------------

:mod:`analysis.compare_methods` compares the geometric model, ray tracing and radiometry point by
point on a common 0.5 cm grid and reports the bias, the root mean square error and the largest local
deviation of every pair of methods (Supporting Information Section S6). The further analyses of the
Supporting Information and all figures and tables of the publication are generated by
``create_manuscript_plots_tables.py``; :doc:`reproduce` lists the outputs and the functions behind them.

.. code-block:: python

    from analysis.compare_methods import compare_all

    residuals, statistics = compare_all("results/")
    print(statistics[["LED", "configuration", "comparison", "RMSE / %"]].round(1))
