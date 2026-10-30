"""
Direct Irradiation Module Optimization - Usage Example

This script demonstrates the main functionality of the project:

1. Setting up the grid and the LED
2./3. Sampling LED configurations systematically or by Monte Carlo
4. Simulating the irradiance distribution of a single configuration
5. Simulating the irradiance of a batch of configurations
6. Finding the best configuration, i.e. the highest product of irradiance ratio and homogeneity
7. Creating the Pareto front for the multi-objective analysis
8. Analyzing and visualizing ray-tracing and radiometry data

All steps run one after the other and write their results to ``results/examples/`` (the example
files of the repository are overwritten). The script can be started from any folder on Windows,
Linux and macOS; it changes to the repository root itself, because the paths are relative to it.
Each figure is displayed until its window is closed.

See Also
--------
For a narrative walkthrough of each step with explanatory prose, see the tutorial
in the documentation: ``doc/source/tutorial.rst``.
"""

from utils.systematic_sampler import Sampler
from utils.systematic_sampler import save_configurations as save_sys_configurations
from utils.monte_carlo_sampler import MonteCarloSampler
from utils.monte_carlo_sampler import save_configurations as save_mc_configurations
from geometric_model import GeometricModel, save_results
from optimizer import Optimizer, save_optimization_results

from utils.grid import Grid
from utils.led import LED

import numpy as np
from matplotlib import pyplot as plt
import os


def main():
    """Demonstration of the radiation field optimization workflow."""
    
    # ===== STEP 1: Define the grid and LED parameters =====
    print("Step 1: Initializing grid and LED configuration...")
    G = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
    led = LED(
        name="LST1-01G01-UV01-00",
        optical_power=0.875,  # @ 500 mA in W
        intensity_data_file="utils/angle_intensity_data/LST1-01G01-UV01-00.csv"
    )
    
    # ===== STEP 2: Generate configurations using systematic sampler =====
    print("\nStep 2: Generating configurations with systematic sampler...")
    # Create a systematic sampler for 8 LEDs at height 13 cm
    sampler = Sampler(G=G, led_count=8, height=13)
    # Generate and save all symmetric configurations of 8 LEDs
    samples_sys = sampler.generate_all_configurations()
    save_sys_configurations(samples_sys, sampler, "results/examples/example_systematic_8_leds_13_cm.npz")
    
    # ===== STEP 3: Generate configurations using Monte Carlo sampler =====
    print("\nStep 3: Generating configurations with Monte Carlo sampler...")
    # Create a Monte Carlo sampler for flexible LED counts and heights
    mc_sampler = MonteCarloSampler(
        G=G, led_min=1, led_max=16, 
        height_min=5, height_max=15, height_step=1, 
        seed=42, sampling_mode="uniform", verbose=False
    )
    # Generate and save 10,000 distinct Monte Carlo samples
    num_samples = 10000
    samples_mc = mc_sampler.generate_unique_samples(num_samples=num_samples, use_parallel=True)
    save_mc_configurations(samples_mc, mc_sampler, 
        f"results/examples/example_MonteCarlo_{mc_sampler.led_min}-{mc_sampler.led_max}_leds_"
        f"{mc_sampler.height_min}-{mc_sampler.height_max}_cm_{num_samples}_samples.npz")
    
    # ===== STEP 4: Simulate irradiance for a single configuration =====
    print("\nStep 4: Simulating irradiance distribution of an example configuration...")
    model = GeometricModel(G=G, led=led)
    sample_config, z = sampler.read_configuration_from_json("results/sampled_configs/near_optimal_configs/LST1-01G01-UV01-00/UV3.json")
    irr = model.simulate(config=sample_config, z=z)
    
    # Plot the result
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle("Irradiance Distribution of UV3 (Geometric Model)", fontsize=14)
    model.plot_irradiance(ax=ax, irr=irr)
    plt.show()
    plt.close(fig)
    
    # ===== STEP 5: Batch simulate configurations =====
    print("\nStep 5: Batch simulating all configurations...")
    # Simulate all configurations of the Monte Carlo samples
    configs, z = model.read_configurations(
        f"results/examples/example_MonteCarlo_{mc_sampler.led_min}-{mc_sampler.led_max}_leds_"
        f"{mc_sampler.height_min}-{mc_sampler.height_max}_cm_{num_samples}_samples.npz")
    sim_results = model.simulate_batch(configurations=configs, z=z, batch_size=500)
    save_results(sim_results, model, "results/examples/example_MonteCarlo_simulation_results.npz")
    
    # ===== STEP 6: Optimize for best configuration =====
    print("\nStep 6: Finding optimal configuration...")
    irradiance_bound = 150  # upper limit of the local irradiance in W/m^2
    optimizer = Optimizer(irradiance_bound=irradiance_bound)
    configs, z, means, maxs, stds = optimizer.read_simulations(
        "results/examples/example_MonteCarlo_simulation_results.npz")
    best_config, best_z, best_ratio, best_homogeneity, best_score = \
        optimizer.find_best_configuration(configs, z, means, maxs, stds)
    print(f"Best configuration: {best_config}")
    print(f"Best height: {best_z} cm")
    
    # ===== STEP 7: Create and visualize Pareto front =====
    print("\nStep 7: Creating Pareto front analysis...")
    # Pareto front over 100 weighting factors t (uses the simulation results of step 6)
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

    # ===== STEP 8: Analyze ray-tracing and radiometry data =====
    print("\nStep 8: Analyzing ray-tracing and radiometry data...")
    from analysis.process_radiometry_data import calculate_correction_factor, process_radiometry_file, plot_irradiance_radiometry
    from analysis.process_raytracing_data import process_raytracing_file, plot_irradiance_raytracing
    
    # Configuration parameters
    N_led = 16  # number of LEDs in the module
    
    # Process radiometry data
    radiometry_file = f"./results/radiometry/{led.name}/UV3_no_reflector.csv"
    reference_file = (f"./results/radiometry/{led.name}/UV3_reflector.csv")
    expected_value = N_led * led.total_power  # emitted radiant flux of all LEDs in W, to which the reference scan is scaled
    correction_factor = calculate_correction_factor(reference_file, expected_value)
    irr, mean_irradiance, homogeneity = process_radiometry_file(radiometry_file, correction_factor)
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle(f"Radiometry Data Analysis for {os.path.basename(radiometry_file)}", fontsize=14)
    plot_irradiance_radiometry(ax=ax, irr=irr)
    plt.show()
    plt.close(fig)
    
    # Process ray-tracing data
    filename = "UV3_no_reflector.txt"
    raytracing_path = f"./results/raytracing/{led.name}/"
    data = process_raytracing_file(filename=filename, path=raytracing_path)
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle(f"Ray-Tracing Data Analysis for {filename}", fontsize=14)
    plot_irradiance_raytracing(ax=ax, data=data)
    plt.show()
    plt.close(fig)
    
    print("\nDemo complete!")

if __name__ == "__main__":
    # the paths of this script are relative to the repository root
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    main()
