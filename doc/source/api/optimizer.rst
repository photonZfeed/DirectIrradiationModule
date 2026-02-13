Optimizer
========================

.. automodule:: optimizer
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
-----------------

    .. code-block:: python

      # Initialize the optimizer with an irradiance bound of 150.0 W/m²
      optimizer = Optimizer(irradiance_bound=150.0)
      
      # Load simulation results from file
      configs, z, means, max_vals, stds = optimizer.read_simulations('results\\examples\\example_MonteCarlo_simulation_results.npz')
      
      # Calculate metrics for all configurations
      ratio_to_bound, homogeneity = optimizer.calc_metrics(means, max_vals, stds)
      
      # Find the best configuration using product metric
      best_config, best_z, best_ratio, best_homogeneity, best_score = optimizer.find_best_configuration(
         configs, z, means, max_vals, stds)
      
      # Generate Pareto front by sweeping weighting factor t from 0 to 1
      t_values = np.linspace(0, 1, 21)
      pareto_configs, pareto_z, pareto_ratios, pareto_homogeneities, pareto_scores = optimizer.create_pareto_front(
         configs, z, means, max_vals, stds, t_values)
      
      # Visualize the Pareto front
      fig, ax = plt.subplots()
      optimizer.plot_pareto_front(
         ax,
         pareto_ratios,
         pareto_homogeneities,
         other_points={
               'Best Configuration': ('x', best_ratio, best_homogeneity),
         },
      )
      ax.legend()
      plt.tight_layout()
      plt.show()
      
      # Save optimization results to file
      save_optimization_results(
         filename="results\\examples\\example_optimization_results.npz",
         pareto_configs=pareto_configs,
         pareto_z=pareto_z,
         pareto_ratios=pareto_ratios,
         pareto_homogeneities=pareto_homogeneities,
         pareto_scores=pareto_scores,
         best_config=best_config,
         best_z=best_z,
         best_ratio=best_ratio,
         best_homogeneity=best_homogeneity,
         best_score=best_score
      )
      
      # Load optimization results from file
      results = read_optimization_results("results\\examples\\example_optimization_results.npz")