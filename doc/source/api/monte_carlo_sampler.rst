Monte-Carlo Sampler
========================

.. automodule:: monte_carlo_sampler
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
----------------

   .. code-block:: python
    
    from utils.grid import Grid
    G = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)

    # Create a Monte Carlo sampler for flexible LED counts and heights
    mc_sampler = MonteCarloSampler(
        G=G, led_min=1, led_max=16, 
        height_min=5, height_max=15, height_step=1, 
        seed=42, sampling_mode="uniform", verbose=False
    )
    # Generate and save 10000 Monte Carlo samples
    num_samples = 10000
    samples_mc = mc_sampler.generate_unique_samples(num_samples=num_samples, use_parallel=True)
    save_mc_configurations(samples_mc, mc_sampler, 
        f"results/examples/example_MonteCarlo_{mc_sampler.led_min}-{mc_sampler.led_max}_leds_"
        f"{mc_sampler.height_min}-{mc_sampler.height_max}_cm_{num_samples}_samples.npz")