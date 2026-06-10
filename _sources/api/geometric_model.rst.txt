Geometric Model 
========================

.. automodule:: geometric_model
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
-----------------

    .. code-block:: python

        from utils.grid import Grid
        from utils.led import LST1_01G01_UV01_00
        from geometric_model import GeometricModel, save_results

        # Define grid and LED
        G = Grid(width=20, height=20, step=2.5)
        led = LST1_01G01_UV01_00
        model = GeometricModel(G, led, resolution_xy=0.5)

        # Define LED configuration and height
        config = [(-5, 0), (5, 0)]  # two LEDs at (-5,0) and (5,0) cm
        z = 10  # height above grid in cm

        # Simulate irradiance
        irr = model.simulate(config, z)

        # Visualize
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        model.plot_irradiance(ax, irr)
        plt.show()

        # simulate batch of configurations
        configs, z = model.read_configurations("results\\examples\\example_MonteCarlo_1-16_leds_5-15_cm_10000_samples.npz") # load configurations and heights from file
        sim_results = model.simulate_batch(configurations=configs, z=z, batch_size=500)
        save_results(sim_results, model, "results/examples/example_MonteCarlo_simulation_results.npz")