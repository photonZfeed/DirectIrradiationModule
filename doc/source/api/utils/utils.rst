LED 
====

.. automodule:: utils.led
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
---------------

.. code-block:: python

    led = LED(name='LST1-01G01-UV01-00', optical_power=0.875,  intensity_data_file='utils/angle_intensity_data/LST1-01G01-UV01-00.csv')
    abs_intensity = led.calc_radiant_intensity()
    print(abs_intensity)
    
Grid
=====

.. automodule:: utils.grid
    :members:
    :undoc-members:
    :show-inheritance:

Example Usage
--------------

.. code-block:: python

    grid = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.0)
    print(grid.grid_points)

Raytracing Data Processing
==========================

.. automodule:: utils.process_raytracing_data
    :members:
    :undoc-members:
    :show-inheritance:

Example Usage
-------------

.. code-block:: python

    # Load raytracing data
    df = process_raytracing_file('best_manual_no_reflector.txt', 'results/raytracing/LST1-01G01-UV01-00/')
    
    # Calculate metrics
    mean_irr, homogeneity = calc_metrics(df['Irradiance'])
    print(f'Mean Irradiance: {mean_irr:.2f} W/m²')
    print(f'Homogeneity Index: {homogeneity:.4f}')
    
    # Create visualization
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_irradiance_raytracing(ax, df, title='Irradiance Distribution')
    plt.show()

Radiometry Data Processing
==========================

Direct Irradiation Module
--------------------------
.. automodule:: utils.process_radiometry_data
    :members:
    :undoc-members:
    :show-inheritance:

Example Usage
^^^^^^^^^^^^^^

    .. code-block:: python

        # Calculate correction factor from a reference measurement
        factor = calculate_correction_factor('results/radiometry/LST1-01G01-UV01-00/UV3_reflector.csv', expected_value=16*0.875)
        
        # Process all candidate files in a folder
        results = process_all_cand_files('results/radiometry/LST1-01G01-UV01-00/', correction_factor=factor)
        
        # Visualize first result
        first_result = list(results.values())[0]
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_irradiance_radiometry(ax, first_result, title='Corrected Irradiance Distribution')
        plt.show()


Indirect Irradiation Module
----------------------------
.. automodule:: utils.process_radiometry_indirect
    :members:
    :undoc-members:
    :show-inheritance:

Example Usage
^^^^^^^^^^^^^^

    .. code-block:: python

        # calculate mean irradiance and homogeneity index from an indirect irradiation measurement
        irr, mean_irradiance, homogeneity = process_radiometry_indirect('results/radiometry/LST1-01G01-UV01-00/indirect_irradiation.csv')
        print(f'Mean Irradiance: {mean_irradiance:.2f} W/m²')
        print(f'Homogeneity Index: {homogeneity:.4f}')

        # Visualize the irradiance distribution
        from utils.process_radiometry_data import plot_irradiance_radiometry
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 6))
        plot_irradiance_radiometry(ax, irr, title='Indirect Irradiation - Irradiance Distribution')
        plt.show()