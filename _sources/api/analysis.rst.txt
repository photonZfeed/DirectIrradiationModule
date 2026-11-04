Analysis
========

The ``analysis`` package evaluates the radiometry and ray-tracing data and contains the analyses
behind the numbers and tables of the publication. :doc:`../reproduce` lists which function
produces which figure, table or number.

Raytracing Data Processing
--------------------------

.. automodule:: analysis.process_raytracing_data
    :members:
    :undoc-members:
    :show-inheritance:

Example Usage
^^^^^^^^^^^^^

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

Radiometry Data Processing: Direct Irradiation Module
-----------------------------------------------------

.. automodule:: analysis.process_radiometry_data
    :members:
    :undoc-members:
    :show-inheritance:

Example Usage
^^^^^^^^^^^^^

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

Radiometry Data Processing: Indirect Irradiation Module
-------------------------------------------------------

.. automodule:: analysis.process_radiometry_indirect
    :members:
    :undoc-members:
    :show-inheritance:

Example Usage
^^^^^^^^^^^^^

.. code-block:: python

    # calculate mean irradiance and homogeneity index from an indirect irradiation measurement
    irr, mean_irradiance, homogeneity = process_radiometry_data('results/radiometry/LST1-01G01-UV01-00/indirect_irradiation.csv')
    print(f'Mean Irradiance: {mean_irradiance:.2f} W/m²')
    print(f'Homogeneity Index: {homogeneity:.4f}')

    # Visualize the irradiance distribution using the raytracing plot helper
    # (irr is a DataFrame with columns ['x', 'y', 'z'])
    from analysis.process_raytracing_data import plot_irradiance_raytracing
    from matplotlib import pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))
    plot_irradiance_raytracing(ax, irr.rename(columns={'z': 'Irradiance'}), title='Indirect Irradiation - Irradiance Distribution')
    plt.show()

Method Comparison
-----------------

.. automodule:: analysis.compare_methods
    :members:
    :undoc-members:
    :show-inheritance:

Example Usage
^^^^^^^^^^^^^

.. code-block:: python

    from analysis.compare_methods import CASES, common_grid, compare_all, compare_case
    from utils.led import LST1_01G01_UV01_00

    # compare all configurations of both LED types on the common evaluation grid
    residuals, statistics = compare_all('results/')
    print(statistics[['LED', 'configuration', 'comparison', 'RMSE / %']].round(1))

    # evaluate a single configuration
    x, y = common_grid()
    case = CASES[LST1_01G01_UV01_00.name][0]  # the best manual arrangement
    residual_fields, rows = compare_case(LST1_01G01_UV01_00, case, x, y)
    delta = residual_fields['geometric model - ray tracing']
    print(f"RMSE: {rows[0]['RMSE / %']:.1f} % of the mean irradiance")

Fast Evaluation of Many Configurations
--------------------------------------

.. automodule:: analysis.fast_evaluation
    :members:
    :undoc-members:
    :show-inheritance:

Example Usage
^^^^^^^^^^^^^

.. code-block:: python

    import numpy as np
    from analysis.fast_evaluation import FieldBank, LEDS, configs_from_object_array

    d = np.load('results/sampled_configs/systematic_8_leds_13_cm_482931_samples.npz', allow_pickle=True)
    pos = configs_from_object_array(d['configurations'])
    z = np.full(len(pos), 13)
    mean, std, max_irradiance = FieldBank(LEDS['UV'], heights=(13,)).evaluate(pos, z)

Point-Source Assumption
-----------------------

.. automodule:: analysis.source_model
    :members:
    :undoc-members:
    :show-inheritance:

Metrics, Reflector Effects and Uncertainties
--------------------------------------------

.. automodule:: analysis.uncertainty
    :members:
    :undoc-members:
    :show-inheritance:

Irradiance Bound
----------------

.. automodule:: analysis.irradiance_bound
    :members:
    :undoc-members:
    :show-inheritance:

Monte-Carlo Samples of the Convergence Study
--------------------------------------------

.. automodule:: analysis.mc_sampling
    :members:
    :undoc-members:
    :show-inheritance:

Convergence of the Monte-Carlo Optimization
-------------------------------------------

.. automodule:: analysis.mc_convergence
    :members:
    :undoc-members:
    :show-inheritance:
