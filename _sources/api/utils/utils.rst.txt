LED 
====

.. automodule:: utils.led
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
---------------

.. code-block:: python

    from utils.led import LED

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

    from utils.grid import Grid

    grid = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.0)
    print(grid.grid_points)
