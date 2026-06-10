Systematic Sampler
========================

.. automodule:: systematic_sampler
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
-----------------

   .. code-block:: python

      from utils.grid import Grid
      from utils.led import LED
      from systematic_sampler import Sampler, save_configurations

      # initialize the grid and LED
      G = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
      led = LED(
         name="LST1-01G01-UV01-00",
         optical_power=0.875,  # @ 500 mA in W
         intensity_data_file="utils/angle_intensity_data/LST1-01G01-UV01-00.csv"
      )

      # Create a systematic sampler for 8 LEDs at height 13 cm
      sampler = Sampler(G=G, led_count=8, height=13)
      # Generate all configurations and save them to a file
      samples_sys = sampler.generate_all_configurations()
      save_configurations(samples_sys, sampler, "results/examples/example_systematic_8_leds_13_cm.npz")