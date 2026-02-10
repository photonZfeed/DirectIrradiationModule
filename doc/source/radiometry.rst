Radiometric Measurements
========================

Overview
--------

Radiometric measurements were conducted to experimentally characterize the
spatial irradiance distribution direct irradiation modules equipped with different
LED configurations.
All measurements were performed using the two-dimensional (2D) radiometry
procedure introduced by Sender et al. [1]_.

The measurements were used to validate geometric radiation field modelling and
ray tracing simulations and to assess the influence of LED configuration, led type, irradiation distance, 
and the use of reflectors.

Investigated LED Types
----------------------


Two commercially available high-power LEDs were investigated:

**UV LED**
    - *Peak emission wavelength*: 365 nm
    - *Model*: LST1-01G01-UV01-00
    - *Datasheet*: `Luminus Horticulture Starboards Industry Leading High Powered LED Starboards <https://www.mouser.de/datasheet/3/1525/1/NewEnergy_StarBoard_Horticulture_Luminus_DataSheet.pdf>`_

**Green LED**
    - *Peak emission wavelength*: 530 nm
    - *Model*: LST1-01F06-GRN1-00
    - *Datasheet*: `OSRAM Horticulture Starboards Industry Leading High Powered LED Starboards <https://www.mouser.de/datasheet/3/1525/1/NewEnergy_StarBoard_Horticulture_Osram_DataSheet.pdf>`_
  
Operating Conditions
--------------------

Radiometric measurements were carried out under the following operating conditions:

* UV LEDs were operated at a constant current of **500 mA** with an integration
  time of **2.88 ms**.
* Green LEDs were operated at a constant current of **350 mA** with an
  integration time of **15.78 ms**.

These operating points were chosen to maximize signal quality while avoiding
detector saturation.

LED Power Supply
-----------------

Radiometric scans were performed for multiple LED configurations comprising a
total of **8** or **16 LEDs**, both **with and without diffuse reflectors**.

For each wavelength, two LED rows were assembled, each consisting of
``2 × 8`` LEDs connected in series. The LED rows were operated independently
using two programmable DC power supplies (Korad Lab KA3005P) operated in
constant-current mode.

Measurement Distances
---------------------

Measurements were conducted at three defined vertical distances between the LED
mounting plane and the cosine corrector of the spectrometer.

Based on geometric radiation field modelling, nominal optimal distances of

* 6 cm,
* 8 cm, and
* 9 cm

were identified. To avoid mechanical collision between the irradiation module
and the spectrometer, an additional safety offset of **1 cm** was applied during
all measurements. Consequently, the effective measurement distances were:

* **7 cm**,
* **9 cm**, and
* **10 cm**.

Data Processing and Spectral Trimming
-------------------------------------

For data evaluation, the recorded spectral data were trimmed to the relevant
emission ranges of the respective LEDs prior to integration:

* **UV LEDs:** 300 nm – 450 nm  
* **Green LEDs:** 450 nm – 600 nm


Photon Flux Correction
----------------------

Using datasheet values for current-dependent photon flux [2]_, [3]_, the photon fluxes obtained from 2D radiometry were corrected
based on reference measurements acquired at an effective distance of **7 cm**
for direct irradiation **without reflectors**.

The resulting correction factors were subsequently applied to all measurements
recorded at the corresponding irradiation distances, including configurations
both with and without reflectors.

This procedure ensured consistent comparability of radiometric data across all
evaluated module configurations.

References
----------

.. [1] M. Sender, B. Wriedt, D. Ziegenbalg, React. Chem. Eng. 2021, 6, 1601–1613.
.. [2] New-EnergyLLC, “Luminus Horticulture Starboards Industry Leading High Powered LED Starboards,” https://www.mouser.de/datasheet/3/1525/1/NewEnergy_StarBoard_Horticulture_Luminus_DataSheet.pdf, 2025.
.. [3] New-EnergyLLC, “OSRAM Horticulture Starboards Industry Leading High Powered LED Starboards,” https://www.mouser.de/datasheet/3/1525/1/NewEnergy_StarBoard_Horticulture_Osram_DataSheet.pdf, 2025.