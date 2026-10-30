Radiometric Measurements
========================

Overview
--------

Radiometric measurements were conducted to experimentally characterize the
spatial irradiance distribution of direct irradiation modules equipped with
different LED configurations. All measurements were performed using the
two-dimensional (2D) radiometry procedure introduced by Sender et al. [1]_.

The measurements were used to validate the geometric radiation field model and
the ray-tracing simulations and to assess the influence of the LED configuration,
the LED type, the irradiation distance and the use of reflectors.

Investigated LED Types
----------------------

Two commercially available high-power LEDs on star boards were investigated:

**UV LED**
    - *Peak emission wavelength*: 365 nm
    - *Model*: LST1-01G01-UV01-00 (emitter: Luminus SST-10-UV)
    - *Datasheets*: star board [2]_, emitter [4]_

**Green LED**
    - *Peak emission wavelength*: 530 nm
    - *Model*: LST1-01F06-GRN1-00 (emitter: OSRAM OSLON SSL 120, GT CSSPM1.13)
    - *Datasheets*: star board [3]_, emitter [5]_

The relative angular intensity distributions used by the geometric model
(``utils/angle_intensity_data/``) were taken from the emitter datasheets.

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
----------------

Radiometric scans were performed for multiple LED configurations comprising a
total of **8** or **16 LEDs**, both **with and without diffuse reflectors**.

For each wavelength, two LED rows were assembled, each consisting of
``2 × 8`` LEDs connected in series. The LED rows were operated independently
using two programmable DC power supplies (Korad Lab KA3005P) operated in
constant-current mode.

Measurement Distances
---------------------

The distances were measured vertically between the LED mounting plane and the
cosine corrector of the spectrometer. For the configurations UV1–UV3 and
GRN1–GRN3, the geometric radiation field model gave the optimal distances of

* 6 cm,
* 8 cm, and
* 9 cm.

To avoid a mechanical collision between the irradiation module and the
spectrometer, an additional safety offset of **1 cm** was applied during all
measurements. Consequently, the effective measurement distances were:

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

The 2D radiometry is not calibrated in absolute units. Using the datasheet values
for the current-dependent photon flux [2]_ [3]_, the photon fluxes obtained from
2D radiometry were therefore corrected with a reference measurement: the reference
scan is scaled such that its integral equals the radiant flux emitted by all LEDs
of the module, and the resulting correction factor is applied to all measurements
of the same group, both without and with reflectors.

* For the **best manual configuration**, its own scan **without reflectors** served
  as the reference.
* For the **optimized configurations** UV1–UV3 and GRN1–GRN3, the scans of **UV3
  and GRN3 with reflectors** at 7 cm served as the references.

This procedure ensures that the radiometric data of all evaluated module
configurations can be compared with each other. Because it assumes that the
entire emitted flux reaches the scanned area, it overestimates the absolute mean
irradiance; Supporting Information Section S10 and
:func:`analysis.uncertainty.radiometry_corrections` estimate this effect and the
effect of the 1 cm safety offset with the geometric model. The evaluation is
implemented in :mod:`analysis.process_radiometry_data`; the processed scans are
stored in ``results/radiometry/``.

References
----------

.. [1] M. Sender, B. Wriedt, D. Ziegenbalg, Radiometric measurement techniques for in-depth characterization of photoreactors – part 1: 2 dimensional radiometry, *React. Chem. Eng.* **2021**, 6, 1601–1613, DOI: `10.1039/D0RE00456A <https://doi.org/10.1039/D0RE00456A>`_.
.. [2] New-Energy LLC, “Luminus Horticulture Starboards Industry Leading High Powered LED Starboards,” https://www.mouser.de/datasheet/3/1525/1/NewEnergy_StarBoard_Horticulture_Luminus_DataSheet.pdf, 2025.
.. [3] New-Energy LLC, “OSRAM Horticulture Starboards Industry Leading High Powered LED Starboards,” https://www.mouser.de/datasheet/3/1525/1/NewEnergy_StarBoard_Horticulture_Osram_DataSheet.pdf, 2025.
.. [4] Luminus, “SST-10-UV Product Datasheet,” https://led-tech-shop.s3.eu-central-1.amazonaws.com/files/208/Luminus_SST-10-UV_Datasheet.pdf, 2023.
.. [5] OSRAM, “OSLON SSL 120, GT CSSPM1.13,” https://look.ams-osram.com/m/3d194621964549f/original/GT-CSSPM1-13.pdf, 2020.
