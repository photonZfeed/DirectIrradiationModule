"""
Evaluation of measurements and simulations for the Direct Irradiation Module project.

This package processes the radiometry and ray-tracing data and contains the analyses behind the
numbers and tables of the manuscript and its Supporting Information:

- :mod:`analysis.process_radiometry_data`: Processing and visualization of direct-irradiation
  radiometry measurements.
- :mod:`analysis.process_radiometry_indirect`: Processing of indirect-irradiation radiometry
  measurements.
- :mod:`analysis.process_raytracing_data`: Parsing and visualization of Speos ray-tracing output files.
- :mod:`analysis.compare_methods`: Point-by-point comparison of geometric model, ray tracing and
  radiometry (Supporting Information S6).
- :mod:`analysis.fast_evaluation`: Exact evaluation of millions of LED configurations with the
  geometric model and the configuration selection of the manuscript.
- :mod:`analysis.source_model`: Point-source assumption of the geometric model (Supporting
  Information S9, Table S3).
- :mod:`analysis.uncertainty`: Metrics of all methods, reflector effects and the uncertainties of
  ray tracing and radiometry (Supporting Information S10, Table S4).
- :mod:`analysis.irradiance_bound`: Influence of the irradiance bound (Supporting Information S11,
  Table S5).
- :mod:`analysis.mc_sampling`: Monte-Carlo samples of the convergence study and their evaluation.
- :mod:`analysis.mc_convergence`: Convergence of the Monte-Carlo optimization (Supporting
  Information S12).

All functions use paths relative to the repository root, like the rest of the repository; run them
from there.
"""
