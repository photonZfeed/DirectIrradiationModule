"""
Visualization utilities for the Direct Irradiation Module project.

This package contains publication-quality plotting scripts used to generate
the figures presented in the associated manuscript and its Supporting Information:

- :mod:`visualization.style`: Plot style shared by all figures.
- :mod:`visualization.cand_plots`: LED configurations (Figures 1a and 6a).
- :mod:`visualization.plots_validation`: Validation of the geometric model (Figure 3).
- :mod:`visualization.pareto_plots_H_vs_E_ratio`: Pareto fronts (Figures 4 and 5).
- :mod:`visualization.bar_plot_screened_cands`: Configurations UV1-UV3 and GRN1-GRN3 with and
  without reflectors (Figures 7 and S17).
- :mod:`visualization.bar_plot_final_comparison`: Comparison with the indirect irradiation
  module (Figure 8).
- :mod:`visualization.ray_independence_plots`: Ray-count independence plots
  assessing the statistical convergence of Speos raytracing simulations (Figure S5).
- :mod:`visualization.irradiance_plots_SI`: Irradiance comparison figures
  (geometric model vs. raytracing vs. radiometry) in SI units (Figures S6-S13).
- :mod:`visualization.difference_plots`: Residual maps and statistics (Figures S14-S16).
- :mod:`visualization.mc_convergence_plots`: Convergence of the Monte-Carlo optimization
  (Figures S18 and S19).

``create_manuscript_plots_tables.py`` in the repository root generates all of them.
"""
