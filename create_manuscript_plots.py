from visualization.irradiance_plots_SI import create_irradiance_plots_SI
from visualization.plots_validation import create_plots_validation
from visualization.bar_plot_final_comparison import create_bar_plot_final_comparison
from visualization.bar_plot_screened_cands import create_bar_plot_screened_cands
from visualization.ray_independence_plots import create_ray_independence_plots_SI
from visualization.cand_plots import plot_cand
from visualization.pareto_plots_H_vs_E_ratio import create_pareto_8_leds_13cm, create_pareto_1_16_leds_5_15cm

if __name__ == "__main__":
    # plot_cand("best_manual") # uncomment to generate plot of best manual configuration
    # plot_cand("GRN3") # uncomment to generate plot of configuration GRN3
    # create_plots_validation() # uncomment to generate validation plots for best manual configuration
    # create_pareto_1_16_leds_5_15cm() # uncomment to generate pareto plot for 1-16 LEDs and 5-15 cm. For the final plot, the subplots were arranged manually in Inkscape.
    # create_pareto_8_leds_13cm() # uncomment to generate pareto plot for 8 LEDs and 13 cm. For the final plot, the subplots were arranged manually in Inkscape.
    # create_bar_plot_screened_cands() # uncomment to generate screened candidates bar plot
    create_bar_plot_final_comparison() # uncomment to generate final comparison bar plot
    # create_ray_independence_plots_SI() # uncomment to generate SI ray independence plots
    # create_irradiance_plots_SI() # uncomment to generate SI irradiance plots
