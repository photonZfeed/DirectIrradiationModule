import os

import matplotlib.pyplot as plt
import numpy as np
from ICIW_Plots import cm2inch, make_square_subplots

from analysis.compare_methods import CASES, LEDS, common_grid, compare_all, compare_case
from utils.led import LST1_01F06_GRN1_00, LST1_01G01_UV01_00
from visualization.style import finish, use_publication_style

#: Column titles of the three pairwise comparisons.
COLUMN_TITLES = ["geometric model\n$-$ ray tracing",
                 "geometric model\n$-$ radiometry",
                 "ray tracing\n$-$ radiometry"]

#: Legend entries of the three pairwise comparisons, together with the transparency used for them
#: in the bar charts. The same shading scheme as in the other bar charts of the manuscript is used.
COMPARISON_LEGEND = [("model $-$ ray tracing", 1.0),
                     ("model $-$ radiometry", 0.67),
                     ("ray tracing $-$ radiometry", 0.33)]


def _plot_residual(ax, delta: np.ndarray, x: np.ndarray, y: np.ndarray, limit: float):
    """
    Draw a single residual map with a diverging colour scale centred at zero.

    :param ax: Axes to draw on.
    :type ax: :class:`matplotlib.axes.Axes`
    :param delta: Residual field in W m\\ :sup:`-2`.
    :type delta: numpy.ndarray
    :param x: x coordinates of the common grid in centimeters.
    :type x: numpy.ndarray
    :param y: y coordinates of the common grid in centimeters.
    :type y: numpy.ndarray
    :param limit: Symmetric limit of the colour scale in W m\\ :sup:`-2`.
    :type limit: float

    :returns: The image created by :meth:`matplotlib.axes.Axes.imshow`.
    :rtype: :class:`matplotlib.image.AxesImage`
    """
    return ax.imshow(delta, origin="lower", cmap="RdBu_r", vmin=-limit, vmax=limit,
                     extent=(x.min(), x.max(), y.min(), y.max()))


def create_residual_statistics_plot(statistics, figure_path: str = "figures/", show: bool = True):
    r"""
    Plot the root mean square error of the pairwise residuals as a grouped bar chart.

    Mean irradiance and homogeneity, which are compared in the manuscript, characterize each
    irradiance field on its own and are insensitive to deviations that cancel in the spatial
    average: two fields can share both values and still disagree at every position. The RMSE of the
    residual measures this local agreement between two methods and is therefore reported here in
    addition. The mean offset between two methods is not shown separately, because it is identical
    to the difference of the mean irradiances already reported in the manuscript.

    :param statistics: Table of residual statistics as returned by
        :func:`analysis.compare_methods.compare_all`.
    :type statistics: :class:`pandas.DataFrame`
    :param figure_path: Directory the figure is written to.
    :type figure_path: str, optional
    :param show: Display the figure(s) after saving; with False they are closed instead. Default is True.
    :type show: bool, optional

    :returns: None

    :side effects:
        - Saves ``difference_statistics.svg`` to ``figure_path``.
        - Displays the figure with :func:`matplotlib.pyplot.show`.
        - Prints the largest absolute residual of each LED type, which is quoted in the text of the
          supporting information.

    **Usage Example**::

        from analysis.compare_methods import compare_all
        from visualization.difference_plots import create_residual_statistics_plot
        _, statistics = compare_all()
        create_residual_statistics_plot(statistics)
    """
    use_publication_style()
    os.makedirs(figure_path, exist_ok=True)

    comparisons = [f"{first} - {second}" for first, second in
                   [("geometric model", "ray tracing"), ("geometric model", "radiometry"),
                    ("ray tracing", "radiometry")]]
    led_groups = [(LST1_01G01_UV01_00, "365 nm LEDs"), (LST1_01F06_GRN1_00, "530 nm LEDs")]
    # the bars are drawn edge to edge, so their width also sets the distance between the value
    # labels above them
    width = 0.29

    # square axes, as in the other multi-panel figures of the manuscript
    fig = plt.figure(figsize=(19 * cm2inch, 10.5 * cm2inch))
    axs = make_square_subplots(fig=fig, ax_width=8 * cm2inch, ax_layout=(1, len(led_groups)),
                               h_sep=1.5 * cm2inch, v_sep=0 * cm2inch, sharex=False, sharey=True,
                               sharelabel=False, xlabel="", ylabel="RMSE / %")
    for index, (led, led_title) in enumerate(led_groups):
        ax = axs[0, index]
        subset = statistics[statistics["LED"] == led.name]
        configurations = list(dict.fromkeys(subset["configuration"]))
        positions = range(len(configurations))
        for offset, (comparison, (legend_label, alpha)) in zip(
                (-width, 0.0, width), zip(comparisons, COMPARISON_LEGEND)):
            values = [float(subset[(subset["configuration"] == configuration)
                                   & (subset["comparison"] == comparison)]["RMSE / %"].iloc[0])
                      for configuration in configurations]
            ax.bar([position + offset for position in positions], values, width=width,
                   label=legend_label if index == 0 else None, color="C0", alpha=alpha)
        ax.set_xticks(list(positions))
        ax.set_xticklabels([configuration.replace(" ", "\n") for configuration in configurations])
        ax.set_title(f"{'ab'[index]}) {led_title}", y=-0.3)
        for rect in ax.patches:
            height = rect.get_height()
            # rotated because three values per configuration do not fit next to each other
            ax.annotate(f"{height:.1f}", xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points", ha="center", va="bottom",
                        fontsize=8, rotation=90)
    axs[0, 0].set_ylim(0, 1.35 * max(rect.get_height() for ax in axs.flat for rect in ax.patches))

    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0))
    fig.savefig(os.path.join(figure_path, "difference_statistics.svg"), dpi=300,
                bbox_inches="tight")
    finish(fig, show)

    for led, led_title in led_groups:
        largest = statistics[statistics["LED"] == led.name]["max |dE| / %"].max()
        print(f"largest absolute residual, {led_title}: {largest:.1f} % of the mean irradiance")


def create_difference_plots_SI(results_path: str = "results/", figure_path: str = "figures/", show: bool = True,
                               table_path: str | None = None):
    r"""
    Create the residual maps and the residual statistics of all compared configurations.

    One composite figure per LED type is generated, with one row per configuration (best manual
    arrangement and the three near-Pareto-optimal configurations) and one column per pairwise
    comparison between geometric model, ray tracing and radiometry. Because the irradiance level
    differs between the configurations, each row carries its own colour bar. In addition, the
    statistics of all residual fields are written as a CSV file.

    :param results_path: Path of the results directory.
    :type results_path: str, optional
    :param figure_path: Directory the figures are written to.
    :type figure_path: str, optional
    :param show: Display the figure(s) after saving; with False they are closed instead. Default is True.
    :type show: bool, optional
    :param table_path: Directory ``difference_summary.csv`` is written to. Default None: ``<results_path>/difference_maps/``.
    :type table_path: str | None, optional

    :returns: None

    :raises FileNotFoundError: If a required radiometry, ray tracing or configuration file is
        missing.

    :side effects:
        - Saves ``difference_maps_<365nm|530nm>.svg`` to ``figure_path``.
        - Saves ``difference_summary.csv`` to ``table_path`` (default ``<results_path>/difference_maps/``).
        - Displays the figures with :func:`matplotlib.pyplot.show`.

    .. note::
        The ray tracing fields are block-averaged onto the common evaluation grid instead of being
        interpolated. Without this averaging, the Monte-Carlo noise of the 1 mm mesh dominates the
        residual maps.

    **Usage Example**::

        from visualization.difference_plots import create_difference_plots_SI
        create_difference_plots_SI()
    """
    use_publication_style()

    x, y = common_grid()
    residuals_all, statistics = compare_all(results_path)

    export_path = table_path if table_path is not None else os.path.join(results_path, "difference_maps")
    os.makedirs(export_path, exist_ok=True)
    os.makedirs(figure_path, exist_ok=True)
    statistics.to_csv(os.path.join(export_path, "difference_summary.csv"), index=False)

    for led, tag in [(LST1_01G01_UV01_00, "365nm"), (LST1_01F06_GRN1_00, "530nm")]:
        cases = CASES[led.name]
        # one common colour scale for all configurations of an LED type, so that the panels can be
        # compared with each other directly
        limit = float(np.nanmax([np.nanmax(np.abs(delta))
                                 for case in cases
                                 for delta in residuals_all[led.name][case.label].values()]))
        fig, axs = plt.subplots(len(cases), 3, figsize=(19 * cm2inch, 5.5 * len(cases) * cm2inch),
                                constrained_layout=True)
        for row, case in enumerate(cases):
            residuals = residuals_all[led.name][case.label]
            for col, (label, delta) in enumerate(residuals.items()):
                ax = axs[row, col]
                image = _plot_residual(ax, delta, x, y, limit)
                if row == 0:
                    ax.set_title(COLUMN_TITLES[col])
                if col == 0:
                    ax.set_ylabel(f"{case.label}\n$y$ / cm")
                else:
                    ax.set_yticklabels([])
                if row == len(cases) - 1:
                    ax.set_xlabel("$x$ / cm")
                else:
                    ax.set_xticklabels([])
        fig.colorbar(image, ax=axs, shrink=0.6, label=r"$\Delta E$ / W m$^{-2}$")

        fig.savefig(os.path.join(figure_path, f"difference_maps_{tag}.svg"), dpi=300,
                    bbox_inches="tight")
        finish(fig, show)

    create_residual_statistics_plot(statistics, figure_path, show)

    with pd_option_context():
        print(statistics.round(2).to_string(index=False))


def pd_option_context():
    """
    Return a pandas option context with a wide console output.

    :returns: Context manager that widens the console output for the summary tables.
    :rtype: pandas.option_context
    """
    import pandas as pd
    return pd.option_context("display.width", 220, "display.max_columns", 20)
