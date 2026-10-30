"""
Figures S18 and S19: convergence of the Monte-Carlo optimization.

* :func:`create_mc_convergence_plot` (Figure S18): (a, b) non-dominated configurations of the first
  5 million candidates of every seed; (c, d) highest H E_mean / E_bound versus the number of
  candidates.
* :func:`create_mc_best_configurations_plot` (Figure S19): the configuration with the highest
  H E_mean / E_bound of every seed at 5 million candidates, for both LED types.

The figures are drawn from the tables of :func:`analysis.mc_convergence.convergence_tables` in the
style of Figures 1a, 4 and 5 of the manuscript: square axes, sub-captions below the panels, dotted
lines between the markers and a red cross for UV3 / GRN3.

**Usage Example**::

    from analysis.mc_convergence import convergence_tables
    from visualization.mc_convergence_plots import create_mc_convergence_plot, create_mc_best_configurations_plot
    best, fronts, _ = convergence_tables()
    create_mc_convergence_plot(best, fronts)
    create_mc_best_configurations_plot(best)
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ICIW_Plots import cm2inch, make_square_subplots

from analysis.fast_evaluation import E_BOUND, HALF, STEP
from analysis.mc_convergence import SIZE_FRONT
from visualization.style import finish, use_publication_style

WAVELENGTH = {"UV": "365 nm LEDs", "GRN": "530 nm LEDs"}
PUBLISHED = {"UV": "UV3", "GRN": "GRN3"}
#: Detection plane width x height in cm (the module footprint).
PLANE = (33.0, 34.0)
#: Figure width in cm, as Figures 4 and 5 of the manuscript.
FIG_WIDTH_CM = 19.0
MARKERS = ("o", "^", "v", "p", "h", "<", ">")
X_LABEL = r"$E_\mathrm{mean}\,/\,E_\mathrm{bound}$ / 1"
Y_LABEL = "Homogeneity / 1"
OBJ_LABEL = r"$H \cdot E_\mathrm{mean}\,/\,E_\mathrm{bound}$ / 1"
N_LABEL = r"Number of candidates / $10^6$"
#: Position of the sub-captions below the axes (axes fraction).
SUBCAPTION_Y = -0.255


def _colors() -> list[str]:
    """Colour cycle without red, which marks the selected configuration (as in Figures 4 and 5)."""
    c = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    return [c[0]] + c[1:3] + c[4:]


def _red() -> str:
    return plt.rcParams["axes.prop_cycle"].by_key()["color"][3]


def _mark_published(ax, x: float, y: float, name: str, dx: float, dy: float, ha: str) -> None:
    ax.scatter(x, y, s=64, color=_red(), marker="X", zorder=5, label="UV3 / GRN3")
    ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, fontweight="bold")


def create_mc_convergence_plot(best: pd.DataFrame, fronts: pd.DataFrame, figure_path: str = "figures/",
                               show: bool = True) -> str:
    """
    Draw Figure S18 and save it as ``mc_convergence.svg``.

    :param best: Best configurations per LED type, seed and sample size
        (:func:`analysis.mc_convergence.convergence_tables`).
    :type best: pandas.DataFrame
    :param fronts: Fronts at 5 million candidates in long format (same function).
    :type fronts: pandas.DataFrame
    :param figure_path: Directory the figure is written to.
    :type figure_path: str, optional
    :param show: Display the figure instead of closing it after saving.
    :type show: bool, optional
    :returns: Path of the figure.
    :rtype: str
    """
    use_publication_style()
    seeds = list(dict.fromkeys(best.seed))
    style = {s: (c, m) for s, c, m in zip(seeds, _colors(), MARKERS)}
    # vertical layout in cm, compact so that the figure and its caption fit on one page
    ax_w, v_sep, bottom, top = 7.0, 1.80, 1.60, 1.11
    fig = plt.figure(figsize=(FIG_WIDTH_CM * cm2inch, (bottom + 2 * ax_w + v_sep + top) * cm2inch))
    axs = make_square_subplots(fig, ax_width=ax_w * cm2inch, ax_layout=(2, 2), h_sep=2.5 * cm2inch,
                               v_sep=v_sep * cm2inch, left_h=2.0 * cm2inch, bottom_v=bottom * cm2inch,
                               xlabel=[[X_LABEL, X_LABEL], [N_LABEL, N_LABEL]],
                               ylabel=[[Y_LABEL, Y_LABEL], [OBJ_LABEL, OBJ_LABEL]])
    h_min = []
    for col, key in enumerate(("UV", "GRN")):
        ax = axs[0, col]
        for s in seeds:
            c, mk = style[s]
            f = fronts[(fronts.led == key) & (fronts.seed == s)]
            h_min.append(f.H.min())
            ax.plot(f.E_ratio, f.H, ":", marker=mk, color=c, label=f"Seed {s}")
        p = best[(best.led == key) & (best.seed == 42) & (best.n == SIZE_FRONT)].iloc[0]
        _mark_published(ax, p.E_mean / E_BOUND, p.H, PUBLISHED[key], -6, -16, "right")
        ax.set_title(f"{'ab'[col]}) {WAVELENGTH[key]}", y=SUBCAPTION_Y)

        ax2 = axs[1, col]
        for s in seeds:
            c, mk = style[s]
            g = best[(best.led == key) & (best.seed == s)].sort_values("n")
            ax2.plot(g.n / 1e6, g.score, ":", marker=mk, color=c, label=f"Seed {s}")
        _mark_published(ax2, SIZE_FRONT / 1e6, p.score, PUBLISHED[key], 0, 7, "center")
        ax2.set_xlim(0, 10.5)
        ax2.set_xticks([0, 2.5, 5, 7.5, 10], labels=["0", "2.5", "5", "7.5", "10"])
        lo, hi = best[best.led == key].score.agg(["min", "max"])
        pad = 0.1 * (hi - lo)
        ax2.set_ylim(lo - pad, hi + 2.5 * pad)
        ax2.set_title(f"{'cd'[col]}) {WAVELENGTH[key]}", y=SUBCAPTION_Y)
    lo = np.floor((min(h_min) - 0.02) * 10) / 10
    for ax in axs[0]:
        ax.set_ylim(lo, 1.0)                      # shared homogeneity range, as in Figure 5
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(handles), bbox_to_anchor=(0.5, 1.0))
    os.makedirs(figure_path, exist_ok=True)
    path = os.path.join(figure_path, "mc_convergence.svg")
    fig.savefig(path)
    finish(fig, show)
    return path


def create_mc_best_configurations_plot(best: pd.DataFrame, figure_path: str = "figures/", show: bool = True) -> str:
    """
    Draw Figure S19 and save it as ``mc_best_configurations.svg``.

    :param best: Best configurations per LED type, seed and sample size
        (:func:`analysis.mc_convergence.convergence_tables`).
    :type best: pandas.DataFrame
    :param figure_path: Directory the figure is written to.
    :type figure_path: str, optional
    :param show: Display the figure instead of closing it after saving.
    :type show: bool, optional
    :returns: Path of the figure.
    :rtype: str
    """
    use_publication_style()
    seeds = list(dict.fromkeys(best.seed))
    w, h = PLANE
    n = len(seeds)
    ax_w, h_sep, left = 3.5, 0.75, 2.05
    fig = plt.figure(figsize=(FIG_WIDTH_CM * cm2inch, 15.0 * cm2inch))
    axs = make_square_subplots(fig, ax_width=ax_w * cm2inch, ax_layout=(2, n), h_sep=h_sep * cm2inch,
                               v_sep=3.7 * cm2inch, left_h=left * cm2inch, bottom_v=2.6 * cm2inch,
                               sharey=True, xlabel=r"$X$ / cm", ylabel=r"$Y$ / cm")
    grid = np.arange(-HALF, HALF + 1) * STEP
    for row, key in enumerate(("UV", "GRN")):
        for col, s in enumerate(seeds):
            ax = axs[row, col]
            b = best[(best.led == key) & (best.seed == s) & (best.n == SIZE_FRONT)].iloc[0]
            xy = np.array([[float(v) for v in p.split()] for p in b.positions_cm.split(";")])
            ax.scatter(xy[:, 0], xy[:, 1], c="blue", s=20, zorder=3)     # as Figure 1a
            ax.set_xlim(-w / 2, w / 2)
            ax.set_ylim(-h / 2, h / 2)
            ax.set_xticks([-10, 0, 10])
            ax.set_yticks([-10, 0, 10])
            ax.set_xticks(grid, minor=True)
            ax.set_yticks(grid, minor=True)
            ax.grid(True, which="both", linestyle="--", linewidth=0.5)
            ax.set_axisbelow(True)
            ax.set_aspect("equal")
            if col > 0:
                ax.tick_params(labelleft=False)
            name = f" ({b.is_published})" if isinstance(b.is_published, str) and b.is_published else ""
            ax.text(0.5, 1.17, f"Seed {s}{name}", transform=ax.transAxes, ha="center", va="bottom",
                    fontweight="bold")
            ax.text(0.5, 1.03, f"{int(b.N)} LEDs, {int(b.z_cm)} cm", transform=ax.transAxes, ha="center",
                    va="bottom")
        # sub-caption centred under the row, as in Figures 4 and 5
        x0 = left / FIG_WIDTH_CM
        x1 = (left + n * ax_w + (n - 1) * h_sep) / FIG_WIDTH_CM
        y = axs[row, 0].get_position().y0 - 1.75 * cm2inch / fig.get_figheight()
        fig.text((x0 + x1) / 2, y, f"{'ab'[row]}) {WAVELENGTH[key]}", ha="center", va="top",
                 fontsize=plt.rcParams["axes.titlesize"])
    os.makedirs(figure_path, exist_ok=True)
    path = os.path.join(figure_path, "mc_best_configurations.svg")
    fig.savefig(path)
    finish(fig, show)
    return path
