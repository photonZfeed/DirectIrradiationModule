"""Figures for Reviewer 3, comment 3 (SI Section S11), in the style of the manuscript figures.

Reads the sampled configurations and metrics of ``r3_mc_runs.py`` (``--results``) and the tables of
``r3_3_convergence.py`` (``--tables``) and writes into ``--out``:

* ``r3_3_convergence.svg/.png`` (Figure S18): (a, b) non-dominated configurations of the first
  5 million candidates of every seed; (c, d) highest H * E_mean / E_bound versus the number of
  candidates (linear axis).
* ``r3_3_best_configurations.svg/.png`` (Figure S19): the configuration with the highest
  H * E_mean / E_bound of every seed at 5 million candidates, for both LED types.

The style is that of the manuscript figures, which the repository draws with the group's ICIW-Plots
package: styles ``ICIWstyle`` and the repository's ``visualization/publication_style.mplstyle``
(10 pt labels, no grid), square axes from ``ICIW_Plots.make_square_subplots``, a 19 cm wide figure,
sub-captions "a) 365 nm LEDs" below the panels, a framed legend on top, dotted lines between the
markers and the red cross for the selected configuration as in Figures 4 and 5
(``visualization/pareto_plots_H_vs_E_ratio.py``, ``optimizer.Optimizer.plot_pareto_front``); the
configuration maps follow Figure 1a (``systematic_sampler.Sampler.plot_configuration``: blue
points, dashed lines through the 2.5 cm mounting grid, axes "$X$ / cm" and "$Y$ / cm").
Every seed keeps its colour and marker in both figures. Tool metadata is stripped from the files.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import ICIW_Plots  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from ICIW_Plots import cm2inch, make_square_subplots  # noqa: E402

from r2_common import DEFAULT_OUT, HALF, STEP, common_args, index_to_xy, pareto_mask_fast, repo_path  # noqa: E402

BOUND = 150.0
SIZE_FRONT = 5_000_000
WAVELENGTH = {"UV": "365 nm LEDs", "GRN": "530 nm LEDs"}
PUBLISHED = {"UV": "UV3", "GRN": "GRN3"}
PLANE = (33.0, 34.0)                      # detection plane width x height in cm (the module footprint)
FIG_WIDTH = 19 * cm2inch                  # as Figures 4 and 5 of the manuscript
MARKERS = ("o", "^", "v", "p", "h", "<", ">")
X_LABEL = r"$E_\mathrm{mean}\,/\,E_\mathrm{bound}$ / 1"            # Figures 4 and 5
Y_LABEL = "Homogeneity / 1"                                          # Figures 4 and 5
OBJ_LABEL = r"$H \cdot E_\mathrm{mean}\,/\,E_\mathrm{bound}$ / 1"
N_LABEL = r"Number of candidates / $10^6$"
_PNG_META = {"Software": None}
_SVG_META = {"Creator": None, "Date": None, "Format": None, "Type": None}


def use_style(repo: Path) -> None:
    """ICIWstyle plus the repository's publication style, applied by file path (works on every OS)."""
    pkg = Path(ICIW_Plots.__file__).parent
    plt.style.use([str(pkg / "ICIWstyle.mplstyle"), str(repo / "visualization" / "publication_style.mplstyle")])


def colors() -> list[str]:
    """ICIW colour cycle without red, which marks the selected configuration (as in Figures 4 and 5)."""
    c = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    return [c[0]] + c[1:3] + c[4:]


def style_of(seeds: list[int]) -> dict[int, tuple[str, str]]:
    cl = colors()
    return {s: (cl[i], MARKERS[i]) for i, s in enumerate(seeds)}


def red() -> str:
    return plt.rcParams["axes.prop_cycle"].by_key()["color"][3]


def save(fig, out: Path, name: str) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.svg", metadata=_SVG_META)
    png = out / f"{name}.png"
    fig.savefig(png, dpi=300, metadata=_PNG_META)
    plt.close(fig)
    return png


def config_file(res: Path, seed: int) -> Path:
    return res / f"mc_seed{seed}_10000000_configs.npz"


def front(res: Path, seed: int, key: str, n: int) -> tuple[np.ndarray, np.ndarray]:
    m = np.load(Path(str(config_file(res, seed)).replace("_configs.npz", f"_metrics_{key}.npz")))
    mean, std, mx = m["mean"][:n], m["std"][:n], m["max"][:n]
    feas = mx <= BOUND
    e, h = np.where(feas, mean / BOUND, -1.0), np.where(feas, 1 - std / mean, -1.0)
    nd = pareto_mask_fast(e, h) & feas
    o = np.argsort(e[nd])
    return e[nd][o], h[nd][o]


def mark_published(ax, x: float, y: float, name: str, dx: float, dy: float, ha: str) -> None:
    ax.scatter(x, y, s=64, color=red(), marker="X", zorder=5, label="UV3 / GRN3")
    ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, fontweight="bold")


def convergence_figure(res: Path, tables: Path, out: Path, seeds: list[int]) -> Path:
    conv = pd.read_csv(tables / "r3_3_convergence.csv")
    best = pd.read_csv(tables / "r3_3_best_configurations.csv")
    sty = style_of(seeds)
    fig = plt.figure(figsize=(FIG_WIDTH, 21.5 * cm2inch))
    axs = make_square_subplots(fig, ax_width=7 * cm2inch, ax_layout=(2, 2), h_sep=2.5 * cm2inch,
                               v_sep=3.3 * cm2inch, left_h=2.0 * cm2inch, bottom_v=2.4 * cm2inch,
                               xlabel=[[X_LABEL, X_LABEL], [N_LABEL, N_LABEL]],
                               ylabel=[[Y_LABEL, Y_LABEL], [OBJ_LABEL, OBJ_LABEL]])
    h_min = []
    for col, key in enumerate(("UV", "GRN")):
        ax = axs[0, col]
        for s in seeds:
            c, mk = sty[s]
            e, h = front(res, s, key, SIZE_FRONT)
            h_min.append(h.min())
            ax.plot(e, h, ":", marker=mk, color=c, label=f"Seed {s}")
        p = best[(best.led == key) & (best.seed == 42) & (best.n == SIZE_FRONT)].iloc[0]
        mark_published(ax, p.E_mean / BOUND, p.H, PUBLISHED[key], -6, -16, "right")
        ax.set_title(f"{'ab'[col]}) {WAVELENGTH[key]}", y=-0.33)

        ax2 = axs[1, col]
        for s in seeds:
            c, mk = sty[s]
            g = conv[(conv.led == key) & (conv.seed == s)].sort_values("n")
            ax2.plot(g.n / 1e6, g.best_score, ":", marker=mk, color=c, label=f"Seed {s}")
        g = conv[(conv.led == key) & (conv.seed == 42) & (conv.n == SIZE_FRONT)].iloc[0]
        mark_published(ax2, g.n / 1e6, g.best_score, PUBLISHED[key], 0, 7, "center")
        ax2.set_xlim(0, 10.5)
        ax2.set_xticks([0, 2.5, 5, 7.5, 10], labels=["0", "2.5", "5", "7.5", "10"])
        lo, hi = conv[conv.led == key].best_score.agg(["min", "max"])
        pad = 0.1 * (hi - lo)
        ax2.set_ylim(lo - pad, hi + 2.5 * pad)
        ax2.set_title(f"{'cd'[col]}) {WAVELENGTH[key]}", y=-0.33)
    lo = np.floor((min(h_min) - 0.02) * 10) / 10
    for ax in axs[0]:
        ax.set_ylim(lo, 1.0)                      # shared homogeneity range, as in Figure 5
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(handles), bbox_to_anchor=(0.5, 0.995))
    return save(fig, out, "r3_3_convergence")


def configuration_figure(tables: Path, out: Path, seeds: list[int]) -> Path:
    best = pd.read_csv(tables / "r3_3_best_configurations.csv")
    w, h = PLANE
    n = len(seeds)
    ax_w, h_sep, left = 3.5, 0.75, 2.4
    fig = plt.figure(figsize=(FIG_WIDTH, 15.0 * cm2inch))
    axs = make_square_subplots(fig, ax_width=ax_w * cm2inch, ax_layout=(2, n), h_sep=h_sep * cm2inch,
                               v_sep=3.7 * cm2inch, left_h=left * cm2inch, bottom_v=2.6 * cm2inch,
                               sharey=True, xlabel=r"$X$ / cm",
                               ylabel=[f"{WAVELENGTH['UV']}\n$Y$ / cm", f"{WAVELENGTH['GRN']}\n$Y$ / cm"])
    grid = np.arange(-HALF, HALF + 1) * STEP
    for row, key in enumerate(("UV", "GRN")):
        for col, s in enumerate(seeds):
            ax = axs[row, col]
            b = best[(best.led == key) & (best.seed == s) & (best.n == SIZE_FRONT)].iloc[0]
            x, y = index_to_xy(np.array([int(i) for i in b.position_indices.split()]))
            ax.scatter(x, y, c="blue", s=20, zorder=3)        # Figure 1a (Sampler.plot_configuration)
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
        x0 = left * cm2inch / FIG_WIDTH
        x1 = (left + n * ax_w + (n - 1) * h_sep) * cm2inch / FIG_WIDTH
        y = axs[row, 0].get_position().y0 - 1.75 * cm2inch / fig.get_figheight()
        fig.text((x0 + x1) / 2, y, f"{'ab'[row]}) {WAVELENGTH[key]}", ha="center", va="top",
                 fontsize=plt.rcParams["axes.titlesize"])
    return save(fig, out, "r3_3_best_configurations")


def main() -> None:
    ap = common_args(__doc__)
    ap.add_argument("--results", default=str(DEFAULT_OUT), help="folder with the sampled configurations and metrics")
    ap.add_argument("--tables", default=None, help="folder with r3_3_convergence.csv (default: --results)")
    ap.add_argument("--seeds", default="42,1,2,3")
    args = ap.parse_args()
    use_style(repo_path(args.repo))
    res, out = Path(args.results), Path(args.out)
    tables = Path(args.tables) if args.tables else res
    seeds = [int(s) for s in args.seeds.split(",")]
    for p in (convergence_figure(res, tables, out, seeds), configuration_figure(tables, out, seeds)):
        print(f"written {p} and {p.with_suffix('.svg')}")


if __name__ == "__main__":
    main()
