"""Figure for Reviewer 3, comment 3: Pareto fronts per seed and convergence with sample size.

Reads the runs of ``r3_mc_runs.py`` and ``r3_3_convergence.csv`` and writes
``<out>/r3_3_convergence.png`` (2 x 2 panels: fronts at 5 million candidates per seed for
365 nm and 530 nm; best product score H * E_mean/E_bound versus sample size).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from r2_common import DEFAULT_OUT, common_args, pareto_mask_fast

BOUND = 150.0
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]   # fixed categorical order


def front(cfg: Path, met: Path, n: int):
    m = np.load(met)
    mean, std, mx = m["mean"][:n], m["std"][:n], m["max"][:n]
    H = 1 - std / mean
    feas = mx <= BOUND
    e, h = np.where(feas, mean / BOUND, -1), np.where(feas, H, -1)
    nd = pareto_mask_fast(e, h) & feas
    o = np.argsort(e[nd])
    return e[nd][o], h[nd][o]


def main() -> None:
    ap = common_args(__doc__)
    ap.add_argument("--results", default=str(DEFAULT_OUT))
    ap.add_argument("--seeds", default="1,2,3,4,5")
    args = ap.parse_args()
    res, out = Path(args.results), Path(args.out)
    seeds = [int(s) for s in args.seeds.split(",")]
    conv = pd.read_csv(out / "r3_3_convergence.csv")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6})
    fig, axs = plt.subplots(2, 2, figsize=(17 / 2.54, 13 / 2.54), constrained_layout=True)
    for col, key in enumerate(("UV", "GRN")):
        ax = axs[0, col]
        runs = [("published, seed 42", res / "mc_published_seed42_5000000_configs.npz", 5_000_000)]
        runs += [(f"seed {s}", res / f"mc_seed{s}_10000000_configs.npz", 5_000_000) for s in seeds]
        for i, (label, cpath, n) in enumerate(runs):
            mpath = Path(str(cpath).replace("_configs.npz", f"_metrics_{key}.npz"))
            if not mpath.exists():
                continue
            e, h = front(cpath, mpath, n)
            ax.plot(e, h, "-o", ms=2.5, lw=1.2, color=COLORS[i], label=label)
        ax.set_xlabel(r"$E_\mathrm{mean}/E_\mathrm{bound}$")
        ax.set_ylabel("homogeneity $H$")
        ax.set_title(f"{'ab'[col]}) {'365' if key == 'UV' else '530'} nm, 5 million candidates", fontsize=9, loc="left")
        ax2 = axs[1, col]
        d = conv[(conv.led == key) & conv.run.str.startswith("seed")]
        for i, (label, g) in enumerate(d.groupby("run", sort=False), start=1):
            ax2.plot(g.n / 1e6, g.best_score, "-o", ms=3, lw=1.2, color=COLORS[i], label=label)
        pub = conv[(conv.led == key) & conv.run.str.startswith("published")]
        ax2.plot(pub.n / 1e6, pub.best_score, "D", ms=5, color=COLORS[0], label="published, seed 42")
        ax2.set_xscale("log")
        ax2.set_xticks([0.5, 1, 2.5, 5, 10])
        ax2.set_xticklabels(["0.5", "1", "2.5", "5", "10"])
        ax2.set_xlabel("distinct candidates / million")
        ax2.set_ylabel(r"best $H\cdot E_\mathrm{mean}/E_\mathrm{bound}$")
        ax2.set_title(f"{'cd'[col]}) {'365' if key == 'UV' else '530'} nm", fontsize=9, loc="left")
    axs[0, 0].legend(fontsize=7, frameon=False, loc="lower left")
    fig.savefig(out / "r3_3_convergence.png", dpi=300)
    print(f"written {out / 'r3_3_convergence.png'}")


if __name__ == "__main__":
    main()
