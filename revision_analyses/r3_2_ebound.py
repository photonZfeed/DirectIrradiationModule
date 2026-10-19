"""Reviewer 3, comment 2: sensitivity of the Pareto analysis to the irradiance bound E_bound.

In ``optimizer.Optimizer`` E_bound plays two roles:
  1. normalisation of the mean irradiance in the objective (E_mean / E_bound), and
  2. a hard constraint: configurations whose maximum local irradiance exceeds E_bound get a
     score of zero.
Role 1 cannot change which configurations are Pareto-optimal or which one maximises
H * E_mean / E_bound (a positive rescaling of one axis preserves dominance and the arg max of the
product). Only role 2 can, so the analysis re-runs the selection of the manuscript for a range of
bounds on the published 5 million-sample set and reports what changes.

The selection logic of ``visualization/pareto_plots_H_vs_E_ratio.process_led`` is reproduced:
t is swept over linspace(0, 1, 100), each t contributes the arg max of
t * E_mean/E_bound + (1 - t) * H (infeasible = 0), the best configuration maximises
H * E_mean / E_bound, and the three configurations exported for ray tracing are the best one and
the two Pareto points preceding it in t.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from r2_common import DEFAULT_OUT, common_args, hypervolume_2d, pareto_mask_fast, position_index, repo_path

BOUNDS = [100, 110, 120, 135, 150, 175, 200, 250, 300, np.inf]
NAMES = {"UV": ("LST1-01G01-UV01-00", "UV"), "GRN": ("LST1-01F06-GRN1-00", "GRN")}


def find_config(pos, z, json_path: Path) -> int:
    d = json.loads(json_path.read_text())
    xy = np.asarray(d["(x, y) Positions in [cm]"], dtype=float)
    target = np.sort(position_index(xy[:, 0], xy[:, 1]))
    zt = int(d["Height in [cm]"])
    cand = np.flatnonzero((z == zt) & ((pos >= 0).sum(axis=1) == len(target)))
    for i in cand:
        p = np.sort(pos[i][pos[i] >= 0])
        if np.array_equal(p, target):
            return int(i)
    return -1


def selection(mean, std, mx, bound):
    H = 1 - std / mean
    feas = mx <= bound
    ratio = np.where(feas, mean / (bound if np.isfinite(bound) else 150.0), 0.0)
    Hf = np.where(feas, H, 0.0)
    front = []
    for t in np.linspace(0, 1, 100):
        i = int(np.argmax(t * ratio + (1 - t) * Hf))
        if i not in front:
            front.append(i)
    best = int(np.argmax(ratio * Hf))
    b = front.index(best) if best in front else None
    exported = front[max(0, b - 2): b + 1] if b is not None else []
    return front, best, exported, feas, H


def main() -> None:
    ap = common_args(__doc__)
    ap.add_argument("--results", default=str(DEFAULT_OUT), help="folder with the mc_published_* files")
    args = ap.parse_args()
    repo = repo_path(args.repo)
    res = Path(args.results)
    cfg = np.load(res / "mc_published_seed42_5000000_configs.npz")
    pos, z, n_leds = cfg["pos"], cfg["z"].astype(int), cfg["n_leds"].astype(int)
    rows, fronts = [], []
    for key, (led_name, prefix) in NAMES.items():
        m = np.load(res / f"mc_published_seed42_5000000_metrics_{key}.npz")
        mean, std, mx = m["mean"], m["std"], m["max"]
        published = {f"{prefix}{k}": find_config(pos, z, repo / f"results/sampled_configs/near_optimal_configs/{led_name}/{prefix}{k}.json")
                     for k in (1, 2, 3)}
        inv = {v: k for k, v in published.items()}
        print(key, "published configurations found at sample indices:", published,
              "| max irradiance there:", {k: round(float(mx[v]), 1) for k, v in published.items()})
        for bound in BOUNDS:
            front, best, exported, feas, H = selection(mean, std, mx, bound)
            nd = pareto_mask_fast(np.where(feas, mean, -1.0), np.where(feas, H, -1.0)) & feas
            hv = hypervolume_2d(mean[nd], H[nd])
            row = dict(led=key, E_bound=bound, feasible_fraction=float(feas.mean()),
                       n_supported_front=len(front), n_nondominated=int(nd.sum()), hypervolume_W_m2=hv,
                       best_N=int(n_leds[best]), best_z_cm=int(z[best]), best_E_mean=float(mean[best]),
                       best_H=float(H[best]), best_E_max=float(mx[best]),
                       best_is=inv.get(best, "other"),
                       exported=" / ".join(inv.get(i, f"N{n_leds[i]} z{z[i]} E{mean[i]:.1f} H{H[i]:.3f}") for i in exported))
            rows.append(row)
            for i in np.flatnonzero(nd):
                fronts.append(dict(led=key, E_bound=bound, E_mean=float(mean[i]), H=float(H[i]),
                                   N=int(n_leds[i]), z_cm=int(z[i]), E_max=float(mx[i])))
    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "r3_2_ebound_sensitivity.csv", index=False, float_format="%.6g")
    pd.DataFrame(fronts).to_csv(out / "r3_2_nondominated_fronts.csv", index=False, float_format="%.6g")
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(df.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
