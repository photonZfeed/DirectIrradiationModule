"""Reviewer 3, comment 3: convergence of the Monte-Carlo Pareto optimisation.

Uses the runs produced by ``r3_mc_runs.py`` (default seeds 42, 1, 2 and 3; 10 million distinct
configurations each, drawn with the repository's ``MonteCarloSampler``). The first 5 million
distinct configurations of seed 42 are the published sample of the manuscript
(``check_seed42_published.py``). For every seed and every sample size n in {0.5, 1, 2.5, 5, 10}
million (the first n distinct configurations of a run, i.e. nested samples), it reports for both
LED types at E_bound = 150 W m-2:

* the configuration that maximises H * E_mean / E_bound and whether it is the published UV3/GRN3,
* the three configurations that the manuscript's procedure exports for ray tracing,
* the hypervolume of the feasible non-dominated set in (E_mean / E_bound, H), reference (0, 0),
* the additive epsilon indicator against the reference front (union of all runs): the smallest
  shift by which the run's front weakly dominates every reference point.

The configuration with the highest H * E_mean / E_bound of every run and sample size is written
with its LED positions to ``r3_3_best_configurations.csv`` (input of Figure S19).

Optional (``--local-search``, not part of the publication): a best-improvement local search
started from the 20 best configurations of every run estimates how far the sampled optimum is
from a local optimum of the full (unsampled) design space.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from r2_common import DEFAULT_OUT, FieldBank, common_args, hypervolume_2d, index_to_xy, pareto_mask_fast, position_index, repo_path

SIZES = [500_000, 1_000_000, 2_500_000, 5_000_000, 10_000_000]
BOUND = 150.0
LED_NAME = {"UV": "LST1-01G01-UV01-00", "GRN": "LST1-01F06-GRN1-00"}
CENTER = 84


def run_list(res: Path, seeds: list[int]) -> list[tuple[str, int, Path, list[int]]]:
    """(label, seed, configuration file, nested sample sizes) of every run, in the order of ``seeds``."""
    return [(f"seed {s}", s, res / f"mc_seed{s}_10000000_configs.npz", SIZES) for s in seeds]


def published_index_key(repo: Path, key: str, k: int):
    prefix = "UV" if key == "UV" else "GRN"
    d = json.loads((repo / f"results/sampled_configs/near_optimal_configs/{LED_NAME[key]}/{prefix}{k}.json").read_text())
    xy = np.asarray(d["(x, y) Positions in [cm]"], dtype=float)
    return tuple(np.sort(position_index(xy[:, 0], xy[:, 1]))), int(d["Height in [cm]"])


def cfg_key(pos_row, z):
    return tuple(np.sort(pos_row[pos_row >= 0])), int(z)


def selection(mean, std, mx):
    H = 1 - std / mean
    feas = mx <= BOUND
    ratio = np.where(feas, mean / BOUND, 0.0)
    Hf = np.where(feas, H, 0.0)
    front = []
    for t in np.linspace(0, 1, 100):
        i = int(np.argmax(t * ratio + (1 - t) * Hf))
        if i not in front:
            front.append(i)
    best = int(np.argmax(ratio * Hf))
    if best not in front:
        # the product optimum is not hit by any of the 100 t values (the repository would export
        # nothing); insert it at its place along the front, which is ordered by increasing E_mean
        front = sorted(front + [best], key=lambda i: ratio[i])
    b = front.index(best)
    return best, front[max(0, b - 2): b + 1], ratio, Hf, feas


def eps_indicator(run_e, run_h, ref_e, ref_h):
    """Additive epsilon (maximisation): max_r min_p max(r_e - p_e, r_h - p_h)."""
    d = np.maximum(ref_e[:, None] - run_e[None, :], ref_h[:, None] - run_h[None, :])
    return float(np.max(np.min(d, axis=1)))


def pair_class(p):
    ix, iy = divmod(int(p), 13)
    if ix == 6 or iy == 6:
        return "axis"
    return "Q13" if (ix > 6) == (iy > 6) else "Q24"


def local_search(bank: FieldBank, start_pos, start_z, max_iter=50):
    """Best-improvement search: swap one point-symmetric pair for another of the same class, or move z by 1 cm."""
    cand = {c: sorted({min(p, 168 - p) for p in range(169) if p != CENTER and pair_class(p) == c}) for c in ("Q13", "Q24", "axis")}

    def score(configs, zs):
        pos = np.full((len(configs), 16), -1, dtype=np.int16)
        for i, c in enumerate(configs):
            pos[i, : len(c)] = c
        m, s, x = bank.evaluate(pos, np.asarray(zs))
        return np.where(x <= BOUND, (1 - s / m) * m / BOUND, 0.0), m, 1 - s / m

    cur = sorted(int(p) for p in start_pos if p >= 0)
    z = int(start_z)
    cur_score = score([cur], [z])[0][0]
    for _ in range(max_iter):
        reps = sorted({min(p, 168 - p) for p in cur if p != CENTER})
        neigh, zs = [], []
        for r in reps:
            for q in cand[pair_class(r)]:
                if q in reps:
                    continue
                new = sorted([p for p in cur if p not in (r, 168 - r)] + [q, 168 - q])
                neigh.append(new); zs.append(z)
        for dz in (-1, 1):
            if 5 <= z + dz <= 15:
                neigh.append(cur); zs.append(z + dz)
        sc, _, _ = score(neigh, zs)
        j = int(np.argmax(sc))
        if sc[j] <= cur_score + 1e-12:
            break
        cur, z, cur_score = neigh[j], zs[j], sc[j]
    _, m, h = score([cur], [z])
    return cur, z, float(cur_score), float(m[0]), float(h[0])


def main() -> None:
    ap = common_args(__doc__)
    ap.add_argument("--results", default=str(DEFAULT_OUT), help="folder with the sampled configurations and metrics")
    ap.add_argument("--seeds", default="42,1,2,3", help="comma-separated seeds (42: the published sample)")
    ap.add_argument("--local-search", action="store_true", help="also run the local search (not reported)")
    args = ap.parse_args()
    repo = repo_path(args.repo)
    res = Path(args.results)
    seeds = [int(s) for s in args.seeds.split(",")]
    runs = run_list(res, seeds)
    Path(args.out).mkdir(parents=True, exist_ok=True)
    rows, fronts, best_rows = [], {}, []
    for key in ("UV", "GRN"):
        pub = {f"{k}": published_index_key(repo, key, k) for k in (1, 2, 3)}
        pub_inv = {v: f"{'UV' if key == 'UV' else 'GRN'}{k}" for k, v in pub.items()}
        union_e, union_h, union_key = [], [], []
        per_run = []
        for label, seed, cpath, sizes in runs:
            mpath = Path(str(cpath).replace("_configs.npz", f"_metrics_{key}.npz"))
            if not (cpath.exists() and mpath.exists()):
                raise FileNotFoundError(f"missing {cpath} or {mpath}")
            c = np.load(cpath)
            m = np.load(mpath)
            pos, z = c["pos"], c["z"].astype(int)
            mean, std, mx = m["mean"], m["std"], m["max"]
            for n in sizes:
                sl = slice(0, n)
                best, exported, ratio, Hf, feas = selection(mean[sl], std[sl], mx[sl])
                nd = pareto_mask_fast(np.where(feas, ratio, -1.0), np.where(feas, Hf, -1.0)) & feas
                e_nd, h_nd = ratio[nd], Hf[nd]
                per_run.append((label, n, e_nd, h_nd))
                union_e.append(e_nd); union_h.append(h_nd)
                union_key += [cfg_key(pos[i], z[i]) for i in np.flatnonzero(nd)]
                bk = cfg_key(pos[best], z[best])
                bx, by = index_to_xy(np.asarray(bk[0]))
                best_rows.append(dict(led=key, run=label, seed=seed, n=n, N=len(bk[0]), z_cm=bk[1],
                                      E_mean=float(mean[best]), H=float(1 - std[best] / mean[best]),
                                      max_irradiance=float(mx[best]), score=float(ratio[best] * Hf[best]),
                                      is_published=pub_inv.get(bk, ""), sample_index=int(best),
                                      position_indices=" ".join(str(int(i)) for i in bk[0]),
                                      positions_cm="; ".join(f"{x:g} {y:g}" for x, y in zip(bx, by))))
                rows.append(dict(led=key, run=label, seed=seed, n=n, best_N=len(bk[0]), best_z=bk[1],
                                 best_E_mean=float(mean[best]), best_H=float(1 - std[best] / mean[best]),
                                 best_score=float(ratio[best] * Hf[best]), best_is=pub_inv.get(bk, "other"),
                                 exported=" / ".join(pub_inv.get(cfg_key(pos[i], z[i]),
                                                                 f"N{(pos[i] >= 0).sum()} z{z[i]} E{mean[i]:.1f} H{1 - std[i] / mean[i]:.3f}")
                                                     for i in exported),
                                 hypervolume=hypervolume_2d(e_nd, h_nd), n_nondominated=int(nd.sum())))
        ue, uh = np.concatenate(union_e), np.concatenate(union_h)
        um = pareto_mask_fast(ue, uh)
        ref_e, ref_h = ue[um], uh[um]
        fronts[key] = (ref_e, ref_h)
        hv_ref = hypervolume_2d(ref_e, ref_h)
        for r, (label, n, e_nd, h_nd) in zip([r for r in rows if r["led"] == key], per_run):
            r["eps_vs_reference_front"] = eps_indicator(e_nd, h_nd, ref_e, ref_h)
            r["hv_shortfall_pct"] = 100 * (1 - r["hypervolume"] / hv_ref)
        if not args.local_search:
            continue
        # local search from the 20 best product scores among all runs (largest run of every seed)
        bank = FieldBank(repo, key)
        starts = {}
        for label, seed, cpath, sizes in runs:
            mpath = Path(str(cpath).replace("_configs.npz", f"_metrics_{key}.npz"))
            if not (cpath.exists() and mpath.exists()):
                continue
            c, m = np.load(cpath), np.load(mpath)
            H = 1 - m["std"] / m["mean"]
            sc = np.where(m["max"] <= BOUND, H * m["mean"] / BOUND, 0.0)
            for i in np.argsort(-sc)[:20]:
                starts[cfg_key(c["pos"][i], c["z"][i])] = sc[i]
        top = sorted(starts.items(), key=lambda kv: -kv[1])
        ls = []
        for (p, zz), s0 in top:
            cfg, z1, s1, m1, h1 = local_search(bank, np.asarray(p), zz)
            ls.append(dict(led=key, start_N=len(p), start_z=zz, start_score=s0, end_N=len(cfg), end_z=z1,
                           end_score=s1, end_E_mean=m1, end_H=h1, gain_pct=100 * (s1 / s0 - 1),
                           end_is=pub_inv.get((tuple(cfg), z1), "other")))
        pd.DataFrame(ls).to_csv(Path(args.out) / f"r3_3_local_search_{key}.csv", index=False, float_format="%.6g")
        print(pd.DataFrame(ls).to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    df = pd.DataFrame(rows)
    out = Path(args.out)
    df.to_csv(out / "r3_3_convergence.csv", index=False, float_format="%.6g")
    pd.DataFrame(best_rows).to_csv(out / "r3_3_best_configurations.csv", index=False, float_format="%.6g")
    np.savez_compressed(out / "r3_3_reference_fronts.npz", **{f"{k}_E_ratio": v[0] for k, v in fronts.items()},
                        **{f"{k}_H": v[1] for k, v in fronts.items()})
    with pd.option_context("display.width", 300, "display.max_columns", 30, "display.max_rows", 200):
        print(df.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
