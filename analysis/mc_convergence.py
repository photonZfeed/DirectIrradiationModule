"""
Convergence of the Monte-Carlo optimization (Supporting Information S12, Figures S18 and S19).

For every LED type, seed (42, 1, 2, 3) and sample size n (the first 0.5, 1, 2.5, 5 and 10 million
distinct configurations of a run of :mod:`analysis.mc_sampling`) the optimization of the manuscript
is evaluated with E_bound = 150 W m-2:

* :func:`best_configurations`: the configuration with the highest H E_mean / E_bound, its LED
  positions and whether it is UV3 / GRN3 (Figure S18c, d and Figure S19);
* :func:`fronts`: the non-dominated configurations of the first 5 million candidates of every seed
  (Figure S18a, b);
* :func:`front_deviation`: how closely the fronts of the seeds agree. ``dH(E)`` is the spread
  (maximum minus minimum over the seeds) of the homogeneity that each front reaches at the same mean
  irradiance E, i.e. the highest H of a front point with E_mean >= E. It is evaluated on 4000
  equidistant values of E over the range covered by all fronts, only where every front has
  H >= 0.5. This cut removes the steep end of the 530 nm front, where H drops from 0.5 to below 0
  within a few hundredths of a W m-2.

**Usage Example**::

    from analysis.mc_convergence import convergence_tables
    best, fronts_5m, deviation = convergence_tables()
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from analysis.compare_methods import read_configuration
from analysis.fast_evaluation import E_BOUND, LEDS, PREFIX, configuration_file, index_to_xy, nondominated_front
from analysis.irradiance_bound import find_configuration
from analysis.mc_sampling import CONVERGENCE_PATH, SEEDS, load_run

#: Sample sizes of the study.
SIZES = (500_000, 1_000_000, 2_500_000, 5_000_000, 10_000_000)
#: Sample size of the fronts of Figure S18a, b.
SIZE_FRONT = 5_000_000
#: Only front segments where every front has at least this homogeneity are compared.
H_MIN = 0.5
#: Number of mean-irradiance values at which the fronts are compared.
N_COMPARE = 4000


def _h_at(front: tuple[np.ndarray, np.ndarray], e_query: np.ndarray) -> np.ndarray:
    """Highest H of a front point with E_mean / E_bound >= e_query (front sorted by increasing E)."""
    e, h = front
    idx = np.searchsorted(e, e_query, side="left")
    out = np.full(e_query.shape, np.nan)
    ok = idx < len(e)
    out[ok] = h[idx[ok]]
    return out


def convergence_tables(path: str = CONVERGENCE_PATH, seeds=SEEDS, sizes=SIZES,
                       results_path: str = "results/") -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Evaluate all runs: best configurations, fronts at 5 million candidates and front agreement.

    :param path: Folder with the samples and metrics of :mod:`analysis.mc_sampling`.
    :param seeds: Seeds of the runs.
    :param sizes: Sample sizes (nested prefixes of every run).
    :param results_path: Path of the results directory (configuration files of UV3 / GRN3).
    :returns: Table of the best configurations per LED type, seed and sample size; the fronts at
        5 million candidates in long format (``led``, ``seed``, ``E_ratio``, ``H``); the agreement
        of the fronts per LED type and sample size.
    :rtype: tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]
    :raises FileNotFoundError: If a file of a run is missing.
    """
    best_rows, front_rows, dev_rows = [], [], []
    for key in LEDS:
        name = f"{PREFIX[key]}3"
        selected = read_configuration(configuration_file(key, name, results_path))
        all_fronts = {}
        for seed in seeds:
            pos, z, mean, std, mx = load_run(path, seed, key, max(sizes))
            h = 1 - std / mean
            score = np.where(mx <= E_BOUND, h * mean / E_BOUND, 0.0)
            selected_row = find_configuration(pos, z, *selected)
            for n in sizes:
                best = int(np.argmax(score[:n]))
                p = np.sort(pos[best][pos[best] >= 0])
                x, y = index_to_xy(p)
                best_rows.append(dict(led=key, seed=seed, n=n, N=len(p), z_cm=int(z[best]),
                                      E_mean=float(mean[best]), H=float(h[best]), max_irradiance=float(mx[best]),
                                      score=float(score[best]), is_published=name if best == selected_row else "",
                                      sample_index=best, position_indices=" ".join(str(int(i)) for i in p),
                                      positions_cm="; ".join(f"{a:g} {b:g}" for a, b in zip(x, y))))
                all_fronts[(seed, n)] = nondominated_front(mean[:n], std[:n], mx[:n])
            e, hh = all_fronts[(seed, SIZE_FRONT)]
            front_rows += [dict(led=key, seed=seed, E_ratio=a, H=b) for a, b in zip(e, hh)]
        for n in sizes:
            fr = [all_fronts[(s, n)] for s in seeds]
            e_lo, e_hi = max(f[0].min() for f in fr), min(f[0].max() for f in fr)
            eq = np.linspace(e_lo, e_hi, N_COMPARE)
            hq = np.array([_h_at(f, eq) for f in fr])
            dh = hq.max(0) - hq.min(0)
            use = hq.min(0) >= H_MIN
            k = int(np.argmax(np.where(use, dh, -1.0)))
            dev_rows.append(dict(led=key, n=n, front_points=" / ".join(str(len(f[0])) for f in fr),
                                 E_mean_min=e_lo * E_BOUND, E_mean_max=e_hi * E_BOUND,
                                 dH_median=float(np.median(dh[use])), dH_p90=float(np.percentile(dh[use], 90)),
                                 dH_max=float(dh[k]), E_mean_at_dH_max=float(eq[k] * E_BOUND)))
    return pd.DataFrame(best_rows), pd.DataFrame(front_rows), pd.DataFrame(dev_rows)
