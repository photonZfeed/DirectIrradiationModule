"""
Influence of the irradiance bound E_bound on the configuration selection (Supporting Information S11).

In :class:`optimizer.Optimizer` the bound plays two roles: it normalises the mean irradiance in the
objective (E_mean / E_bound) and it excludes configurations whose maximum local irradiance exceeds
it. The normalisation changes neither the Pareto-optimal set nor the configuration with the highest
H E_mean / E_bound, so only the exclusion can change the selection. This module repeats the
selection of the manuscript for a range of bounds and reports the highest local irradiance of the
investigated configurations.

The 5 million Monte-Carlo candidates of the manuscript are the first 5 million distinct
configurations of seed 42 in ``results/geometric/monte_carlo_convergence/``
(see :mod:`analysis.mc_sampling`; checked row by row against the published sample in
``tests/test_seed42_published.py``).

**Usage Example**::

    from analysis.irradiance_bound import irradiance_bound_table, max_irradiance_8_leds
    table_s5 = irradiance_bound_table()
    maxima = max_irradiance_8_leds()
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

from analysis.compare_methods import read_configuration
from analysis.fast_evaluation import (E_BOUND, GRID_KW, LEDS, PREFIX, RESOLUTION, FieldBank, configs_from_object_array,
                                      configuration_file, position_index, select_configurations)
from analysis.mc_sampling import CONVERGENCE_PATH, load_run
from geometric_model import GeometricModel
from utils.grid import Grid

#: Irradiance bounds of Table S5 in W m-2 (``numpy.inf``: no bound).
BOUNDS = [100, 110, 120, 135, 150, 175, 200, 250, 300, np.inf]
#: Number of candidates of the manuscript.
N_CANDIDATES = 5_000_000
#: Seed of the published sample.
SEED = 42
#: Systematic 8-LED sample at 13 cm (Figure 4 of the manuscript).
SYSTEMATIC_SAMPLES = "sampled_configs/systematic_8_leds_13_cm_482931_samples.npz"


def find_configuration(pos: np.ndarray, z: np.ndarray, positions: list, height: float) -> int:
    """
    Row of a configuration in a sample.

    :param pos: Position indices of the sample, shape (n, 16), padded with -1.
    :param z: Heights of the sample in cm.
    :param positions: LED positions (x, y) in cm of the configuration to find.
    :param height: Height in cm of the configuration to find.
    :returns: Row index, or -1 if the configuration is not in the sample.
    :rtype: int
    """
    xy = np.asarray(positions, dtype=float)
    target = np.sort(position_index(xy[:, 0], xy[:, 1]))
    cand = np.flatnonzero((z == int(height)) & ((pos >= 0).sum(axis=1) == len(target)))
    for i in cand:
        p = np.sort(pos[i][pos[i] >= 0])
        if np.array_equal(p, target):
            return int(i)
    return -1


def max_irradiance_8_leds(results_path: str = "results/") -> pd.DataFrame:
    """
    Highest local irradiance of all 482,931 symmetric 8-LED configurations at 13 cm.

    :param results_path: Path of the results directory.
    :type results_path: str, optional
    :returns: One row per LED type with the largest local irradiance of all configurations and the
        number of configurations above the bound of the manuscript.
    :rtype: pandas.DataFrame
    """
    d = np.load(os.path.join(results_path, SYSTEMATIC_SAMPLES), allow_pickle=True)
    height = int(d["height"])
    pos = configs_from_object_array(d["configurations"])
    z = np.full(len(pos), height, dtype=int)
    rows = []
    for key, led in LEDS.items():
        mean, _, mx = FieldBank(led, heights=(height,)).evaluate(pos, z)
        i = int(np.nanargmax(mx))
        rows.append(dict(led=key, n_configs=len(pos), n_leds=int(d["led_count"]), z_cm=height,
                         E_max_all_configs=float(mx[i]), E_mean_of_that_config=float(mean[i]),
                         n_above_bound=int((mx > E_BOUND).sum())))
    return pd.DataFrame(rows)


def max_irradiance_selected(results_path: str = "results/") -> pd.DataFrame:
    """
    Mean and highest local irradiance of UV1-UV3 and GRN1-GRN3 according to the geometric model.

    :param results_path: Path of the results directory.
    :type results_path: str, optional
    :rtype: pandas.DataFrame
    """
    rows = []
    for key, led in LEDS.items():
        model = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        for k in (1, 2, 3):
            name = f"{PREFIX[key]}{k}"
            positions, z = read_configuration(configuration_file(key, name, results_path))
            E = model.simulate(positions, z)
            rows.append(dict(led=key, config=name, n_leds=len(positions), z_cm=z, E_mean=E.mean(),
                             H=1 - E.std() / E.mean(), E_max=E.max()))
    return pd.DataFrame(rows)


def irradiance_bound_table(convergence_path: str = CONVERGENCE_PATH, results_path: str = "results/") -> pd.DataFrame:
    """
    Compute Table S5: the selection of the manuscript repeated for different irradiance bounds.

    :param convergence_path: Folder with the Monte-Carlo samples and metrics of
        :mod:`analysis.mc_sampling` (seed 42).
    :type convergence_path: str, optional
    :param results_path: Path of the results directory (configuration files).
    :type results_path: str, optional
    :returns: One row per LED type and bound with the best configuration (number of LEDs, height,
        E_mean, H, maximum local irradiance), whether it is one of UV1-UV3 / GRN1-GRN3, and the
        three configurations selected for ray tracing.
    :rtype: pandas.DataFrame
    :raises FileNotFoundError: If the Monte-Carlo files of seed 42 are missing.
    """
    rows = []
    for key in LEDS:
        pos, z, mean, std, mx = load_run(convergence_path, SEED, key, N_CANDIDATES)
        n_leds = (pos >= 0).sum(axis=1)
        published = {}
        for k in (1, 2, 3):
            name = f"{PREFIX[key]}{k}"
            published[find_configuration(pos, z, *read_configuration(configuration_file(key, name, results_path)))] = name
        published.pop(-1, None)
        for bound in BOUNDS:
            best, selected, _, _ = select_configurations(mean, std, mx, bound)
            h = 1 - std / mean
            rows.append(dict(led=key, E_bound=bound, best_N=int(n_leds[best]), best_z_cm=int(z[best]),
                             best_E_mean=float(mean[best]), best_H=float(h[best]), best_E_max=float(mx[best]),
                             best_is=published.get(best, "other"),
                             selected=" / ".join(published.get(i, f"N{n_leds[i]} z{z[i]} E{mean[i]:.1f} H{h[i]:.3f}")
                                                 for i in selected)))
    return pd.DataFrame(rows)
