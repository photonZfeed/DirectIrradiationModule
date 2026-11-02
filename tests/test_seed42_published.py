"""
Seed 42 of the convergence study reproduces the published 5 million-sample Monte-Carlo set.

The published sample of the manuscript (``results/sampled_configs/MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz``)
is compared with the first 5 million distinct configurations of seed 42 in
``results/geometric/monte_carlo_convergence/`` (:mod:`analysis.mc_sampling`):

1. row by row with the sampled set;
2. as a set with the configurations of both
   ``results/geometric/<LED>/simulation_results_MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz``
   (their rows are stored in a different order);
3. the mean, standard deviation and maximum irradiance of the convergence study with the values
   stored in those files, matched by configuration.

The published files and the convergence data are part of the Zenodo dataset and not of the git
repository; the test is skipped if they are missing. Reading the pickled configurations of the
three published files takes several minutes.

Run from the repository root with ``python -m pytest tests`` or ``python tests/test_seed42_published.py``.
"""
import os
import pickletools
import sys
import zipfile

if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402
import numpy.lib.format as npformat  # noqa: E402
import pytest  # noqa: E402

from analysis.fast_evaluation import LEDS, position_index  # noqa: E402
from analysis.mc_sampling import CONVERGENCE_PATH, config_keys, configs_file, metrics_file  # noqa: E402

N_PUBLISHED = 5_000_000
SAMPLED = "results/sampled_configs/MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz"
SIMULATED = "results/geometric/{}/simulation_results_MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz"
TOLERANCE = 1e-6


def read_pickled_configurations(path: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Read the LED positions of a published Monte-Carlo file without unpickling it.

    The files store the configurations as a pickled object array (members ``configs`` and ``z`` of
    the simulation results, ``configurations`` and ``height`` of the sampled set), which needs more
    than 5 GB of memory as Python objects. The pickle opcodes are walked as a stream instead.

    :returns: Position indices of shape (n, 16), padded with -1, and the heights in cm.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    zf = zipfile.ZipFile(path)
    names = set(zf.namelist())
    cfg_member, z_member = ("configs.npy", "z.npy") if "configs.npy" in names else ("configurations.npy", "height.npy")
    with zf.open(cfg_member) as f:
        version = npformat.read_magic(f)
        shape, _, _ = (npformat.read_array_header_1_0(f) if version == (1, 0) else npformat.read_array_header_2_0(f))
        n = shape[0]
        xs = np.full((n, 16), np.nan, dtype=np.float32)
        ys = np.full((n, 16), np.nan, dtype=np.float32)
        count = np.zeros(n, dtype=np.int8)
        cur = -2                      # -2: before the outer list, -1: outer list opened
        pending: list[float] = []
        for op, arg, _ in pickletools.genops(f):
            if op.name == "EMPTY_LIST":
                cur += 1
                pending = []
            elif op.name == "BINFLOAT" and cur >= 0:
                pending.append(arg)
                if len(pending) == 2:
                    k = count[cur]
                    xs[cur, k], ys[cur, k] = pending
                    count[cur] = k + 1
                    pending = []
            elif op.name == "STOP":
                break
        if cur != n - 1:
            raise RuntimeError(f"parsed {cur + 1} configurations, expected {n}")
    with zf.open(z_member) as f:
        z = np.load(f).astype(np.int16)
    pos = np.full((n, 16), -1, dtype=np.int16)
    valid = ~np.isnan(xs)
    pos[valid] = position_index(xs[valid], ys[valid])
    return pos, z


def _required_files() -> list[str]:
    return ([SAMPLED, configs_file(CONVERGENCE_PATH, 42)] + [SIMULATED.format(led.name) for led in LEDS.values()]
            + [metrics_file(CONVERGENCE_PATH, 42, key) for key in LEDS])


def _seed42_keys() -> np.ndarray:
    with np.load(configs_file(CONVERGENCE_PATH, 42)) as d:
        return config_keys(d["pos"][:N_PUBLISHED], d["z"][:N_PUBLISHED].astype(np.int16))


def check_sampled_set() -> int:
    """Number of rows in which the published sampled set differs from seed 42."""
    pos, z = read_pickled_configurations(SAMPLED)
    return int(np.sum(_seed42_keys() != config_keys(pos, z)))


def check_simulation_file(key: str) -> dict:
    """Compare a published simulation file with seed 42 (as a set, and the stored metrics)."""
    ours = _seed42_keys()
    order = np.argsort(ours, kind="stable")
    ours_sorted = ours[order]
    pos, z = read_pickled_configurations(SIMULATED.format(LEDS[key].name))
    k = config_keys(pos, z)
    loc = np.clip(np.searchsorted(ours_sorted, k), 0, len(ours) - 1)
    found = ours_sorted[loc] == k
    row = order[loc]                                    # row of every published entry in the seed-42 run
    out = dict(rows=len(k), all_in_seed42=bool(found.all()), distinct=len(np.unique(k)) == len(k))
    with np.load(metrics_file(CONVERGENCE_PATH, 42, key)) as m, np.load(SIMULATED.format(LEDS[key].name)) as ref:
        for a, b in (("mean", "means"), ("std", "stds"), ("max", "max")):
            x, y = m[a][row], np.asarray(ref[b], dtype=float)
            out[f"{a}_max_rel_dev"] = float(np.max(np.abs(x - y) / np.abs(y)))
    return out


@pytest.mark.skipif(not all(os.path.exists(f) for f in _required_files()),
                    reason="published Monte-Carlo files or convergence data (Zenodo) not present")
def test_sampled_set_row_by_row():
    assert check_sampled_set() == 0


@pytest.mark.skipif(not all(os.path.exists(f) for f in _required_files()),
                    reason="published Monte-Carlo files or convergence data (Zenodo) not present")
@pytest.mark.parametrize("key", list(LEDS))
def test_simulation_files(key):
    r = check_simulation_file(key)
    assert r["rows"] == N_PUBLISHED and r["all_in_seed42"] and r["distinct"]
    assert max(r[f"{a}_max_rel_dev"] for a in ("mean", "std", "max")) < TOLERANCE


if __name__ == "__main__":
    print("sampled set, rows that differ from seed 42:", check_sampled_set())
    for led_key in LEDS:
        print(led_key, check_simulation_file(led_key))
