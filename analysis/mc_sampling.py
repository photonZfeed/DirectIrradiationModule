"""
Monte-Carlo samples of the convergence study (Supporting Information S12) and their evaluation.

The study draws 10 million distinct configurations (1-16 LEDs, heights 5-15 cm in 1 cm steps) for
each of the seeds 42, 1, 2 and 3 with the repository's :class:`utils.monte_carlo_sampler.MonteCarloSampler`
(uniform LED-count sampling, sequential mode, duplicates removed in the order of the first draw, as
in :meth:`MonteCarloSampler.generate_unique_samples`). With seed 42 the first 5 million
configurations are the published sample of the manuscript, row by row
(``tests/test_seed42_published.py``). The first n configurations of a run are exactly the sample
that drawing until n distinct configurations would give, so the prefixes of one run serve as the
smaller sample sizes.

Every configuration is evaluated for both LED types with :class:`analysis.fast_evaluation.FieldBank`.
The files are written to ``results/geometric/monte_carlo_convergence/`` (about 2.1 GB for four
seeds; also available from the Zenodo dataset):

* ``mc_seed<seed>_10000000_configs.npz``: ``pos`` (position indices, int16, padded with -1),
  ``z`` (height in cm), ``n_leds``, ``seed``, ``n_drawn``;
* ``mc_seed<seed>_10000000_metrics_<UV|GRN>.npz``: ``mean``, ``std`` and ``max`` of the irradiance
  in W m-2.

Sampling all four seeds takes about two hours on two cores; a sampler needs about 2 GB of memory at
the end of a run. Finished files are not recomputed, so an interrupted run can be restarted.

**Usage Example** (from the repository root)::

    python -m analysis.mc_sampling --seeds 42 1 2 3
"""
from __future__ import annotations

import argparse
import contextlib
import io
import os
import time

import numpy as np

from analysis.fast_evaluation import GRID_KW, HALF, LEDS, FieldBank
from utils.monte_carlo_sampler import MonteCarloSampler
from utils.grid import Grid

#: Folder of the samples and metrics of the convergence study.
CONVERGENCE_PATH = os.path.join("results", "geometric", "monte_carlo_convergence")
#: Seeds of the convergence study; 42 is the seed of the published sample.
SEEDS = (42, 1, 2, 3)
#: Number of distinct configurations per seed.
N_PER_SEED = 10_000_000


def configs_file(path: str, seed: int, n: int = N_PER_SEED) -> str:
    """Path of the configuration file of a run."""
    return os.path.join(path, f"mc_seed{seed}_{n}_configs.npz")


def metrics_file(path: str, seed: int, key: str, n: int = N_PER_SEED) -> str:
    """Path of the metrics file of a run for the LED type ``key`` (``"UV"`` or ``"GRN"``)."""
    return os.path.join(path, f"mc_seed{seed}_{n}_metrics_{key}.npz")


def config_keys(pos: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    One sortable key per configuration: the 169-bit position mask (3 x uint64) and the height.

    :param pos: Position indices of shape (n, 16), padded with -1.
    :param z: Heights in cm.
    :returns: Structured array with the fields ``w0``, ``w1``, ``w2`` and ``z``.
    :rtype: numpy.ndarray
    """
    words = np.zeros((len(pos), 3), dtype=np.uint64)
    for col in range(pos.shape[1]):
        p = pos[:, col].astype(np.int64)
        v = p >= 0
        w, b = np.divmod(np.where(v, p, 0), 64)
        bit = np.where(v, np.left_shift(np.uint64(1), b.astype(np.uint64)), np.uint64(0))
        for k in range(3):
            words[:, k] |= np.where(w == k, bit, np.uint64(0))
    keys = np.empty(len(pos), dtype=[("w0", "u8"), ("w1", "u8"), ("w2", "u8"), ("z", "i2")])
    keys["w0"], keys["w1"], keys["w2"], keys["z"] = words[:, 0], words[:, 1], words[:, 2], z
    return keys


def sample(seed: int, n: int = N_PER_SEED, path: str = CONVERGENCE_PATH, chunk: int = 500_000) -> str:
    """
    Draw configurations with :class:`MonteCarloSampler` until ``n`` distinct ones exist.

    The configurations are kept in the order of their first draw. Drawing in chunks does not change
    the random stream.

    :param seed: Random seed of the sampler.
    :param n: Number of distinct configurations.
    :param path: Output folder.
    :param chunk: Number of configurations drawn at once.
    :returns: Path of the configuration file (an existing file is kept).
    :rtype: str
    """
    dst = configs_file(path, seed, n)
    if os.path.exists(dst):
        print(f"exists: {dst}")
        return dst
    os.makedirs(path, exist_ok=True)
    sampler = MonteCarloSampler(Grid(**GRID_KW), led_min=1, led_max=16, height_min=5, height_max=15,
                                height_step=1, seed=seed, sampling_mode="uniform", verbose=False)
    pos_parts, z_parts = [], []
    seen = np.empty(0, dtype=[("w0", "u8"), ("w1", "u8"), ("w2", "u8"), ("z", "i2")])
    n_unique, n_drawn = 0, 0
    t0 = time.time()
    while n_unique < n:
        with contextlib.redirect_stdout(io.StringIO()):
            batch = sampler.generate_samples(chunk, use_parallel=False)
        if len(batch) != chunk:
            raise RuntimeError("sampler returned fewer configurations than requested")
        pos = np.full((chunk, 16), -1, dtype=np.int16)
        z = np.empty(chunk, dtype=np.int16)
        for j, s in enumerate(batch):
            c = np.asarray(s["configuration"], dtype=np.int16)      # integer grid units
            pos[j, : len(c)] = (c[:, 0] + HALF) * (2 * HALF + 1) + (c[:, 1] + HALF)
            z[j] = s["height"]
        n_drawn += chunk
        keys = config_keys(pos, z)
        _, first = np.unique(keys, return_index=True)            # first occurrence within the chunk
        first.sort()
        new = first[~np.isin(keys[first], seen)]                   # drop keys seen in earlier chunks
        pos_parts.append(pos[new])
        z_parts.append(z[new])
        seen = np.sort(np.concatenate([seen, keys[new]]))
        n_unique += len(new)
        print(f"seed {seed}: {n_drawn} drawn, {n_unique} distinct ({time.time() - t0:.0f} s)", flush=True)
    pos = np.concatenate(pos_parts)[:n]
    z = np.concatenate(z_parts)[:n]
    n_leds = (pos >= 0).sum(axis=1).astype(np.int8)
    np.savez_compressed(dst, pos=pos, z=z, n_leds=n_leds, seed=seed, n_drawn=n_drawn)
    return dst


def evaluate(seed: int, n: int = N_PER_SEED, path: str = CONVERGENCE_PATH) -> list[str]:
    """
    Mean, standard deviation and maximum irradiance of every configuration of a run, both LED types.

    :returns: Paths of the metrics files (existing files are kept).
    :rtype: list[str]
    """
    d = np.load(configs_file(path, seed, n))
    pos, z = d["pos"], d["z"].astype(int)
    written = []
    for key, led in LEDS.items():
        dst = metrics_file(path, seed, key, n)
        if not os.path.exists(dst):
            t0 = time.time()
            mean, std, mx = FieldBank(led).evaluate(pos, z)
            np.savez_compressed(dst, mean=mean, std=std, max=mx)
            print(f"seed {seed}, {key}: {len(pos)} configurations evaluated in {time.time() - t0:.0f} s", flush=True)
        written.append(dst)
    return written


def run(seeds=SEEDS, n: int = N_PER_SEED, path: str = CONVERGENCE_PATH) -> None:
    """Sample and evaluate all runs of the convergence study (finished files are skipped)."""
    for seed in seeds:
        sample(seed, n, path)
        evaluate(seed, n, path)


def available(seeds=SEEDS, path: str = CONVERGENCE_PATH, n: int = N_PER_SEED) -> bool:
    """True if the configuration and metrics files of all runs exist."""
    files = [configs_file(path, s, n) for s in seeds] + [metrics_file(path, s, k, n) for s in seeds for k in LEDS]
    return all(os.path.exists(f) for f in files)


def load_run(path: str, seed: int, key: str, n: int | None = None):
    """
    Read the first ``n`` configurations of a run and their metrics for one LED type.

    :param path: Folder of the run files.
    :param seed: Seed of the run.
    :param key: LED type, ``"UV"`` or ``"GRN"``.
    :param n: Number of configurations to read (all if None).
    :returns: Position indices, heights, mean, standard deviation and maximum irradiance.
    :rtype: tuple[numpy.ndarray, ...]
    :raises FileNotFoundError: If a file of the run is missing.
    """
    sl = slice(0, n)
    with np.load(configs_file(path, seed)) as c:
        pos, z = c["pos"][sl], c["z"][sl].astype(int)
    with np.load(metrics_file(path, seed, key)) as m:
        mean, std, mx = (np.asarray(m[k][sl], dtype=float) for k in ("mean", "std", "max"))
    return pos, z, mean, std, mx


def main() -> None:
    """Command line interface: sample and evaluate the given seeds."""
    ap = argparse.ArgumentParser(description="Monte-Carlo runs of the convergence study (SI Section S12).")
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    ap.add_argument("--n", type=int, default=N_PER_SEED, help="distinct configurations per seed")
    ap.add_argument("--path", default=CONVERGENCE_PATH, help="output folder")
    args = ap.parse_args()
    run(args.seeds, args.n, args.path)


if __name__ == "__main__":
    main()
