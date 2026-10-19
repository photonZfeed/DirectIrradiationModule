"""Check that MonteCarloSampler(seed=42) reproduces the published 5 million-sample set.

The published sample of the manuscript (``results/sampled_configs/MonteCarlo_1-16_leds_5-15_cm_
5000000_samples.npz``, written by the sampler under its former name
VariableLEDHeightMonteCarloGenerator) is compared with the first 5 million distinct
configurations of ``r3_mc_runs.py sample --seed 42`` (sequential mode, duplicates removed in the
order of the first draw, i.e. what ``MonteCarloSampler.generate_unique_samples`` returns):

1. row by row with the sampled set;
2. as a set with the configurations of both ``simulation_results_MonteCarlo_..._5000000_samples.npz``
   files (their rows are in a different order);
3. the recomputed mean, std and max irradiance with the values stored in those files, matched by
   configuration.

The published files are read with ``extract_mc_configs.py`` beforehand (positions only, no
unpickling); their file names are passed on the command line.

Usage: python check_seed42_published.py <results dir> <sampled.npz> <UV sim.npz> <UV sim extracted.npz>
                                        <GRN sim.npz> <GRN sim extracted.npz>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from r3_mc_runs import config_keys

N_PUB = 5_000_000


def keys_of(npz: Path, n: int | None = None) -> np.ndarray:
    d = np.load(npz)
    sl = slice(0, n)
    return config_keys(d["pos"][sl], d["z"][sl].astype(np.int16))


def main(results: Path, sampled: Path, sims: dict[str, tuple[Path, Path]]) -> None:
    ours = keys_of(results / "mc_seed42_10000000_configs.npz", N_PUB)
    pub = keys_of(sampled)
    print(f"sampled set: {len(pub)} rows; identical to the first {N_PUB} distinct draws of seed 42 "
          f"row by row: {bool(np.array_equal(ours, pub))} ({int(np.sum(ours != pub))} rows differ)")
    order = np.argsort(ours, kind="stable")
    ours_sorted = ours[order]
    for key, (sim, extracted) in sims.items():
        k = keys_of(extracted)
        loc = np.clip(np.searchsorted(ours_sorted, k), 0, len(ours) - 1)
        found = ours_sorted[loc] == k
        print(f"{key} simulation file: {len(k)} rows, all in the seed-42 set: {bool(found.all())}, "
              f"distinct: {len(np.unique(k)) == len(k)}")
        row = order[loc]                                   # row of every published entry in the seed-42 run
        m = np.load(results / f"mc_seed42_10000000_metrics_{key}.npz")
        ref = np.load(sim, allow_pickle=False)
        out = []
        for a, b in (("mean", "means"), ("std", "stds"), ("max", "max")):
            x, y = m[a][row], np.asarray(ref[b], dtype=float)
            out.append(f"{a} max rel. dev. {np.max(np.abs(x - y) / np.abs(y)):.2e}")
        print(f"  stored vs recomputed metrics: " + ", ".join(out))


if __name__ == "__main__":
    a = [Path(x) for x in sys.argv[1:]]
    main(a[0], a[1], {"UV": (a[2], a[3]), "GRN": (a[4], a[5])})
