"""Reviewer 3, comments 2 and 3: Monte-Carlo sample generation and evaluation.

Stage 1 (``sample``): draws N configurations with the repository's own ``MonteCarloSampler``
(sequential mode, so the result depends on the seed only and not on the number of CPU cores;
uniform LED-count sampling, 1-16 LEDs, heights 5-15 cm in 1 cm steps, exactly as in the
manuscript) and stores them compactly as position indices.

Stage 2 (``evaluate``): mean, std and max irradiance of every configuration for both LED types
with ``r2_common.FieldBank`` (validated against the repository's stored results).

Every stage skips work whose output file already exists, so the script can be re-run after an
interruption.

Example:
    python r3_mc_runs.py sample   --repo <repo> --out <dir> --seed 1 --n 10000000
    python r3_mc_runs.py evaluate --repo <repo> --out <dir> --configs <dir>/mc_seed1_10000000_configs.npz
"""
from __future__ import annotations

import contextlib
import io
import time
from pathlib import Path

import numpy as np

from r2_common import FieldBank, GRID_KW, HALF, common_args, import_repo, repo_path


def config_keys(pos: np.ndarray, z: np.ndarray) -> np.ndarray:
    """One structured key per configuration: 169-bit position mask (3 x uint64) plus height."""
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


def sample(repo: Path, out: Path, seed: int, n: int, chunk: int = 500_000) -> Path:
    """Draw configurations until ``n`` distinct ones exist; keep them in order of first draw.

    The draws and the deduplication are those of ``MonteCarloSampler.generate_unique_samples``
    (same random stream, first occurrence kept); drawing in chunks does not change the stream.
    For seed 42 the first 5 million entries are the published sample of the manuscript, row by
    row (``check_seed42_published.py``). The first k entries of the result are exactly the sample
    that "draw until k distinct configurations" would give, which makes the prefixes usable as
    smaller sample sizes.
    """
    dst = out / f"mc_seed{seed}_{n}_configs.npz"
    if dst.exists():
        print(f"exists: {dst}")
        return dst
    import_repo(repo)
    from monte_carlo_sampler import MonteCarloSampler
    from utils.grid import Grid
    sampler = MonteCarloSampler(Grid(**GRID_KW), led_min=1, led_max=16, height_min=5, height_max=15,
                                height_step=1, seed=seed, sampling_mode="uniform", verbose=False)
    pos_parts, z_parts, key_parts = [], [], []
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
        pos_parts.append(pos[new]); z_parts.append(z[new]); key_parts.append(keys[new])
        seen = np.sort(np.concatenate([seen, keys[new]]))
        n_unique += len(new)
        print(f"seed {seed}: {n_drawn} drawn, {n_unique} distinct ({time.time() - t0:.0f} s)", flush=True)
    pos = np.concatenate(pos_parts)[:n]
    z = np.concatenate(z_parts)[:n]
    n_leds = (pos >= 0).sum(axis=1).astype(np.int8)
    np.savez_compressed(dst, pos=pos, z=z, n_leds=n_leds, seed=seed, n_drawn=n_drawn)
    return dst


def evaluate(repo: Path, out: Path, configs: Path) -> list[Path]:
    d = np.load(configs)
    pos, z = d["pos"], d["z"].astype(int)
    written = []
    for key in ("UV", "GRN"):
        dst = out / configs.name.replace("_configs.npz", f"_metrics_{key}.npz")
        if dst.exists():
            print(f"exists: {dst}")
            written.append(dst)
            continue
        t0 = time.time()
        bank = FieldBank(repo, key)
        mean, std, mx = bank.evaluate(pos, z)
        np.savez_compressed(dst, mean=mean, std=std, max=mx)
        print(f"{key}: {len(pos)} configurations evaluated in {time.time() - t0:.0f} s -> {dst}", flush=True)
        written.append(dst)
    return written


def main() -> None:
    ap = common_args(__doc__)
    ap.add_argument("stage", choices=["sample", "evaluate"])
    ap.add_argument("--seed", type=int)
    ap.add_argument("--n", type=int, default=10_000_000)
    ap.add_argument("--configs")
    args = ap.parse_args()
    repo = repo_path(args.repo)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.stage == "sample":
        sample(repo, out, args.seed, args.n)
    else:
        evaluate(repo, out, Path(args.configs))


if __name__ == "__main__":
    main()
