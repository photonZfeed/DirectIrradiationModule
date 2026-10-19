"""Extract the LED positions of the published 5 million-sample Monte-Carlo set without unpickling.

``simulation_results_MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz`` (members ``configs`` and
``z``) and the sampled set ``MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz`` (members
``configurations`` and ``height``) store the configurations as a pickled object array
(1.2 GB uncompressed, > 5 GB as Python objects).
This script walks the pickle opcodes as a stream and writes a compact file with the
position indices (int16, padded with -1), the LED count and the height of every sample.

Usage: python extract_mc_configs.py <simulation_results_...npz | MonteCarlo_..._samples.npz> <output.npz>
"""
from __future__ import annotations

import pickletools
import sys
import zipfile

import numpy as np
import numpy.lib.format as fmt

from r2_common import position_index


def main(src: str, dst: str) -> None:
    zf = zipfile.ZipFile(src)
    names = set(zf.namelist())
    cfg_member, z_member = ("configs.npy", "z.npy") if "configs.npy" in names else ("configurations.npy", "height.npy")
    with zf.open(cfg_member) as f:
        version = fmt.read_magic(f)
        shape, _, _ = (fmt.read_array_header_1_0(f) if version == (1, 0) else fmt.read_array_header_2_0(f))
        n = shape[0]
        xs = np.full((n, 16), np.nan, dtype=np.float32)
        ys = np.full((n, 16), np.nan, dtype=np.float32)
        count = np.zeros(n, dtype=np.int8)
        cur = -2                      # -2: before outer list, -1: outer list opened
        pending: list[float] = []
        for op, arg, _ in pickletools.genops(f):
            name = op.name
            if name == "EMPTY_LIST":
                cur += 1
                pending = []
            elif name == "BINFLOAT" and cur >= 0:
                pending.append(arg)
                if len(pending) == 2:
                    k = count[cur]
                    xs[cur, k], ys[cur, k] = pending
                    count[cur] = k + 1
                    pending = []
            elif name == "STOP":
                break
        if cur != n - 1:
            raise RuntimeError(f"parsed {cur + 1} configurations, expected {n}")
    with zf.open(z_member) as f:
        z = np.load(f).astype(np.int16)
    pos = np.full((n, 16), -1, dtype=np.int16)
    valid = ~np.isnan(xs)
    pos[valid] = position_index(xs[valid], ys[valid])
    np.savez_compressed(dst, pos=pos, n_leds=count, z=z)
    print(f"{n} configurations written to {dst}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
