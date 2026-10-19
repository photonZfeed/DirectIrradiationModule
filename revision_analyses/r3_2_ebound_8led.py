"""Reviewer 3, comment 2: is the irradiance bound active in the systematic 8-LED study at 13 cm?

Evaluates the highest local irradiance of all 482 931 symmetric 8-LED configurations at 13 cm
(``results/sampled_configs/systematic_8_leds_13_cm_482931_samples.npz`` of the repository) for both
LED types and writes the maximum over all configurations to ``r3_2_ebound_8led.csv``.

Usage: python revision_analyses/r3_2_ebound_8led.py [--repo <DirectIrradiationModule>] [--out <dir>]
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from r2_common import FieldBank, common_args, configs_from_object_array, repo_path

SAMPLES = "results/sampled_configs/systematic_8_leds_13_cm_482931_samples.npz"


def main() -> None:
    args = common_args(__doc__.splitlines()[0]).parse_args()
    repo = repo_path(args.repo)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    d = np.load(repo / SAMPLES, allow_pickle=True)
    height = int(d["height"])
    pos = configs_from_object_array(d["configurations"])
    z = np.full(len(pos), height, dtype=int)
    rows = []
    for led in ("UV", "GRN"):
        bank = FieldBank(repo, led, heights=(height,))
        mean, std, mx = bank.evaluate(pos, z, with_max=True)
        i = int(np.nanargmax(mx))
        rows.append(dict(led=led, n_configs=len(pos), n_leds=int(d["led_count"]), z_cm=height,
                         E_max_all_configs=float(mx[i]), E_mean_of_that_config=float(mean[i]),
                         n_above_150=int((mx > 150).sum())))
    df = pd.DataFrame(rows)
    df.to_csv(out / "r3_2_ebound_8led.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
