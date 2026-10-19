"""Check that the fast evaluator reproduces the stored repository results exactly.

Compares mean, std and max irradiance of FieldBank.evaluate with
``results/geometric/<LED>/simulation_results_*.npz`` of the repository.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from r2_common import FieldBank, common_args, configs_from_object_array, repo_path


def compare(bank: FieldBank, npz: Path, limit: int | None = None) -> dict:
    d = np.load(npz, allow_pickle=True)
    configs = d["configs"][:limit]
    z = np.asarray(d["z"][:limit]).astype(int)
    pos = configs_from_object_array(configs)
    mean, std, mx = bank.evaluate(pos, z)
    ref_mean, ref_std, ref_max = (np.asarray(d[k][:limit], dtype=float) for k in ("means", "stds", "max"))
    rel = lambda a, b: float(np.max(np.abs(a - b) / np.abs(b)))
    return {
        "file": npz.name, "n": len(pos),
        "max_rel_dev_mean": rel(mean, ref_mean),
        "max_rel_dev_std": rel(std, ref_std),
        "max_rel_dev_max": rel(mx, ref_max),
    }


def main() -> None:
    ap = common_args(__doc__)
    ap.add_argument("--mc-uv", help="optional path of the 5 million-sample UV simulation file")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    repo = repo_path(args.repo)
    rows = []
    for key, led_name in (("UV", "LST1-01G01-UV01-00"), ("GRN", "LST1-01F06-GRN1-00")):
        bank = FieldBank(repo, key)
        rows.append(compare(bank, repo / f"results/geometric/{led_name}/"
                            "simulation_results_systematic_8_leds_13_cm_482931_samples.npz", args.limit))
        if key == "UV" and args.mc_uv:
            rows.append(compare(bank, Path(args.mc_uv), args.limit))
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
