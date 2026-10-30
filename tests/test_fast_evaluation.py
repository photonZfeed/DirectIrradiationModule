"""
The fast evaluation of :mod:`analysis.fast_evaluation` equals the geometric model of the repository.

The mean, standard deviation and maximum irradiance of all 482,931 configurations of the systematic
8-LED study at 13 cm are compared with the values stored by :mod:`geometric_model` in
``results/geometric/<LED>/simulation_results_systematic_8_leds_13_cm_482931_samples.npz``.
The maximum is evaluated in single precision, hence the tolerance of 1e-6.

Run from the repository root with ``python -m pytest tests`` or ``python tests/test_fast_evaluation.py``.
"""
import os
import sys

if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402

from analysis.fast_evaluation import LEDS, FieldBank, configs_from_object_array  # noqa: E402

STORED = "results/geometric/{}/simulation_results_systematic_8_leds_13_cm_482931_samples.npz"
TOLERANCE = 1e-6


def max_relative_deviation(key: str) -> dict[str, float]:
    """Largest relative deviation of mean, std and max from the stored results for one LED type."""
    d = np.load(STORED.format(LEDS[key].name), allow_pickle=True)
    pos = configs_from_object_array(d["configs"])
    z = np.asarray(d["z"]).astype(int)
    mean, std, mx = FieldBank(LEDS[key], heights=tuple(np.unique(z))).evaluate(pos, z)
    ref = {k: np.asarray(d[k], dtype=float) for k in ("means", "stds", "max")}
    return {name: float(np.max(np.abs(value - ref[k]) / np.abs(ref[k])))
            for name, value, k in (("mean", mean, "means"), ("std", std, "stds"), ("max", mx, "max"))}


def test_fast_evaluation_uv():
    assert max(max_relative_deviation("UV").values()) < TOLERANCE


def test_fast_evaluation_green():
    assert max(max_relative_deviation("GRN").values()) < TOLERANCE


if __name__ == "__main__":
    for key in LEDS:
        print(key, max_relative_deviation(key))
