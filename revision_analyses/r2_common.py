"""Shared helpers for the analyses of the second review round (see README.md in this folder).

All analyses import the published code of the DirectIrradiationModule repository, so the
irradiance of a single LED is computed by exactly the function used for the manuscript
(``geometric_model.calculate_irradiance``). The repository is the one that contains this folder;
``--repo`` or the environment variable ``DIM_REPO`` override it. Results are written to
``results/revision_analyses`` of the repository unless ``--out`` is given; nothing is hard-coded.

Fast evaluation of many configurations
--------------------------------------
For a fixed height h the irradiance field of a configuration is the linear superposition of
single-LED fields, E = m^T F_h, where m is the 0/1 indicator vector over the 169 mounting
positions and F_h (169 x n_grid) holds the single-LED fields. Hence

* E_mean = m^T mu_h                with mu_h = mean of every single-LED field,
* sum(E^2) = m^T G_h m             with the Gram matrix G_h = F_h F_h^T,
* std(E) = sqrt(sum(E^2)/n - E_mean^2)  (population std, as numpy.std),
* max(E) needs the field itself; it is evaluated by a matrix product on half of the grid,
  which is exact because every sampled configuration is point-symmetric and the evaluation
  grid is symmetric about the module centre.

``validate`` reproduces the stored repository results (means, stds, max) to rounding error.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

GRID_KW = dict(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.0)
RESOLUTION = 0.5
HEIGHTS = tuple(range(5, 16))          # 5 ... 15 cm, 1 cm steps, as in the manuscript
STEP = 2.5                             # mounting-grid pitch in cm
HALF = 6                               # mounting grid runs from -6 to +6 steps
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "results" / "revision_analyses"


def repo_path(cli_value: str | None = None) -> Path:
    """Return the repository root: the CLI value, ``DIM_REPO`` or the repository of this folder."""
    value = cli_value or os.environ.get("DIM_REPO")
    path = Path(value).resolve() if value else REPO_ROOT
    if not (path / "geometric_model.py").exists():
        raise SystemExit(f"{path} is not the DirectIrradiationModule repository.")
    return path


def import_repo(repo: Path):
    """Make the repository importable and return (GeometricModel, Grid, LED, leds)."""
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from geometric_model import GeometricModel  # noqa: E402
    from utils.grid import Grid  # noqa: E402
    from utils.led import LED  # noqa: E402

    leds = {
        "UV": LED("LST1-01G01-UV01-00", 0.875,
                  str(repo / "utils/angle_intensity_data/LST1-01G01-UV01-00.csv")),
        "GRN": LED("LST1-01F06-GRN1-00", 0.281,
                   str(repo / "utils/angle_intensity_data/LST1-01F06-GRN1-00.csv")),
    }
    return GeometricModel, Grid, LED, leds


def position_index(x_cm: np.ndarray, y_cm: np.ndarray) -> np.ndarray:
    """Map mounting coordinates in cm to the index 0..168 (row-major in x, then y)."""
    ix = np.rint(np.asarray(x_cm) / STEP).astype(int) + HALF
    iy = np.rint(np.asarray(y_cm) / STEP).astype(int) + HALF
    if np.any((ix < 0) | (ix > 2 * HALF) | (iy < 0) | (iy > 2 * HALF)):
        raise ValueError("position outside the 13 x 13 mounting grid")
    return ix * (2 * HALF + 1) + iy


def index_to_xy(idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    idx = np.asarray(idx)
    return (idx // (2 * HALF + 1) - HALF) * STEP, (idx % (2 * HALF + 1) - HALF) * STEP


class FieldBank:
    """Single-LED fields, their means and Gram matrices for one LED type and a set of heights."""

    def __init__(self, repo: Path, led_key: str, heights=HEIGHTS):
        GeometricModel, Grid, _, leds = import_repo(repo)
        self.led = leds[led_key]
        self.model = GeometricModel(Grid(**GRID_KW), self.led, resolution_xy=RESOLUTION)
        ny, nx = self.model.X.shape
        self.shape = (ny, nx)
        self.n_grid = nx * ny
        self.half_rows = ny // 2              # rows 0..ny/2-1 form one half of a symmetric grid
        xs, ys = index_to_xy(np.arange((2 * HALF + 1) ** 2))
        self.fields, self.mu, self.gram, self.half = {}, {}, {}, {}
        for h in heights:
            F = np.stack([self.model.simulate([(x, y)], float(h)).ravel() for x, y in zip(xs, ys)])
            self.fields[h] = F
            self.mu[h] = F.mean(axis=1)
            self.gram[h] = F @ F.T
            self.half[h] = np.ascontiguousarray(
                F.reshape(-1, ny, nx)[:, : self.half_rows, :].reshape(len(xs), -1).astype(np.float32))

    def evaluate(self, pos: np.ndarray, z: np.ndarray, with_max: bool = True, chunk: int = 20000):
        """Mean, std and max irradiance (W m-2) for configurations.

        :param pos: int array (n, 16) of position indices, padded with -1.
        :param z: int array (n,) of heights in cm.
        """
        n = len(pos)
        mean = np.empty(n)
        std = np.empty(n)
        mx = np.full(n, np.nan)
        valid = pos >= 0
        for h in np.unique(z):
            sel = np.flatnonzero(z == h)
            mu, G = self.mu[int(h)], self.gram[int(h)]
            for s in range(0, len(sel), chunk):
                rows = sel[s:s + chunk]
                p, v = pos[rows], valid[rows]
                pp = np.where(v, p, 0)
                m = np.where(v, mu[pp], 0.0).sum(axis=1)
                gsum = np.where(v[:, :, None] & v[:, None, :], G[pp[:, :, None], pp[:, None, :]], 0.0)
                e2 = gsum.sum(axis=(1, 2)) / self.n_grid
                mean[rows] = m
                std[rows] = np.sqrt(np.maximum(e2 - mean[rows] ** 2, 0.0))
                if with_max:
                    M = np.zeros((len(rows), self.half[int(h)].shape[0]), dtype=np.float32)
                    r_idx = np.repeat(np.arange(len(rows)), v.sum(axis=1))
                    M[r_idx, p[v]] = 1.0
                    mx[rows] = (M @ self.half[int(h)]).max(axis=1)
        return mean, std, mx


def configs_from_object_array(configs) -> np.ndarray:
    """Convert the repository's object array of (x, y) lists (cm) to padded index arrays."""
    out = np.full((len(configs), 16), -1, dtype=np.int16)
    for i, c in enumerate(configs):
        c = np.asarray(c, dtype=float)
        out[i, : len(c)] = position_index(c[:, 0], c[:, 1])
    return out


def pareto_mask(e: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Non-dominated set (maximise both e and h). Ties in both are kept once."""
    order = np.lexsort((-h, -e))          # e descending, then h descending
    best_h = -np.inf
    mask = np.zeros(len(e), dtype=bool)
    for i in order:
        if h[i] > best_h:
            mask[i] = True
            best_h = h[i]
    return mask


def pareto_mask_fast(e: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Vectorised non-dominated filter for large arrays."""
    order = np.lexsort((-h, -e))
    hs = h[order]
    run = np.maximum.accumulate(hs)
    prev = np.concatenate(([-np.inf], run[:-1]))
    keep_sorted = hs > prev
    mask = np.zeros(len(e), dtype=bool)
    mask[order[keep_sorted]] = True
    return mask


def hypervolume_2d(e: np.ndarray, h: np.ndarray, ref=(0.0, 0.0)) -> float:
    """Area dominated by a set of points (maximisation) relative to ``ref``."""
    m = pareto_mask_fast(e, h)
    pe, ph = e[m], h[m]
    order = np.argsort(-pe)
    pe, ph = pe[order], ph[order]
    area, h_prev = 0.0, ref[1]
    for x, y in zip(pe, ph):
        if y > h_prev:
            area += (x - ref[0]) * (y - h_prev)
            h_prev = y
    return area


def common_args(description: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--repo", help="path to the DirectIrradiationModule repository "
                    "(default: DIM_REPO or the repository containing this folder)")
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="output folder (default: results/revision_analyses)")
    return ap
