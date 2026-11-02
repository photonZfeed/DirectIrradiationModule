"""
Fast and exact evaluation of many LED configurations with the geometric model.

The Monte-Carlo studies of the Supporting Information (Sections S11 and S12) evaluate up to
10 million configurations per random seed. Simulating each configuration with
:meth:`geometric_model.GeometricModel.simulate` would take days, so this module uses the linearity
of the model instead.

For a fixed height ``z`` the irradiance field of a configuration is the superposition of the
fields of its LEDs (Eq. (10) of the manuscript), ``E = m^T F_z``, where ``m`` is the 0/1 indicator
vector over the 169 mounting positions and ``F_z`` (169 x number of grid points) holds the
single-LED fields computed with :meth:`GeometricModel.simulate`. Hence

* ``E_mean = m^T mu_z`` with ``mu_z`` the mean of every single-LED field,
* ``sum(E^2) = m^T G_z m`` with the Gram matrix ``G_z = F_z F_z^T``,
* ``std(E) = sqrt(sum(E^2) / n - E_mean^2)`` (population standard deviation, as :func:`numpy.std`),
* ``max(E)`` needs the field itself; it is evaluated by a matrix product on one half of the grid.
  This is exact because every sampled configuration is point-symmetric about the module centre and
  the evaluation grid is symmetric about it with an even number of rows (both are checked).

``tests/test_fast_evaluation.py`` checks that the result equals the stored results of the
systematic 8-LED study (relative deviations below 1e-6; the maximum is evaluated in single
precision).

The module also holds the selection procedure of the manuscript on precomputed metrics
(:func:`select_configurations`) and the non-dominated filter used for the Pareto fronts.
"""
from __future__ import annotations

import os

import numpy as np

from geometric_model import GeometricModel
from utils.grid import Grid
from utils.led import LST1_01F06_GRN1_00, LST1_01G01_UV01_00

#: Grid of the direct irradiation module (33 cm x 34 cm detection plane, 13 x 13 mounting
#: positions with 2.5 cm pitch), as in all figures of the manuscript.
GRID_KW = dict(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.0)
#: Resolution of the evaluation grid of the geometric model in cm.
RESOLUTION = 0.5
#: Heights of the Monte-Carlo optimization in cm (5-15 cm in 1 cm steps).
HEIGHTS = tuple(range(5, 16))
#: Pitch of the mounting grid in cm.
STEP = 2.5
#: Half width of the mounting grid in pitches (positions run from -6 to +6).
HALF = 6
#: Number of mounting positions of the 13 x 13 mounting grid.
N_POSITIONS = (2 * HALF + 1) ** 2
#: Irradiance bound of the manuscript in W m-2.
E_BOUND = 150.0
#: Weighting factors t of the Pareto analysis of the manuscript (Eq. (5)).
T_VALUES = np.linspace(0.0, 1.0, 100)
#: LED objects by short key.
LEDS = {"UV": LST1_01G01_UV01_00, "GRN": LST1_01F06_GRN1_00}
#: Prefix of the configuration names (UV1-UV3, GRN1-GRN3) by short key.
PREFIX = {"UV": "UV", "GRN": "GRN"}


def configuration_file(key: str, name: str, results_path: str = "results/") -> str:
    """
    Path of the JSON file of an investigated configuration.

    :param key: LED type, ``"UV"`` or ``"GRN"``.
    :type key: str
    :param name: ``"best manual"`` or one of UV1-UV3 / GRN1-GRN3.
    :type name: str
    :param results_path: Path of the results directory.
    :type results_path: str, optional
    :returns: Path of the configuration file.
    :rtype: str
    """
    if name == "best manual":
        return os.path.join(results_path, "sampled_configs", "best_manual.json")
    return os.path.join(results_path, "sampled_configs", "near_optimal_configs", LEDS[key].name, f"{name}.json")


def position_index(x_cm: np.ndarray, y_cm: np.ndarray) -> np.ndarray:
    """
    Map mounting coordinates to the index 0 ... 168 of the 13 x 13 mounting grid.

    :param x_cm: x coordinates of the LEDs in cm (multiples of 2.5 cm).
    :type x_cm: numpy.ndarray
    :param y_cm: y coordinates of the LEDs in cm (multiples of 2.5 cm).
    :type y_cm: numpy.ndarray
    :returns: Position indices, row-major in x, then y.
    :rtype: numpy.ndarray
    :raises ValueError: If a position lies outside the mounting grid.
    """
    ix = np.rint(np.asarray(x_cm, dtype=float) / STEP).astype(int) + HALF
    iy = np.rint(np.asarray(y_cm, dtype=float) / STEP).astype(int) + HALF
    if np.any((ix < 0) | (ix > 2 * HALF) | (iy < 0) | (iy > 2 * HALF)):
        raise ValueError("position outside the 13 x 13 mounting grid")
    return ix * (2 * HALF + 1) + iy


def index_to_xy(idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Map position indices 0 ... 168 back to mounting coordinates in cm.

    :param idx: Position indices.
    :type idx: numpy.ndarray
    :returns: x and y coordinates in cm.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    idx = np.asarray(idx)
    return (idx // (2 * HALF + 1) - HALF) * STEP, (idx % (2 * HALF + 1) - HALF) * STEP


def configs_from_object_array(configs) -> np.ndarray:
    """
    Convert configurations stored as lists of (x, y) positions in cm to padded index arrays.

    :param configs: Sequence of configurations, each a sequence of (x, y) tuples in cm (as stored in
        the ``.npz`` files of the samplers and the geometric model).
    :returns: Integer array of shape (n, 16) with the position indices, padded with -1.
    :rtype: numpy.ndarray
    """
    out = np.full((len(configs), 16), -1, dtype=np.int16)
    for i, c in enumerate(configs):
        c = np.asarray(c, dtype=float)
        out[i, : len(c)] = position_index(c[:, 0], c[:, 1])
    return out


class FieldBank:
    """
    Single-LED fields, their means and Gram matrices for one LED type and a set of heights.

    :param led: LED object, e.g. :data:`utils.led.LST1_01G01_UV01_00`.
    :type led: utils.led.LED
    :param heights: Heights in cm for which the fields are precomputed.
    :type heights: tuple[int, ...], optional
    """

    def __init__(self, led, heights=HEIGHTS):
        self.led = led
        self.model = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        ny, nx = self.model.X.shape
        if ny % 2:
            raise ValueError("the half-grid evaluation of the maximum needs an even number of grid rows")
        self.shape = (ny, nx)
        self.n_grid = nx * ny
        self.half_rows = ny // 2              # rows 0 ... ny/2 - 1 form one half of the symmetric grid
        xs, ys = index_to_xy(np.arange(N_POSITIONS))
        self.fields, self.mu, self.gram, self.half = {}, {}, {}, {}
        for h in heights:
            F = np.stack([self.model.simulate([(x, y)], float(h)).ravel() for x, y in zip(xs, ys)])
            self.fields[h] = F
            self.mu[h] = F.mean(axis=1)
            self.gram[h] = F @ F.T
            self.half[h] = np.ascontiguousarray(
                F.reshape(-1, ny, nx)[:, : self.half_rows, :].reshape(len(xs), -1).astype(np.float32))

    def evaluate(self, pos: np.ndarray, z: np.ndarray, with_max: bool = True, chunk: int = 20000):
        """
        Mean, standard deviation and maximum of the irradiance of many configurations.

        :param pos: Position indices of shape (n, 16), padded with -1. Every configuration must be
            point-symmetric about the module centre if ``with_max`` is True.
        :type pos: numpy.ndarray
        :param z: Heights in cm, shape (n,).
        :type z: numpy.ndarray
        :param with_max: Also compute the maximum local irradiance.
        :type with_max: bool, optional
        :param chunk: Number of configurations evaluated at once (memory use).
        :type chunk: int, optional
        :returns: Mean, standard deviation and maximum in W m-2 (the maximum is NaN if
            ``with_max`` is False).
        :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
        :raises ValueError: If ``with_max`` is True and a configuration is not point-symmetric.
        """
        pos = np.asarray(pos)
        z = np.asarray(z).astype(int)
        valid = pos >= 0
        if with_max:
            mirrored = np.where(valid, N_POSITIONS - 1 - pos, -1)
            if not np.array_equal(np.sort(pos, axis=1), np.sort(mirrored, axis=1)):
                raise ValueError("the half-grid maximum is only exact for point-symmetric configurations")
        n = len(pos)
        mean = np.empty(n)
        std = np.empty(n)
        mx = np.full(n, np.nan)
        for h in np.unique(z):
            sel = np.flatnonzero(z == h)
            mu, G = self.mu[int(h)], self.gram[int(h)]
            for s in range(0, len(sel), chunk):
                rows = sel[s:s + chunk]
                p, v = pos[rows], valid[rows]
                pp = np.where(v, p, 0)
                mean[rows] = np.where(v, mu[pp], 0.0).sum(axis=1)
                gsum = np.where(v[:, :, None] & v[:, None, :], G[pp[:, :, None], pp[:, None, :]], 0.0)
                e2 = gsum.sum(axis=(1, 2)) / self.n_grid
                std[rows] = np.sqrt(np.maximum(e2 - mean[rows] ** 2, 0.0))
                if with_max:
                    M = np.zeros((len(rows), self.half[int(h)].shape[0]), dtype=np.float32)
                    r_idx = np.repeat(np.arange(len(rows)), v.sum(axis=1))
                    M[r_idx, p[v]] = 1.0
                    mx[rows] = (M @ self.half[int(h)]).max(axis=1)
        return mean, std, mx


def pareto_mask(e: np.ndarray, h: np.ndarray) -> np.ndarray:
    """
    Non-dominated set of points when both coordinates are maximised.

    :param e: First objective, e.g. E_mean / E_bound.
    :type e: numpy.ndarray
    :param h: Second objective, e.g. the homogeneity H.
    :type h: numpy.ndarray
    :returns: Boolean mask of the non-dominated points (of identical points, one is kept).
    :rtype: numpy.ndarray
    """
    order = np.lexsort((-h, -e))          # e descending, then h descending
    hs = h[order]
    prev = np.concatenate(([-np.inf], np.maximum.accumulate(hs)[:-1]))
    mask = np.zeros(len(e), dtype=bool)
    mask[order[hs > prev]] = True
    return mask


def nondominated_front(mean: np.ndarray, std: np.ndarray, mx: np.ndarray,
                       bound: float = E_BOUND) -> tuple[np.ndarray, np.ndarray]:
    """
    Non-dominated feasible configurations in the plane (E_mean / E_bound, H).

    Configurations whose maximum local irradiance exceeds ``bound`` are excluded, as in
    :class:`optimizer.Optimizer`.

    :returns: E_mean / E_bound and H of the front, sorted by increasing E_mean / E_bound.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    feas = mx <= bound
    e, h = mean[feas] / bound, 1.0 - std[feas] / mean[feas]
    nd = pareto_mask(e, h)
    order = np.argsort(e[nd])
    return e[nd][order], h[nd][order]


def select_configurations(mean: np.ndarray, std: np.ndarray, mx: np.ndarray,
                          bound: float = E_BOUND, t_values: np.ndarray = T_VALUES):
    """
    Selection procedure of the manuscript on precomputed metrics.

    Reproduces :meth:`optimizer.Optimizer.create_pareto_front`,
    :meth:`optimizer.Optimizer.find_best_configuration` and the export of the configurations for
    ray tracing in :func:`visualization.pareto_plots_H_vs_E_ratio.process_led`: every weighting
    factor t contributes the configuration that maximises ``t E_mean/E_bound + (1 - t) H``
    (Eq. (5)); configurations whose local irradiance exceeds the bound score zero. The best
    configuration maximises ``H E_mean / E_bound``, and the configurations exported for ray tracing
    are the best one and the two front points preceding it.

    :param bound: Irradiance bound in W m-2; ``numpy.inf`` disables the bound (the mean
        irradiance is then normalised by :data:`E_BOUND`).
    :returns: Index of the best configuration, indices of the exported configurations (ordered by
        increasing E_mean), the normalised mean irradiance and the homogeneity (both zero for
        infeasible configurations).
    :rtype: tuple[int, list[int], numpy.ndarray, numpy.ndarray]
    """
    feas = mx <= bound
    norm = bound if np.isfinite(bound) else E_BOUND
    ratio = np.where(feas, mean / norm, 0.0)
    hf = np.where(feas, 1.0 - std / mean, 0.0)
    front = []
    for t in t_values:
        i = int(np.argmax(t * ratio + (1 - t) * hf))
        if i not in front:
            front.append(i)
    best = int(np.argmax(ratio * hf))
    if best not in front:
        # the product optimum is not reached by any of the t values; insert it at its place
        # along the front, which is ordered by increasing E_mean
        front = sorted(front + [best], key=lambda i: ratio[i])
    b = front.index(best)
    return best, front[max(0, b - 2): b + 1], ratio, hf
