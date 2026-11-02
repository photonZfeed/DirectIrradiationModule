"""
Sensitivity of the geometric model to the point-source assumption (Supporting Information S9).

The geometric model treats every LED as a point source (Eq. (9) of the manuscript). Here every
LED is replaced by a uniformly emitting planar disk of diameter ``d`` whose surface elements all
radiate with the radiant-intensity profile of Eq. (6), weighted by their share of the disk area.
The contribution of an LED is therefore the point-source contribution averaged over the disk
(Eq. (S10) of the Supporting Information)::

    E_i,disk(x, y) = 4 / (pi d^2) * int_0^{2 pi} int_0^{d/2} E_i(x - rho cos phi, y - rho sin phi) rho drho dphi

It is evaluated with 12 Gauss-Legendre nodes in (2 rho / d)^2 (area-uniform) and 24 equidistant
nodes in phi. ``d`` is the lens diameter given in the datasheets: 2.6 mm for the 365 nm LED
(Luminus SST-10-UV, lens radius R1.30) and 2.15 mm for the 530 nm LED (OSRAM GT CSSPM1.13).

Evaluated are all investigated configurations at their heights (best manual configuration for both
LED types, UV1-UV3, GRN1-GRN3) and UV3 and GRN3 at 5 cm, the smallest height of the optimization.
The result is Table S3 of the Supporting Information.

**Usage Example**::

    from analysis.source_model import point_source_table
    table = point_source_table()
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from analysis.compare_methods import read_configuration
from analysis.fast_evaluation import GRID_KW, LEDS, RESOLUTION, configuration_file
from geometric_model import GeometricModel
from utils.grid import Grid

#: Lens diameters from the datasheets in mm (also the source diameters of the ray tracing).
LENS_DIAMETER_MM = {"UV": 2.6, "GRN": 2.15}
#: Configurations of Table S3: (LED type, name, height override in cm or None for the design height).
CONFIGURATIONS = [(key, name, z) for key in ("UV", "GRN")
                  for name, z in [("best manual", None)] + [(f"{key}{k}", None) for k in (1, 2, 3)] + [(f"{key}3", 5.0)]]


def disk_nodes(d_cm: float, n_r: int = 12, n_phi: int = 24) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Quadrature nodes on a disk: Gauss-Legendre in the area coordinate, uniform in the azimuth.

    :param d_cm: Disk diameter in cm; 0 returns the centre only (point source).
    :type d_cm: float
    :param n_r: Number of radial nodes.
    :type n_r: int, optional
    :param n_phi: Number of azimuthal nodes.
    :type n_phi: int, optional
    :returns: x and y offsets of the nodes in cm and their weights (summing to 1).
    :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
    """
    if d_cm == 0:
        return np.zeros(1), np.zeros(1), np.ones(1)
    a = d_cm / 2
    xg, wg = np.polynomial.legendre.leggauss(n_r)
    r = a * np.sqrt((xg + 1) / 2)                   # area-uniform mapping
    wr = wg / 2
    phi = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    R, P = np.meshgrid(r, phi, indexing="ij")
    W = np.repeat(wr[:, None], n_phi, axis=1) / n_phi
    return (R * np.cos(P)).ravel(), (R * np.sin(P)).ravel(), W.ravel()


def disk_field(X: np.ndarray, Y: np.ndarray, positions: list, z: float, angle_deg: np.ndarray,
               intensity: np.ndarray, d_cm: float = 0.0) -> np.ndarray:
    """
    Irradiance of LEDs modelled as uniformly emitting disks (Eq. (S10)).

    With ``d_cm = 0`` the result equals :meth:`geometric_model.GeometricModel.simulate`.

    :param X: x coordinates of the evaluation grid in cm.
    :param Y: y coordinates of the evaluation grid in cm.
    :param positions: LED positions (x, y) in cm.
    :param z: Height of the LEDs above the detection plane in cm.
    :param angle_deg: Angles of the radiant-intensity profile in degrees.
    :param intensity: Radiant intensity in W sr-1 (after :meth:`utils.led.LED.calc_radiant_intensity`).
    :param d_cm: Disk diameter in cm.
    :returns: Irradiance in W m-2.
    :rtype: numpy.ndarray
    """
    ex, ey, ew = disk_nodes(d_cm)
    E = np.zeros_like(X)
    for (px, py) in positions:
        for dx, dy, w in zip(ex, ey, ew):
            xs, ys = X - (px + dx), Y - (py + dy)
            rho2 = xs ** 2 + ys ** 2
            theta = np.arctan(np.sqrt(rho2) / z)
            radiant_intensity = np.interp(np.degrees(theta), angle_deg, intensity)
            E += w * z / (rho2 + z ** 2) ** 1.5 * radiant_intensity
    return E * 1e4


def point_source_table(results_path: str = "results/") -> pd.DataFrame:
    """
    Compute Table S3: point source versus disk with the datasheet lens diameter.

    :param results_path: Path of the results directory (configuration files).
    :type results_path: str, optional
    :returns: One row per configuration with E_mean and H of the point-source model, their changes
        for the disk source and the largest and root-mean-square change of the local irradiance in
        percent of the mean irradiance.
    :rtype: pandas.DataFrame
    :raises AssertionError: If the point-source evaluation differs from the repository model.
    """
    rows = []
    for key, name, z_override in CONFIGURATIONS:
        led = LEDS[key]
        model = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        positions, z = read_configuration(configuration_file(key, name, results_path))
        z = z_override or z
        E0 = disk_field(model.X, model.Y, positions, z, led.angle_data, led.intensity_data)
        if not np.allclose(E0, model.simulate(positions, z), rtol=1e-9, atol=1e-9):
            raise AssertionError("point-source baseline differs from the repository model")
        d_mm = LENS_DIAMETER_MM[key]
        E1 = disk_field(model.X, model.Y, positions, z, led.angle_data, led.intensity_data, d_cm=d_mm / 10)
        m0, m1 = E0.mean(), E1.mean()
        h0, h1 = 1 - E0.std() / m0, 1 - E1.std() / m1
        rows.append(dict(led=key, config=name, n_leds=len(positions), z_cm=z, d_mm=d_mm,
                         E_mean_point=m0, H_point=h0, E_mean_disk=m1, H_disk=h1,
                         dE_mean_rel_pct=100 * (m1 - m0) / m0, dH_abs=h1 - h0,
                         max_local_dev_pct=100 * np.max(np.abs(E1 - E0)) / m0,
                         rms_local_dev_pct=100 * np.sqrt(np.mean((E1 - E0) ** 2)) / m0))
    return pd.DataFrame(rows)
