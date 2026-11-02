"""
Mean irradiance and homogeneity of all methods, reflector effects and uncertainties.

This module computes the numbers of the validation and reflector sections of the manuscript and of
Supporting Information S10 that are not part of a figure:

* :func:`metrics_all_methods`: E_mean and H of the geometric model, ray tracing and radiometry for
  all investigated configurations with and without reflectors, with the conventions of Figures 3,
  7, 8 and S17 (deviation of the homogeneity in the Model Validation section).
* :func:`reflector_effects` and :func:`reflector_spread`: gain of the mean irradiance by the
  reflectors, the ratio between the highest and the lowest mean irradiance of UV1-UV3 / GRN1-GRN3
  and their range of homogeneity, with and without reflectors.
* :func:`radiometry_corrections`: Table S4, the corrections of the radiometric mean irradiance for
  the share of the emitted flux that reaches the scanned area and for the 1 cm safety offset of the
  scans, both estimated with the geometric model.
* :func:`ray_tracing_noise`: share of the rays that reach the detection plane, the resulting
  statistical uncertainty of the mean irradiance, and the effect of the Monte-Carlo noise on the
  homogeneity estimated from the ray independence study of Supporting Information S3.

**Usage Example**::

    from analysis.uncertainty import metrics_all_methods, radiometry_corrections
    metrics = metrics_all_methods()
    table_s4 = radiometry_corrections()
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

from analysis.compare_methods import read_configuration
from analysis.fast_evaluation import GRID_KW, LEDS, PREFIX, RESOLUTION, configuration_file
from analysis.process_radiometry_data import calculate_correction_factor, process_radiometry_file, read_radiometry_file
from analysis.process_raytracing_data import calc_metrics, process_raytracing_file
from geometric_model import GeometricModel, calculate_irradiance
from utils.grid import Grid

#: Number of rays of every ray-tracing simulation.
N_RAYS = 1.6e7
#: Area of the detection plane in m2.
PLANE_AREA_M2 = 0.33 * 0.34
#: Width and height of the detection plane in cm, as passed to the radiometry processing.
PLANE_HEIGHT_CM, PLANE_WIDTH_CM = 33, 34
#: Safety offset of the radiometric scans in cm.
SCAN_OFFSET_CM = 1.0
#: Height of the reference scans of UV1-UV3 / GRN1-GRN3 (UV3 / GRN3, nominal 6 cm) in cm.
REFERENCE_SCAN_DISTANCE_CM = 7.0
#: Pitch of the fine grid on which the flux on the scanned area is integrated, in cm.
FINE_PITCH_CM = 0.05
#: Ray counts of the ray independence study compared in Supporting Information S10.
RAYS_COMPARED = (6.4e7, 1.2e8)
#: Blocks of ``results/raytracing/ray_independence_test.xlsx``: (LED type, reflector) -> rows, columns.
RAY_INDEPENDENCE_BLOCKS = {("UV", True): (slice(3, 10), slice(0, 3)), ("UV", False): (slice(3, 10), slice(3, 6)),
                           ("GRN", True): (slice(15, 22), slice(0, 3)), ("GRN", False): (slice(15, 22), slice(3, 6))}


def configuration_names(key: str) -> list[str]:
    """Names of the investigated configurations of one LED type: best manual, X1, X2, X3."""
    return ["best manual"] + [f"{PREFIX[key]}{k}" for k in (1, 2, 3)]


def _stem(name: str) -> str:
    return "best_manual" if name == "best manual" else name


def _reference(key: str, name: str) -> tuple[str, int]:
    """Radiometric reference scan and its number of LEDs, as in the manuscript."""
    if name == "best manual":
        return "best_manual_no_reflector", 8
    return f"{PREFIX[key]}3_reflector", 16


def _radiometry_file(key: str, stem: str, results_path: str) -> str:
    return os.path.join(results_path, "radiometry", LEDS[key].name, f"{stem}.csv")


def _model_field(X: np.ndarray, Y: np.ndarray, positions: list, z: float, led) -> np.ndarray:
    """Irradiance of the geometric model in W m-2 on an arbitrary grid (positions and z in cm)."""
    led.calc_radiant_intensity()
    X = np.ascontiguousarray(X, dtype=np.float64)
    Y = np.ascontiguousarray(Y, dtype=np.float64)
    return sum(calculate_irradiance(X - p[0], Y - p[1], float(z), led.angle_data, led.intensity_data)
               for p in positions) * 1e4


def _captured_fraction(scan_file: str, positions: list, z: float, led) -> float:
    """Share of the emitted flux that the geometric model puts onto the square covered by a scan."""
    scan = read_radiometry_file(scan_file)
    c = scan.columns.to_numpy(float)
    side = c.max() - c.min() + (c[1] - c[0])
    fine = np.arange(-side / 2 + FINE_PITCH_CM / 2, side / 2, FINE_PITCH_CM)
    Xf, Yf = np.meshgrid(fine, fine)
    E = _model_field(Xf, Yf, positions, z, led)
    return float(E.sum() * (FINE_PITCH_CM / 100) ** 2 / (len(positions) * led.total_power))


def metrics_all_methods(results_path: str = "results/") -> pd.DataFrame:
    """
    E_mean and H of the geometric model, ray tracing and radiometry for all investigated configurations.

    The radiometric scans are scaled as in the manuscript: the best manual configuration with its
    own reflector-free scan, UV1-UV3 / GRN1-GRN3 with the scan of UV3 / GRN3 with reflectors, each
    set to N x P_LED. The geometric model is given for the reflector-free cases only.

    :param results_path: Path of the results directory.
    :type results_path: str, optional
    :returns: One row per LED type, configuration and reflector setting, including the share of
        the rays that reach the detection plane and the resulting relative standard uncertainty of
        the ray-tracing mean irradiance, 1 / sqrt(n).
    :rtype: pandas.DataFrame
    """
    rows = []
    for key, led in LEDS.items():
        model = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        for name in configuration_names(key):
            positions, z = read_configuration(configuration_file(key, name, results_path))
            Eg = model.simulate(positions, z)
            ref_stem, n_ref = _reference(key, name)
            cf = calculate_correction_factor(_radiometry_file(key, ref_stem, results_path), n_ref * led.total_power)
            for reflector in (False, True):
                stem = f"{_stem(name)}_{'reflector' if reflector else 'no_reflector'}"
                rt = process_raytracing_file(f"{stem}.txt", os.path.join(results_path, "raytracing", led.name))
                rt_E, rt_H = calc_metrics(rt["Irradiance"])
                _, rad_E, rad_H = process_radiometry_file(_radiometry_file(key, stem, results_path), cf,
                                                          PLANE_HEIGHT_CM, PLANE_WIDTH_CM)
                share = rt_E * PLANE_AREA_M2 / (len(positions) * led.total_power)
                rows.append(dict(led=key, config=name, z_cm=z, n_leds=len(positions), reflector=reflector,
                                 model_E=np.nan if reflector else Eg.mean(),
                                 model_H=np.nan if reflector else 1 - Eg.std() / Eg.mean(),
                                 rt_E=rt_E, rt_H=rt_H, rad_E=rad_E, rad_H=rad_H, rad_reference=ref_stem,
                                 rt_rays_on_plane_share=share,
                                 rt_E_mean_stat_rel_pct=100 / np.sqrt(N_RAYS * share)))
    return pd.DataFrame(rows)


def reflector_effects(metrics: pd.DataFrame) -> pd.DataFrame:
    """
    Reflector gains and the deviations between the methods, per configuration.

    :param metrics: Result of :func:`metrics_all_methods`.
    :type metrics: pandas.DataFrame
    :returns: Gain of the mean irradiance by the reflectors (radiometry and ray tracing), the
        homogeneities with and without reflectors, and the relative (%) and absolute deviations of
        the geometric model and of the ray tracing from the radiometry.
    :rtype: pandas.DataFrame
    """
    rows = []
    for (key, name), g in metrics.groupby(["led", "config"], sort=False):
        a, b = g[~g.reflector].iloc[0], g[g.reflector].iloc[0]
        rows.append(dict(led=key, config=name, z_cm=a.z_cm,
                         gain_rad_pct=100 * (b.rad_E / a.rad_E - 1), gain_rt_pct=100 * (b.rt_E / a.rt_E - 1),
                         E_rad_no=a.rad_E, E_rad_refl=b.rad_E, H_rad_no=a.rad_H, H_rad_refl=b.rad_H,
                         H_rt_no=a.rt_H, H_rt_refl=b.rt_H,
                         model_vs_rad_E_rel_pct=100 * (a.model_E / a.rad_E - 1),
                         model_vs_rad_H_rel_pct=100 * (a.model_H / a.rad_H - 1),
                         model_vs_rad_H_abs=a.model_H - a.rad_H,
                         rt_vs_rad_E_rel_pct=100 * (a.rt_E / a.rad_E - 1), rt_vs_rad_H_abs=a.rt_H - a.rad_H,
                         rt_vs_rad_reflector_E_rel_pct=100 * (b.rt_E / b.rad_E - 1),
                         rt_vs_rad_reflector_H_abs=b.rt_H - b.rad_H))
    return pd.DataFrame(rows)


def reflector_spread(metrics: pd.DataFrame) -> pd.DataFrame:
    """
    Spread of the measured mean irradiance and homogeneity between UV1-UV3 / GRN1-GRN3.

    :param metrics: Result of :func:`metrics_all_methods`.
    :type metrics: pandas.DataFrame
    :returns: Per LED type and reflector setting: the measured E_mean and H of the three
        configurations, the ratio of the highest to the lowest E_mean and the range of H.
    :rtype: pandas.DataFrame
    """
    rows = []
    for key in LEDS:
        for reflector in (False, True):
            g = metrics[(metrics.led == key) & (metrics.config != "best manual") & (metrics.reflector == reflector)]
            rows.append(dict(led=key, reflector=reflector,
                             E_rad=" / ".join(f"{v:.1f}" for v in g.rad_E), H_rad=" / ".join(f"{v:.3f}" for v in g.rad_H),
                             E_max_over_min=g.rad_E.max() / g.rad_E.min(), H_range=g.rad_H.max() - g.rad_H.min()))
    return pd.DataFrame(rows)


def radiometry_corrections(results_path: str = "results/") -> pd.DataFrame:
    """
    Compute Table S4: corrections of the radiometric mean irradiance of the reflector-free configurations.

    *Detected share.* A reference scan is scaled such that its integral equals N x P_LED, i.e. it is
    assumed that the entire emitted flux reaches the scanned area. For the best manual configuration
    (its own reflector-free reference, scanned at 14 cm) the share that actually reaches the scanned
    square is integrated from the geometric model at the scan distance. UV1-UV3 / GRN1-GRN3 were
    scaled with the scan of UV3 / GRN3 with reflectors, which the model does not describe; its share
    is the model share of the reflector-free scan of the same configuration at the same distance
    (7 cm) times the ratio of the integrated signals of the two scans.

    *Scan distance.* The scans were recorded 1 cm further from the LEDs than the nominal height. The
    change from z to z + 1 cm is taken from the geometric model.

    The corrected radiometric mean irradiance is the published one times the detected share times
    E_model(z) / E_model(z + 1 cm).

    :param results_path: Path of the results directory.
    :type results_path: str, optional
    :returns: One row per reflector-free configuration.
    :rtype: pandas.DataFrame
    """
    rows = []
    for key, led in LEDS.items():
        model = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        x3 = f"{PREFIX[key]}3"
        positions3, _ = read_configuration(configuration_file(key, x3, results_path))
        scan_no_reflector = _radiometry_file(key, f"{x3}_no_reflector", results_path)
        scan_reflector = _radiometry_file(key, f"{x3}_reflector", results_path)
        model_share_x3 = _captured_fraction(scan_no_reflector, positions3, REFERENCE_SCAN_DISTANCE_CM, led)
        # integrated signal of the reflector scan relative to the reflector-free scan (both in raw units)
        expected = 16 * led.total_power
        signal_ratio = (calculate_correction_factor(scan_no_reflector, expected)
                        / calculate_correction_factor(scan_reflector, expected))
        for name in configuration_names(key):
            positions, z = read_configuration(configuration_file(key, name, results_path))
            ref_stem, n_ref = _reference(key, name)
            if name == "best manual":
                share = _captured_fraction(_radiometry_file(key, ref_stem, results_path), positions,
                                           z + SCAN_OFFSET_CM, led)
                model_share, ratio = share, np.nan
            else:
                share, model_share, ratio = model_share_x3 * signal_ratio, model_share_x3, signal_ratio
            cf = calculate_correction_factor(_radiometry_file(key, ref_stem, results_path), n_ref * led.total_power)
            _, rad_E, _ = process_radiometry_file(_radiometry_file(key, f"{_stem(name)}_no_reflector", results_path),
                                                  cf, PLANE_HEIGHT_CM, PLANE_WIDTH_CM)
            E0, E1 = model.simulate(positions, z), model.simulate(positions, z + SCAN_OFFSET_CM)
            m0, m1 = E0.mean(), E1.mean()
            distance_factor = m0 / m1
            rad_E_corrected = rad_E * share * distance_factor
            rows.append(dict(led=key, config=name, z_cm=z, reference=ref_stem,
                             model_share_of_reference_scan=model_share, signal_ratio_reflector_to_no_reflector=ratio,
                             detected_share=share, implied_overestimate_pct=100 * (1 / share - 1),
                             dE_mean_scan_distance_pct=100 * (m1 / m0 - 1),
                             dH_scan_distance=(1 - E1.std() / m1) - (1 - E0.std() / m0),
                             distance_factor=distance_factor, net_factor=share * distance_factor,
                             rad_E_published=rad_E, rad_E_corrected=rad_E_corrected, model_E=m0,
                             model_vs_rad_published_pct=100 * (m0 / rad_E - 1),
                             model_vs_rad_corrected_pct=100 * (m0 / rad_E_corrected - 1)))
    return pd.DataFrame(rows)


def ray_tracing_noise(metrics: pd.DataFrame, results_path: str = "results/") -> pd.DataFrame:
    """
    Effect of the Monte-Carlo noise of the ray tracing on the homogeneity (Supporting Information S10).

    The ray independence study (Supporting Information S3, Figure S5) gives the standard deviation
    sigma of the irradiance data set of the best manual configuration for 2e6 ... 1.2e8 rays. sigma
    at 1.2e8 rays serves as a rough estimate of the noise-free standard deviation; combined with the
    published ray-tracing mean irradiance it gives the homogeneity without most of the noise. The
    difference to the published ray-tracing homogeneity (1.6e7 rays) is a lower estimate of the
    effect of the noise.

    :param metrics: Result of :func:`metrics_all_methods`.
    :type metrics: pandas.DataFrame
    :param results_path: Path of the results directory.
    :type results_path: str, optional
    :returns: One row per LED type and reflector setting.
    :rtype: pandas.DataFrame
    """
    raw = pd.read_excel(os.path.join(results_path, "raytracing", "ray_independence_test.xlsx"),
                        sheet_name="Planilha1", header=None)
    rows = []
    for (key, reflector), (r, c) in RAY_INDEPENDENCE_BLOCKS.items():
        block = raw.iloc[r, c].to_numpy(dtype=float)
        rays, sigma = block[:, 0], block[:, 2]
        s_lo, s_hi = (float(sigma[np.argmin(np.abs(rays - n))]) for n in RAYS_COMPARED)
        bm = metrics[(metrics.led == key) & (metrics.config == "best manual") & (metrics.reflector == reflector)].iloc[0]
        h_hi = 1 - s_hi / bm.rt_E
        rows.append(dict(led=key, reflector=reflector, sigma_6_4e7=s_lo, sigma_1_2e8=s_hi,
                         sigma_change_pct=100 * (s_hi / s_lo - 1), rt_E_published=bm.rt_E, rt_H_published=bm.rt_H,
                         H_from_sigma_1_2e8=h_hi, dH_noise=h_hi - bm.rt_H))
    return pd.DataFrame(rows)
