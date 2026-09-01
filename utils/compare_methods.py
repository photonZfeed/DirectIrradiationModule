# written by claude cowork powered by Opus 5 (Anthropic, PBC)

import json
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import binned_statistic_2d

from geometric_model import GeometricModel
from utils.grid import Grid
from utils.led import LST1_01F06_GRN1_00, LST1_01G01_UV01_00
from utils.process_radiometry_data import calculate_correction_factor, process_radiometry_file
from utils.process_raytracing_data import process_raytracing_file

#: Width (x extent) of the detection plane in cm.
PLANE_WIDTH = 33.0
#: Height (y extent) of the detection plane in cm.
PLANE_HEIGHT = 34.0
#: Resolution of the common evaluation grid in cm.
RESOLUTION = 0.5


@dataclass(frozen=True)
class Case:
    """
    One configuration for which all three evaluation methods are available.

    :param label: Human-readable name of the configuration, e.g. ``"UV1"``.
    :type label: str
    :param stem: File stem used in ``results/raytracing/<led>/`` and ``results/radiometry/<led>/``.
    :type stem: str
    :param config_json: Path of the JSON file with the LED positions and the module height. ``None``
        selects ``results/sampled_configs/best_manual.json``.
    :type config_json: str | None
    :param n_leds: Number of LEDs of the configuration.
    :type n_leds: int
    :param reference_stem: Radiometric scan used to derive the photon-flux correction factor. The
        convention of the manuscript is reproduced: the best manual configuration serves as its own
        reference, the model-based configurations are referenced to UV3/GRN3 with reflectors.
    :type reference_stem: str
    :param reference_n_leds: Number of LEDs assumed for the reference measurement.
    :type reference_n_leds: int
    """

    label: str
    stem: str
    config_json: str | None
    n_leds: int
    reference_stem: str
    reference_n_leds: int


#: Reflector-free configurations compared in the manuscript and in Supporting Information S8.
#: Configurations with reflectors are excluded because the geometric model does not describe
#: reflective components.
CASES: dict[str, list[Case]] = {
    "LST1-01G01-UV01-00": [
        Case("best manual", "best_manual_no_reflector", None, 8, "best_manual_no_reflector", 8),
        Case("UV1", "UV1_no_reflector",
             "results/sampled_configs/near_optimal_configs/LST1-01G01-UV01-00/UV1.json", 16,
             "UV3_reflector", 16),
        Case("UV2", "UV2_no_reflector",
             "results/sampled_configs/near_optimal_configs/LST1-01G01-UV01-00/UV2.json", 16,
             "UV3_reflector", 16),
        Case("UV3", "UV3_no_reflector",
             "results/sampled_configs/near_optimal_configs/LST1-01G01-UV01-00/UV3.json", 16,
             "UV3_reflector", 16),
    ],
    "LST1-01F06-GRN1-00": [
        Case("best manual", "best_manual_no_reflector", None, 8, "best_manual_no_reflector", 8),
        Case("GRN1", "GRN1_no_reflector",
             "results/sampled_configs/near_optimal_configs/LST1-01F06-GRN1-00/GRN1.json", 16,
             "GRN3_reflector", 16),
        Case("GRN2", "GRN2_no_reflector",
             "results/sampled_configs/near_optimal_configs/LST1-01F06-GRN1-00/GRN2.json", 16,
             "GRN3_reflector", 16),
        Case("GRN3", "GRN3_no_reflector",
             "results/sampled_configs/near_optimal_configs/LST1-01F06-GRN1-00/GRN3.json", 16,
             "GRN3_reflector", 16),
    ],
}

#: The three pairwise comparisons, in the order used in the figures.
PAIRS = [
    ("geometric model", "ray tracing"),
    ("geometric model", "radiometry"),
    ("ray tracing", "radiometry"),
]

LEDS = {"LST1-01G01-UV01-00": LST1_01G01_UV01_00, "LST1-01F06-GRN1-00": LST1_01F06_GRN1_00}


def common_grid(resolution: float = RESOLUTION) -> tuple[np.ndarray, np.ndarray]:
    """
    Create the common evaluation grid on which all three methods are compared.

    :param resolution: Grid resolution in centimeters. Default is :data:`RESOLUTION`.
    :type resolution: float, optional

    :returns: 1D arrays of the x and y coordinates of the grid points in centimeters.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    x = np.arange(-PLANE_WIDTH / 2 + resolution / 2, PLANE_WIDTH / 2, resolution)
    y = np.arange(-PLANE_HEIGHT / 2 + resolution / 2, PLANE_HEIGHT / 2, resolution)
    return x, y


def read_configuration(path: str | None, results_path: str = "results/") -> tuple[list, float]:
    """
    Read LED positions and module height from a configuration JSON file.

    :param path: Path of the configuration file relative to the repository root. ``None`` selects
        the best manual configuration.
    :type path: str | None
    :param results_path: Path of the results directory, used when ``path`` is ``None``.
    :type results_path: str, optional

    :returns: List of (x, y) LED positions in centimeters and the module height in centimeters.
    :rtype: tuple[list, float]

    :raises FileNotFoundError: If the configuration file does not exist.
    """
    file = path if path is not None else os.path.join(results_path, "sampled_configs/best_manual.json")
    with open(file, encoding="utf-8") as handle:
        data = json.load(handle)
    return [tuple(position) for position in data["(x, y) Positions in [cm]"]], float(data["Height in [cm]"])


def regrid_interpolated(x_src: np.ndarray, y_src: np.ndarray, values: np.ndarray,
                        x_dst: np.ndarray, y_dst: np.ndarray) -> np.ndarray:
    """
    Interpolate a field bilinearly onto the destination grid.

    :param x_src: x coordinates of the source grid in centimeters (shape ``(nx,)``).
    :type x_src: numpy.ndarray
    :param y_src: y coordinates of the source grid in centimeters (shape ``(ny,)``).
    :type y_src: numpy.ndarray
    :param values: Field values on the source grid (shape ``(ny, nx)``), in W m\\ :sup:`-2`.
    :type values: numpy.ndarray
    :param x_dst: x coordinates of the destination grid in centimeters.
    :type x_dst: numpy.ndarray
    :param y_dst: y coordinates of the destination grid in centimeters.
    :type y_dst: numpy.ndarray

    :returns: Field on the destination grid; points outside the source grid are ``numpy.nan``.
    :rtype: numpy.ndarray
    """
    interpolator = RegularGridInterpolator((y_src, x_src), values, method="linear",
                                           bounds_error=False, fill_value=np.nan)
    X, Y = np.meshgrid(x_dst, y_dst)
    return interpolator(np.stack([Y.ravel(), X.ravel()], axis=-1)).reshape(Y.shape)


def regrid_averaged(x_src: np.ndarray, y_src: np.ndarray, values: np.ndarray,
                    x_dst: np.ndarray, y_dst: np.ndarray) -> np.ndarray:
    """
    Block-average a fine field onto a coarser destination grid.

    All source cells falling into one destination cell are averaged. This is used for the ray
    tracing fields, whose 1 mm mesh carries Monte-Carlo noise that would otherwise dominate the
    residual maps, while systematic deviations are unaffected by the averaging.

    :param x_src: x coordinates of the source grid in centimeters.
    :type x_src: numpy.ndarray
    :param y_src: y coordinates of the source grid in centimeters.
    :type y_src: numpy.ndarray
    :param values: Field values on the source grid (shape ``(ny, nx)``), in W m\\ :sup:`-2`.
    :type values: numpy.ndarray
    :param x_dst: x coordinates of the destination grid in centimeters.
    :type x_dst: numpy.ndarray
    :param y_dst: y coordinates of the destination grid in centimeters.
    :type y_dst: numpy.ndarray

    :returns: Averaged field on the destination grid.
    :rtype: numpy.ndarray
    """
    X, Y = np.meshgrid(x_src, y_src)
    dx = float(np.mean(np.diff(x_dst)))
    dy = float(np.mean(np.diff(y_dst)))
    x_edges = np.append(x_dst - dx / 2, x_dst[-1] + dx / 2)
    y_edges = np.append(y_dst - dy / 2, y_dst[-1] + dy / 2)
    return binned_statistic_2d(Y.ravel(), X.ravel(), values.ravel(), statistic="mean",
                               bins=[y_edges, x_edges]).statistic


def field_geometric(led, positions: list, height: float, x_dst: np.ndarray,
                    y_dst: np.ndarray) -> np.ndarray:
    """
    Evaluate the geometric radiation field model on the common grid.

    :param led: LED object used for the simulation.
    :type led: :class:`utils.led.LED`
    :param positions: LED positions in centimeters.
    :type positions: list
    :param height: Distance between the LEDs and the detection plane in centimeters.
    :type height: float
    :param x_dst: x coordinates of the common grid in centimeters.
    :type x_dst: numpy.ndarray
    :param y_dst: y coordinates of the common grid in centimeters.
    :type y_dst: numpy.ndarray

    :returns: Irradiance field in W m\\ :sup:`-2` on the common grid, mirrored on the y axis so
        that it is directly comparable with the ray tracing and radiometry fields.
    :rtype: numpy.ndarray

    .. note::
        The geometric model describes the detection plane as seen from the opposite side compared
        with the ray tracing export and the radiometric scan, so the simulated field is mirrored on
        the y axis before it is returned. The same convention is applied in
        :func:`visualization.irradiance_plots_SI.create_irradiance_plots_SI`. Mean irradiance and
        homogeneity are invariant under this operation; only the point-by-point residuals depend
        on it.
    """
    grid = Grid(width=PLANE_WIDTH, height=PLANE_HEIGHT, step=2.5, side_space=1.5, top_bottom_space=2.0)
    model = GeometricModel(grid, led, resolution_xy=RESOLUTION)
    irradiance = model.simulate(list(positions), height)
    # mirror the model field on the y axis for comparability with ray tracing and radiometry
    irradiance = np.flip(irradiance, axis=0)
    return regrid_interpolated(model.X[0, :], model.Y[:, 0], irradiance, x_dst, y_dst)


def field_raytracing(led_name: str, stem: str, x_dst: np.ndarray, y_dst: np.ndarray,
                     results_path: str = "results/") -> np.ndarray:
    """
    Load a ray tracing result and block-average it onto the common grid.

    :param led_name: Name of the LED, used to build the data path.
    :type led_name: str
    :param stem: File stem of the ray tracing export without extension.
    :type stem: str
    :param x_dst: x coordinates of the common grid in centimeters.
    :type x_dst: numpy.ndarray
    :param y_dst: y coordinates of the common grid in centimeters.
    :type y_dst: numpy.ndarray
    :param results_path: Path of the results directory.
    :type results_path: str, optional

    :returns: Irradiance field in W m\\ :sup:`-2` on the common grid.
    :rtype: numpy.ndarray

    :raises FileNotFoundError: If the ray tracing file does not exist.
    """
    data = process_raytracing_file(path=os.path.join(results_path, "raytracing", led_name),
                                   filename=f"{stem}.txt")
    pivot = data.pivot_table(index="y", columns="x", values="Irradiance")
    return regrid_averaged(pivot.columns.to_numpy(float), pivot.index.to_numpy(float),
                           pivot.to_numpy(float), x_dst, y_dst)


def field_radiometry(led, case: Case, x_dst: np.ndarray, y_dst: np.ndarray,
                     results_path: str = "results/") -> np.ndarray:
    """
    Load a radiometric scan, apply the photon-flux correction and interpolate onto the common grid.

    :param led: LED object of the measurement.
    :type led: :class:`utils.led.LED`
    :param case: Configuration to be processed.
    :type case: :class:`Case`
    :param x_dst: x coordinates of the common grid in centimeters.
    :type x_dst: numpy.ndarray
    :param y_dst: y coordinates of the common grid in centimeters.
    :type y_dst: numpy.ndarray
    :param results_path: Path of the results directory.
    :type results_path: str, optional

    :returns: Irradiance field in W m\\ :sup:`-2` on the common grid.
    :rtype: numpy.ndarray

    :raises FileNotFoundError: If the radiometry file does not exist.
    """
    folder = os.path.join(results_path, "radiometry", led.name)
    correction_factor = calculate_correction_factor(
        os.path.join(folder, f"{case.reference_stem}.csv"), case.reference_n_leds * led.total_power)
    frame, _, _ = process_radiometry_file(os.path.join(folder, f"{case.stem}.csv"),
                                          correction_factor, height=PLANE_HEIGHT, width=PLANE_WIDTH)
    return regrid_interpolated(frame.columns.to_numpy(float), frame.index.to_numpy(float),
                               frame.to_numpy(float), x_dst, y_dst)


def residual_statistics(delta: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    """
    Summarize a residual field by bias, root mean square error and maximum deviation.

    :param delta: Residual field in W m\\ :sup:`-2`.
    :type delta: numpy.ndarray
    :param reference: Field of the second method of the pair, used for the relative values.
    :type reference: numpy.ndarray

    :returns: Dictionary with the absolute values in W m\\ :sup:`-2` and the relative values in
        percent of the mean irradiance of the reference.
    :rtype: dict[str, float]
    """
    mask = np.isfinite(delta) & np.isfinite(reference)
    delta, reference = delta[mask], reference[mask]
    reference_mean = float(np.mean(reference))
    bias = float(np.mean(delta))
    rmse = float(np.sqrt(np.mean(delta ** 2)))
    max_abs = float(np.max(np.abs(delta)))
    return {
        "bias / W m-2": bias,
        "RMSE / W m-2": rmse,
        "max |dE| / W m-2": max_abs,
        "bias / %": 100.0 * bias / reference_mean,
        "RMSE / %": 100.0 * rmse / reference_mean,
        "max |dE| / %": 100.0 * max_abs / reference_mean,
    }


def compare_case(led, case: Case, x: np.ndarray, y: np.ndarray,
                 results_path: str = "results/") -> tuple[dict[str, np.ndarray], list[dict]]:
    """
    Evaluate all three methods for one configuration and compute the pairwise residuals.

    :param led: LED object of the configuration.
    :type led: :class:`utils.led.LED`
    :param case: Configuration to be processed.
    :type case: :class:`Case`
    :param x: x coordinates of the common grid in centimeters.
    :type x: numpy.ndarray
    :param y: y coordinates of the common grid in centimeters.
    :type y: numpy.ndarray
    :param results_path: Path of the results directory.
    :type results_path: str, optional

    :returns: Dictionary mapping the name of each comparison to its residual field, and a list of
        one summary dictionary per comparison.
    :rtype: tuple[dict[str, numpy.ndarray], list[dict]]
    """
    positions, height = read_configuration(case.config_json, results_path)
    fields = {
        "geometric model": field_geometric(led, positions, height, x, y),
        "ray tracing": field_raytracing(led.name, case.stem, x, y, results_path),
        "radiometry": field_radiometry(led, case, x, y, results_path),
    }

    residuals, rows = {}, []
    for first, second in PAIRS:
        delta = fields[first] - fields[second]
        residuals[f"{first} - {second}"] = delta
        rows.append({
            "LED": led.name,
            "configuration": case.label,
            "height / cm": height,
            "comparison": f"{first} - {second}",
            "E_mean(A) / W m-2": float(np.nanmean(fields[first])),
            "E_mean(B) / W m-2": float(np.nanmean(fields[second])),
            "H(A) / 1": 1.0 - float(np.nanstd(fields[first])) / float(np.nanmean(fields[first])),
            "H(B) / 1": 1.0 - float(np.nanstd(fields[second])) / float(np.nanmean(fields[second])),
            **residual_statistics(delta, fields[second]),
        })
    return residuals, rows


def compare_all(results_path: str = "results/") -> tuple[dict, pd.DataFrame]:
    """
    Compare all configurations of both LED types and collect residual fields and statistics.

    :param results_path: Path of the results directory.
    :type results_path: str, optional

    :returns: Dictionary ``{led_name: {configuration_label: residual fields}}`` and the table of
        residual statistics.
    :rtype: tuple[dict, pandas.DataFrame]
    """
    x, y = common_grid()
    all_residuals: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    statistics: list[dict] = []
    for led_name, cases in CASES.items():
        led = LEDS[led_name]
        all_residuals[led_name] = {}
        for case in cases:
            residuals, rows = compare_case(led, case, x, y, results_path)
            all_residuals[led_name][case.label] = residuals
            statistics.extend(rows)
    return all_residuals, pd.DataFrame(statistics)
