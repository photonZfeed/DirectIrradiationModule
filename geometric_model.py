import json
import numpy as np

from time import time
from concurrent.futures import ProcessPoolExecutor, as_completed 
from numba import njit, prange
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.colors as mcolors
from utils.grid import Grid
from utils.led import LED, LST1_01F06_GRN1_00, LST1_01G01_UV01_00
from tqdm import tqdm
from pathlib import Path


class GeometricModel:
    """
    Geometric model for simulating direct irradiation from multiple LED sources over a 2D spatial grid.

    This class is used to simulate the irradiance distribution from configurable LED arrangements on a defined grid. It supports single configuration simulations, 
    batch processing of multiple configurations (with parallelization), and visualization of irradiance distribution on the detection plane.

    :param G: Grid object defining the area of interest, including width, height, and spatial discretization (units: cm).
    :type G: Grid
    :param led: LED source to be used in the simulation. Must provide angle and intensity data and implement ``calc_radiant_intensity()``.
    :type led: LED
    :param resolution_xy: Spatial resolution of the grid in centimeters (cm). Determines spacing between grid points. Default is 0.5 cm.
    :type resolution_xy: float, optional

    :ivar h: Planck constant (Joule seconds, J·s).
    :vartype h: float
    :ivar c: Speed of light (meters per second, m/s).
    :vartype c: float
    :ivar Na: Avogadro's number (1/mol).
    :vartype Na: float
    :ivar led: The LED object used for simulation.
    :vartype led: LED
    :ivar X: 2D array of x-coordinates for the simulation grid (units: cm).
    :vartype X: numpy.ndarray
    :ivar Y: 2D array of y-coordinates for the simulation grid (units: cm).
    :vartype Y: numpy.ndarray

    :raises ValueError: If the provided grid or LED object is invalid or missing required attributes.

    :side effects:
        Calls ``calc_radiant_intensity()`` on the provided LED object, which may modify its internal state.

    .. note::
        - The grid is centered at (0, 0) in the x-y plane.
        - All spatial units are in centimeters unless otherwise specified.
        - The LED object must provide angle and intensity data arrays.
        - Designed for use with batch simulation and parallel processing for a large number of configurations.
        - For visualization, use the :meth:`plot_irradiance` method with a Matplotlib Axes object as an argument.
    """

    def __init__(self, G: Grid, led: LED, resolution_xy: float = 0.5):
        """
        Initialize the geometric model with a grid, LED, and spatial resolution.

        :param G: Grid object defining the area of interest (width, height, step size, etc.; units: cm).
        :type G: Grid
        :param led: LED source to be used in the simulation. Must provide angle and intensity data and implement ``calc_radiant_intensity()``.
        :type led: LED
        :param resolution_xy: Spatial resolution of the grid in centimeters (cm). Default is 0.5 cm.
        :type resolution_xy: float, optional

        :raises ValueError: If the provided grid or LED object is invalid or missing required attributes.

        :side effects:
            Calls ``calc_radiant_intensity()`` on the provided LED object, which may modify its internal state.

        .. note::
            The grid is centered at (0, 0) in the x-y plane. All spatial units are in centimeters unless otherwise specified.
        """


        self.h = 6.62607015e-34 # Planck constant (J s)
        self.c = 299792458 # Speed of light (m/s)
        self.Na = 6.02214076e23 # Avogadro's number (1/mol)

        self.led = led
        self.led.calc_radiant_intensity()

        self.X, self.Y = self._make_grid(G.width, G.height, resolution_xy)

    def _make_grid(self, width: float, height: float, resolution_xy: float):
        """
        Construct a 2D meshgrid for the area of interest, centered at (0, 0).

        :param width: Width of the area (cm).
        :type width: float
        :param height: Height of the area (cm).
        :type height: float
        :param resolution_xy: Spatial resolution in centimeters (cm).
        :type resolution_xy: float

        :return: Tuple of meshgrid arrays for X and Y coordinates (both numpy.ndarray, units: cm).
        :rtype: tuple[numpy.ndarray, numpy.ndarray]

        .. note::
            The grid is constructed such that its center is at (0, 0).
        """

        x = np.arange(-width/2 + resolution_xy/2, width/2 + resolution_xy/2, resolution_xy)
        y = np.arange(-height/2 + resolution_xy/2, height/2 + resolution_xy/2, resolution_xy)
        return np.meshgrid(x, y, indexing="xy")

    def read_configurations(self, filename: str):
        """
        Load LED array configurations and source heights from a compressed NumPy (.npz) file.

        :param filename: Path to the .npz file containing configurations and heights.
        :type filename: str

        :return: Tuple containing the list/array of configurations and the corresponding heights (z values).
        :rtype: tuple[list, float or numpy.ndarray]

        :raises FileNotFoundError: If the file does not exist.
        :raises KeyError: If required keys are missing in the file.

        .. note::
            The .npz file must contain arrays named 'configurations' and 'height'.
        """
        data = np.load(filename, allow_pickle=True)
        configurations = data['configurations']
        z = data['height']
        return configurations, z

    def simulate(self, config: list, z: float):
        """
        Compute the total irradiance over the 2D grid from multiple LED sources at specified positions and height.

        :param config: List of tuples representing (x, y) positions (in cm) of the LED sources.
        :type config: list[tuple[float, float]]
        :param z: Height of the sources above the grid (cm).
        :type z: float

        :return: 2D array of total irradiance values over the grid (units: W/m²).
        :rtype: numpy.ndarray

        :side effects:
            Calls the ``calculate_irradiance`` function for each LED source.

        .. note::
            The result is converted from W/cm² to W/m² by multiplying by 1e4.
        """
        irr_total = np.zeros_like(self.X, dtype=np.float64)
        for pos in config:
            X_shift = self.X - pos[0]
            Y_shift = self.Y - pos[1]
            irr_led = calculate_irradiance(X_shift, Y_shift, z, self.led.angle_data, self.led.intensity_data)
            irr_total += irr_led
        return irr_total * 1e4  # convert from W/cm² to W/m²

    def _batch_worker(self, config_batch, z_batch):
        """
        Simulate a batch of configurations, returning summary irradiance metrics for each.

        :param config_batch: List of configurations (each a list of (x, y) tuples) for the batch.
        :type config_batch: list
        :param z_batch: List or array of heights (cm) corresponding to each configuration.
        :type z_batch: list or numpy.ndarray

        :return: Dictionary containing configurations, heights, and irradiance metrics (mean, max, std for each configuration).
        :rtype: dict

        .. note::
            Used internally for parallel batch processing.
        """
        batch_results = {
            "configs": config_batch,
            "z": z_batch,
            "means": [],
            "max": [],
            "stds": [],
        }
        for config, z in zip(config_batch, z_batch):
            irr = self.simulate(config, z)
            batch_results["means"].append(np.mean(irr))
            batch_results["max"].append(np.max(irr))
            batch_results["stds"].append(np.std(irr))
        return batch_results

    def _worker_wrapper(self, args):
        """
        Unpack arguments and call the batch worker for parallel processing.

        :param args: Tuple of (config_batch, z_batch) for the worker.
        :type args: tuple

        :return: Result of :meth:`_batch_worker`.
        :rtype: dict

        .. note::
            Used internally for parallel batch processing.
        """
        return self._batch_worker(*args)

    def simulate_batch(self, configurations, z, batch_size=500, max_workers=None):
        """
        Simulate a batch of configurations in parallel, returning summary irradiance metrics for each configuration.

        :param configurations: List of configurations to simulate (each a list of (x, y) tuples).
        :type configurations: list
        :param z: Height(s) of the sources above the grid (cm). A scalar float is broadcast to all
            configurations; a 1-D array must have one entry per configuration.
        :type z: float or numpy.ndarray
        :param batch_size: Number of configurations per batch. Default is 500.
        :type batch_size: int, optional
        :param max_workers: Maximum number of parallel workers. Default is None (uses os.cpu_count()).
        :type max_workers: int, optional

        :return: Dictionary containing configurations, heights, and irradiance metrics (mean, max, std for each configuration).
        :rtype: dict

        :side effects:
            Uses multiprocessing for parallel execution; may increase memory usage and CPU load.

        :raises Exception: Propagates exceptions from worker processes.

        .. note::
            Designed for efficient batch simulation of large configuration sets. Results are aggregated from all batches.
            A plain Python ``float`` or a zero-dimensional array is accepted for ``z`` and is automatically
            broadcast to the length of ``configurations``.

        .. rubric:: Example

        .. code-block:: python

            from utils.grid import Grid
            from utils.led import LST1_01G01_UV01_00
            from geometric_model import GeometricModel, save_results

            G = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
            model = GeometricModel(G, LST1_01G01_UV01_00, resolution_xy=0.5)

            configs, z = model.read_configurations(
                "results/examples/example_MonteCarlo_1-16_leds_5-15_cm_10000_samples.npz"
            )
            results = model.simulate_batch(configs, z, batch_size=500)
            save_results(results, model, "results/examples/batch_results.npz")
        """
        z = np.atleast_1d(np.asarray(z, dtype=np.float64))
        if z.size == 1:
            z = np.full(len(configurations), z)
        tasks = [
            (configurations[i:i+batch_size], z[i:i+batch_size])
            for i in range(0, len(configurations), batch_size)
        ]
        results = {"configs": [], "z": [], "means": [], "max": [], "stds": []}
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(self._worker_wrapper, task) for task in tasks]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Simulating Batches"):
                batch_result = future.result()
                results["configs"].extend(batch_result["configs"])
                results["z"].extend(batch_result["z"])
                results["means"].extend(batch_result["means"])
                results["max"].extend(batch_result["max"])
                results["stds"].extend(batch_result["stds"])
        return results
    
    def plot_irradiance(self, ax: Axes, irr : np.ndarray, title: str = "", y_title: float = -0.4, v_max: float = -1., extent: tuple = (-2, -1, -2, -1), shrink_cbar: float = 1.0):
        """
        Visualize the irradiance distribution as a contour plot on a given Matplotlib Axes.

        :param ax: Matplotlib Axes object to plot on.
        :type ax: matplotlib.axes.Axes
        :param irr: 2D array containing the irradiance data to plot (units: W/m²).
        :type irr: numpy.ndarray
        :param title: Title of the plot. Default is empty string.
        :type title: str, optional
        :param y_title: Y position of the title. Default is -0.4.
        :type y_title: float, optional
        :param v_max: Maximum value for color scale. Default is -1 (auto scale).
        :type v_max: float, optional
        :param extent: Extent of the plot in the form (xmin, xmax, ymin, ymax). Default is (-2, -1, -2, -1) (auto extent).
        :type extent: tuple, optional
        :param shrink_cbar: Factor to shrink the color bar. Default is 1.0.
        :type shrink_cbar: float, optional

        :side effects:
            Modifies the provided Axes object and creates a colorbar. Does not return a value.

        .. note::
            The color scale and contour levels are automatically determined based on the data unless specified.
        """

        # Set color bar limits
        if v_max != -1.:
            vmax = v_max
        else:
            vmax = (int(irr.max() // 5) + 1) * 5 # round up to next multiple of 5
        norm = mcolors.Normalize(vmin=0, vmax=vmax, clip=False)

        # define contour levels based on vmax
        if vmax <= 30:
            step = 5
        elif vmax <= 50:
            step = 10
        elif vmax <= 100:
            step = 20
        else:
            step = 50
        levels = np.arange(0, vmax + step, step)

        if extent == (-2, -1, -2, -1):
            extent = (
                self.X.min(),
                self.X.max(),
                self.Y.min(),
                self.Y.max()
            )

        # Create colormesh plot
        tc = ax.imshow(
            irr,
            extent=extent,
            origin='lower',
            cmap='viridis',
            norm=norm,
            interpolation='nearest',
            aspect='equal'
        )

        # Create colorbar
        cbar = plt.colorbar(tc, ax=ax, ticks=levels, shrink=shrink_cbar)
        cbar.set_label('$E$ / W m⁻²')

        ax.set_xlabel('$X$ / cm')
        ax.set_ylabel('$Y$ / cm')
        ax.set_aspect('equal', 'box')

        # Set x and y limits
        ax.set_xlim(self.X.min(), self.X.max())
        ax.set_ylim(self.Y.min(), self.Y.max())

        # set subtitle
        if title:
            ax.set_title(title, y=y_title, fontsize=10)

   
@njit(fastmath=True)
def calculate_angle(x, y, z):
    """
    Compute the angle (theta) in radians between the z-axis and the point (x, y, z).

    :param x: X-coordinate (float)
    :type x: float
    :param y: Y-coordinate (float)
    :type y: float
    :param z: Z-coordinate (float)
    :type z: float

    :return: Angle theta in radians between the z-axis and the point (x, y, z).
    :rtype: float

    .. note::
        Used for determining the emission angle from the LED to a grid point.
    """
    return np.arctan(np.sqrt(x**2 + y**2) / z)

@njit(parallel=True, fastmath=True)
def calculate_irradiance(X, Y, z, angle_data, intensity_data):
    """
    Calculate the irradiance at each point in the grid from a single LED source.

    :param X: 2D array of x-coordinates (cm).
    :type X: numpy.ndarray
    :param Y: 2D array of y-coordinates (cm).
    :type Y: numpy.ndarray
    :param z: Height of the LED source above the grid (cm).
    :type z: float
    :param angle_data: 1D array of emission angles (degrees) for the LED.
    :type angle_data: numpy.ndarray
    :param intensity_data: 1D array of intensities corresponding to angle_data.
    :type intensity_data: numpy.ndarray

    :return: 2D array of irradiance values at each grid point (units: W cm⁻²).
    :rtype: numpy.ndarray

    .. note::
        The result is in W cm⁻². For W m⁻², multiply by 1e4.
        Uses linear interpolation for intensity as a function of angle.
    """
    rows, cols = X.shape
    result = np.empty((rows, cols), dtype=np.float64)
    for i in prange(rows):
        for j in range(cols):
            theta = np.degrees(calculate_angle(X[i, j], Y[i, j], z))
            intensity = np.interp(theta, angle_data, intensity_data)
            result[i, j] = z / ((X[i, j]**2 + Y[i, j]**2 + z**2)**1.5) * intensity
    return result

def save_irr(irr, model, filename):
    """
    Save the irradiance data of a single configuration to a compressed NumPy (.npz) file.

    :param irr: 2D array of irradiance values (W m⁻²).
    :type irr: numpy.ndarray
    :param model: The geometric model used for the simulation (provides grid and LED info).
    :type model: GeometricModel
    :param filename: Path to the .npz file to save the data.
    :type filename: str

    :side effects:
        Writes a file to disk at the specified location.

    .. note::
        The saved file contains arrays for X, Y, irr, led_name, and optical_power_per_led.
    """
    np.savez_compressed(
        filename,
        X=model.X,
        Y=model.Y,
        irr=irr,
        led_name=model.led.name,
        optical_power_per_led=model.led.total_power,
    )

def read_irr(filename):
    """
    Load irradiance data and associated metadata from a compressed NumPy (.npz) file.

    :param filename: Path to the .npz file containing irradiance data.
    :type filename: str

    :return: Tuple containing X (grid x-coordinates), Y (grid y-coordinates), irr (irradiance array), led_name (str), and optical_power_per_led (float).
    :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, str, float]

    :raises FileNotFoundError: If the file does not exist.
    :raises KeyError: If required keys are missing in the file.

    .. note::
        The file must contain arrays for X, Y, irr, led_name, and optical_power_per_led.
    """
    data = np.load(filename)
    X = data['X']
    Y = data['Y']
    irr = data['irr']
    led_name = data['led_name'].item()
    optical_power_per_led = data['optical_power_per_led'].item()
    return X, Y, irr, led_name, optical_power_per_led

def save_results(results, model, filename):
    """
    Save simulation results for all configurations to a compressed NumPy (.npz) file.

    :param results: Dictionary containing simulation results (configs, z, means, max, stds).
    :type results: dict
    :param model: The geometric model used for the simulation (provides LED info).
    :type model: GeometricModel
    :param filename: Path to the .npz file to save the data.
    :type filename: str

    .. note::
        The saved file contains arrays for configs, z, means, max, stds, led_name, and optical_power_per_led.
    """
    np.savez_compressed(
        filename,
        configs=np.array(results["configs"], dtype=object),
        z=np.array(results["z"]),
        means=np.array(results["means"]),
        max=np.array(results["max"]),
        stds=np.array(results["stds"]),
        led_name=model.led.name,
        optical_power_per_led=model.led.total_power,
    )

if __name__ == "__main__":

    # create grid, led, and model
    G = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
    led = LST1_01G01_UV01_00
    # led = LST1_01F06_GRN1_00
    model = GeometricModel(G, led, resolution_xy=0.5)

    # simulate batch of configurations and save results
    config_path = "results\\sampled_configs\\systematic_8_leds_13_cm_482931_samples.npz" # 8 leds, 13 cm height
    # config_path = "results\\sampled_configs\\MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz" # 1-16 leds, 5-15 cm height
    
    configs, z = model.read_configurations(config_path)
    configs = configs[:10]
    results = model.simulate_batch(configs, z, batch_size=500, max_workers=None)
    save_results(results, model, f"results/geometric/{led.name}/simulation_results_{config_path.split('\\')[-1]}")

    # simulate single configuration
    # config = [(-7.5,0), (7.5,0), (-10, -10), (10, 10), (-10, 10), (10, -10), (0,12.5), (0,-12.5)]
    # z = 13

    # irr = model.simulate(config, z)
    # print(f"Mean Irradiance: {np.mean(irr):.4f} W/m²")
    # print(f"Homogeneity: {1 - np.std(irr) / np.mean(irr):.4f}")
    # fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))
    # model.plot_irradiance(ax, irr, shrink_cbar=0.8)
    # plt.show()

    # print(f"Mean Irradiance: {np.mean(irr):.2f} W/m²")
    # print(f"Std Dev of Irradiance: {np.std(irr):.2f} W/m²")

    # track runtime of single simulation
    # simulate five times and measure time and print average time and standard deviation
    # sim_times = []
    # for i in range(5):
    #     start_time = time()
    #     model.simulate(config, z)
    #     end_time = time()
    #     sim_times.append(end_time - start_time)

    # # convert to array of milliseconds
    # sim_times = np.array(sim_times) * 1000
    # print(f"Simulation Time: {np.mean(sim_times):.2f} ± {np.std(sim_times):.2f} milliseconds")