import numpy as np


class LED:
    """
    Represents a light-emitting diode (LED) with methods to calculate its absolute radiant intensity distribution
    from relative angular intensity data.

    This class loads relative intensity data from a file, assumes rotational symmetry around the optical axis,
    and computes the absolute radiant intensity distribution based on the total optical power.

    Attributes
    ----------
    name : str
        Name or identifier of the LED.
    total_power : float
        Total optical power of the LED in watts (W).
    intensity_data_file : str
        Path to the CSV file containing angular intensity data.
    angle_data : np.ndarray
        Array of angles (in degrees) loaded from the intensity data file.
    intensity_data : np.ndarray
        Array of absolute radiant intensity values (in W/sr) after calculation.

    Examples
    --------
    >>> led = LED("Example-LED", 0.5, "path/to/intensity.csv")
    >>> abs_intensity = led.calc_radiant_intensity()
    >>> print(abs_intensity)
    """

    def __init__(self, name, optical_power, intensity_data_file):
        """
        Initialize an LED instance with its name, total optical power, and path to intensity data file.

        Parameters
        ----------
        name : str
            Name or identifier for the LED.
        optical_power : float
            Total optical power of the LED in watts (W).
        intensity_data_file : str
            Path to the CSV file containing angular intensity data. The file must have at least two columns:
            angle (degrees) and relative intensity (arbitrary units), with a header row.

        Notes
        -----
        The angle and intensity data are initialized as empty arrays and populated when `calc_radiant_intensity` is called.
        """
        self.name = name
        self.total_power = optical_power  # in watts
        self.intensity_data_file = intensity_data_file  # path to the intensity data file
        self.angle_data = np.array([])  # placeholder for angle data
        self.intensity_data = np.array([])  # placeholder for intensity data

    def calc_radiant_intensity(self):
        """
        Calculate the absolute radiant intensity distribution from relative angular intensity data.

        Loads the relative intensity data from the specified CSV file, assumes rotational symmetry around the optical axis,
        and computes the absolute radiant intensity (in W/sr) for each angle using the total optical power.

        The CSV file must have at least two columns: angle (degrees) and relative intensity (arbitrary units),
        with a header row that is skipped during loading.

        Returns
        -------
        np.ndarray
            Array of absolute radiant intensity values (in W/sr) corresponding to each angle in `self.angle_data`.

        Raises
        ------
        IOError
            If the intensity data file cannot be read.
        ValueError
            If the intensity data file does not have the expected format or contains invalid data.

        Side Effects
        ------------
        Updates `self.angle_data` and `self.intensity_data` with the loaded angles and calculated absolute intensities.

        Notes
        -----
        - Assumes rotational symmetry around the optical axis.
        - Uses the trapezoidal rule for numerical integration.
        - The calculation normalizes the relative intensity data such that the total emitted power matches `self.total_power`.

        Examples
        --------
        >>> led = LED("Example-LED", 0.5, "path/to/intensity.csv")
        >>> abs_intensity = led.calc_radiant_intensity()
        >>> print(abs_intensity)
        """
        try:
            intensity_data = np.loadtxt(self.intensity_data_file, dtype=np.float64, delimiter=",", skiprows=1)
        except OSError as e:
            raise IOError(f"Could not read intensity data file: {self.intensity_data_file}") from e
        except ValueError as e:
            raise ValueError(f"Invalid format in intensity data file: {self.intensity_data_file}") from e

        self.angle_data = intensity_data[:, 0]  # angles in degrees
        theta_rad = np.deg2rad(self.angle_data)

        # Assume rotational symmetry around optical axis
        integrand = intensity_data[:, 1] * np.sin(theta_rad)
        integral = np.trapezoid(integrand, theta_rad)  # numerical integration
        I_max = self.total_power / (2 * np.pi * integral)

        # Calculate absolute intensity
        self.intensity_data = I_max * intensity_data[:, 1]  # update intensity data with absolute values
        return self.intensity_data

# Define specific LED instances
LST1_01F06_GRN1_00 = LED(
    name="LST1-01F06-GRN1-00",
    optical_power=0.281,  # @ 350 mA in W
    intensity_data_file="utils/angle_intensity_data/LST1-01F06-GRN1-00.csv"
)

LST1_01G01_UV01_00 = LED(
    name="LST1-01G01-UV01-00",
    optical_power=0.875,  # @ 500 mA in W
    intensity_data_file="utils/angle_intensity_data/LST1-01G01-UV01-00.csv"
)