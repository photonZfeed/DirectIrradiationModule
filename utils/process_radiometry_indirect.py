import pandas as pd

def process_radiometry_data(filepath: str, height: float = 0.33, width: float = 0.34) -> tuple[pd.DataFrame, float, float]:
    """
    Processes radiometry data for indirect irradiation experiments, following the methodology described in Kowalczyk et al. (2023), DOI: 10.1039/D3RE00398A. This function reads radiometry data from a CSV file, determines the mean irradiance based on published received optical power for specific LED types, and calculates the homogeneity of the irradiance distribution. The function is tailored for datasets with a specific format and LED types as described in the referenced publication.

    The CSV file is expected to contain columns for x, y, and z coordinates (with z representing measured irradiance values), separated by semicolons and without a header row. The function automatically crops the data to a region of interest based on published values and applies a correction factor for certain LED types.

    :param filepath: Path to the CSV file containing radiometry data. The file must use semicolon (';') as a delimiter and contain columns for x, y, and z (irradiance) values. The LED type must be identifiable from the filepath string (see Notes).
    :type filepath: str
    :param height: Height of the detection plane in meters. Default is 0.33 m.
    :type height: float, optional
    :param width: Width of the detection plane in meters. Default is 0.34 m.
    :type width: float, optional

    :returns: A tuple containing:
        - irr (pandas.DataFrame): Cropped irradiance data with columns ['x', 'y', 'z'].
        - mean_irradiance (float): Mean irradiance in W/m², calculated from published optical power and detection plane area.
        - homogeneity (float): Homogeneity of the irradiance distribution, defined as 1 - (standard deviation / mean) of the irradiance values.
    :rtype: tuple[pandas.DataFrame, float, float]

    :raises ValueError: If the LED type cannot be determined from the filepath.

    .. note::
        - Only supports LED types "LST1-01G01-UV01-00" and "LST1-01F06-GRN1-00" as identified in the filepath. For the green LED, a correction factor is applied to estimate optical power at 350 mA from 1 A data.
        - The region of interest is cropped to x in [11, 45] and y in [16, 49] as per the referenced publication.
        - The CSV file must not contain a header row; the first row is skipped.
        - The returned ``irr`` DataFrame has columns ``['x', 'y', 'z']`` (where ``z`` is the irradiance). To visualize it, use :func:`utils.process_raytracing_data.plot_irradiance_raytracing` after renaming the ``z`` column to ``Irradiance``.

    :side effects:
        None. The function only reads from the specified file and returns processed data.

    :references:
        Kowalczyk et al., Reaction Chemistry & Engineering, 2023, DOI: 10.1039/D3RE00398A
    """


    # Read the CSV file into a DataFrame
    irr = pd.read_csv(filepath, delimiter=';', header=None, names=['x', 'y', 'z'], skiprows=1)

    # get led type from filepath
    if "LST1-01G01-UV01-00" in filepath:
        P_opt = 1.88 # W # optical power measured on the detection plane according to doi.org/10.1039/D3RE00398A
    elif "LST1-01F06-GRN1-00" in filepath:
        P_opt = 1.032 # W # optical power measured on the detection plane according to doi.org/10.1039/D3RE00398A
        corr_factor = 0.5 # correction factor for green LED according to estimate optical power at 350 mA from the 1 A data in doi.org/10.1039/D3RE00398A. Taken from led datasheet 
        P_opt *= corr_factor
    else:
        P_opt = 0
        raise ValueError("Unknown LED type in filepath")


    # crop values according to values from doi.org/10.1039/D3RE00398A
    x_min = 11
    x_max = 45
    y_min = 16
    y_max = 49

    irr = irr[(irr['x'] >= x_min) & (irr['x'] <= x_max) & (irr['y'] >= y_min) & (irr['y'] <= y_max)]

    # Calculate mean irradiance I_mean
    mean_irradiance = P_opt / (height * width)  # W/m^2
    homogeneity = 1 - irr['z'].std() / irr['z'].mean()

    return irr, mean_irradiance, homogeneity
