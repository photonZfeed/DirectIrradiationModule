from systematic_sampler import Sampler1
import matplotlib.pyplot as plt

def plot_cand(conf_name):
    """
    Generates and saves a plot visualizing a candidate configuration for a radiation field experiment.

    This function loads a configuration from a JSON file, visualizes it using matplotlib, and saves the resulting figure as an SVG file in the 'figures' directory. The configuration file path is determined by the content of ``conf_name``. The function is intended for visual inspection and documentation of candidate configurations.

    :param conf_name: Name of the candidate configuration (str). Used to determine the configuration file path and output filename. Should contain either 'UV', 'GRN', or neither, which selects the corresponding configuration directory or a default manual configuration.
    :type conf_name: str

    :returns: None

    :raises FileNotFoundError: If the configuration JSON file does not exist at the resolved path.
    :raises ValueError: If the configuration file cannot be parsed or is invalid.

    :side effects:
        - Saves a plot as an SVG file in the 'figures' directory with the name ``<conf_name>.svg``.

    .. note::
        - The function assumes that the directories 'results/sampled_configs/' and 'figures/' exist and are writable.
        - The configuration JSON files must be compatible with ``Sampler1.read_configuration_from_json``.
        - The function does not display the plot; it only saves it to disk.

    **Example usage:**

    .. code-block:: python

        plot_cand('LST1-01G01-UV01-00-001')
        plot_cand('LST1-01F06-GRN1-00-002')
        plot_cand('best_manual')
    """
    if "UV" in conf_name:
        conf_path = f"results/sampled_configs/near_optimal_configs/LST1-01G01-UV01-00/{conf_name}.json"
    elif "GRN" in conf_name:
        conf_path = f"results/sampled_configs/near_optimal_configs/LST1-01F06-GRN1-00/{conf_name}.json"
    else:
        conf_path = f"results/sampled_configs/best_manual.json"
    
        
    filename = f"figures/{conf_name}.svg"
    fig, ax = plt.subplots(figsize=(2,2))
    conf, height = Sampler1.read_configuration_from_json(conf_path)
    ax = Sampler1.plot_configuration(ax, conf, z=height, point_size=20, set_title=False)
    fig.savefig(filename, bbox_inches='tight')