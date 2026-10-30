"""
Plot style shared by all figures of the manuscript.

The figures use the ``ICIWstyle`` of the group's ICIW-Plots package and the repository's
``visualization/publication_style.mplstyle``. Both are applied by file path, so the style works on
every operating system: ICIW-Plots 1.0.42 registers ``ICIWstyle`` by name only on Windows.

**Usage Example**::

    from visualization.style import use_publication_style, finish
    use_publication_style()
    fig, ax = plt.subplots()
    ...
    finish(fig, show=False)
"""
from pathlib import Path

import ICIW_Plots
import matplotlib.pyplot as plt

#: The repository's publication style, next to this module.
PUBLICATION_STYLE = Path(__file__).with_name("publication_style.mplstyle")


def use_publication_style() -> None:
    """Apply ``ICIWstyle`` and the publication style of the repository."""
    plt.style.use([str(Path(ICIW_Plots.__file__).parent / "ICIWstyle.mplstyle"), str(PUBLICATION_STYLE)])


def finish(fig=None, show: bool = True) -> None:
    """
    Display a figure, or close it when it is only written to disk.

    :param fig: Figure to close; the current figure if None.
    :type fig: matplotlib.figure.Figure, optional
    :param show: Display the figure with :func:`matplotlib.pyplot.show` instead of closing it.
    :type show: bool, optional
    """
    if show:
        plt.show()
    else:
        plt.close(fig if fig is not None else plt.gcf())
