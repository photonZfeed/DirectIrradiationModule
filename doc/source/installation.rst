Getting Started: Installation
==============================

The geometric radiation field model, the sampling and optimization routines and the evaluation code
for the radiometry and ray-tracing data form a single Python project. This page describes how to
obtain it and how to set up a working environment.

Requirements
------------

Python **3.12 or newer** is required. The published results were produced with **Python 3.13.3**. The code
runs on Windows, Linux and macOS.

All dependencies are pinned, in ``requirements.txt`` and in ``pyproject.toml``, to the versions of
that reference environment, so a fresh installation reproduces the published figures and statistics.
Both files are kept identical; ``doc/requirements.txt`` holds the two additional packages needed to
build this documentation.

Obtaining the code
------------------

.. code-block:: bash

    git clone https://github.com/photonZfeed/DirectIrradiationModule
    cd DirectIrradiationModule

Setting up a virtual environment
--------------------------------

Working in a virtual environment keeps the pinned versions separate from the rest of the system.

Linux and macOS
^^^^^^^^^^^^^^^

.. code-block:: bash

    python3 -m venv .venv
    source .venv/bin/activate

Windows
^^^^^^^

.. code-block:: doscon

    py -m venv .venv
    .venv\Scripts\activate

Installing the dependencies
---------------------------

.. code-block:: bash

    pip install -r requirements.txt

This is the recommended route for reproducing the results of the publication. The scripts
``usage_example.py`` and ``create_manuscript_plots_tables.py`` resolve the ``results/``, ``figures/``
and ``tables/`` directories relative to the root of the clone and change to it themselves, so they
can be started from any folder (see :doc:`reproduce`). Modules imported in an interactive session or
in a script of your own also use paths relative to the root of the clone; start such sessions there.

Installing the project as a package
-----------------------------------

If the modules are to be imported from outside the repository, the project can also be installed
from its ``pyproject.toml``:

.. code-block:: bash

    pip install .      # regular installation
    pip install -e .   # editable installation for development

This installs the top-level modules :mod:`geometric_model` and :mod:`optimizer` together with the
``utils`` (including the samplers :mod:`utils.systematic_sampler` and :mod:`utils.monte_carlo_sampler`),
``analysis`` and ``visualization`` packages. The two runner scripts above and the tests are deliberately
**not** installed, and the measurement and simulation data under ``results/`` are not part of the
distribution, so reproducing the publication figures still requires a clone. The optional extras
``docs`` and ``test`` install the packages for building the documentation and for running the tests,
e.g. ``pip install -e ".[test]"``.

Building the documentation
--------------------------

.. code-block:: bash

    pip install -r doc/requirements.txt
    cd doc
    make html

On Windows use ``make.bat html`` instead. The rendered pages are written to ``doc/build/html``.
The same build runs automatically on every push to ``main`` and is published to GitHub Pages.

Checking the installation
-------------------------

.. code-block:: bash

    python usage_example.py

The script walks through the full workflow — grid and LED definition, candidate generation,
simulation with the geometric model, optimization and the evaluation of radiometry and ray-tracing
data — and is described step by step in :doc:`tutorial`. It displays its figures one after the other
and writes its results to ``results/examples/``; the example files stored there in the repository are
overwritten.

The tests check the fast evaluation of many configurations against the geometric model and, if the
Zenodo data are present, the reproducibility of the published Monte-Carlo sample (see :doc:`reproduce`):

.. code-block:: bash

    pip install pytest
    python -m pytest tests
