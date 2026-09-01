Getting Started: Installation
==============================

The geometric radiation field model, the sampling and optimization routines and the evaluation code
for the radiometry and raytracing data form a single Python project. This page describes how to
obtain it and how to set up a working environment.

Requirements
------------

Python **3.12 or newer** is required. The published results were produced with **Python 3.13.3**.

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

    python3 -m venv venv
    source venv/bin/activate

Windows
^^^^^^^

.. code-block:: doscon

    py -m venv venv
    venv\Scripts\activate

Installing the dependencies
---------------------------

.. code-block:: bash

    pip install -r requirements.txt

This is the recommended route for reproducing the results of the publication: the scripts
``usage_example.py`` and ``create_manuscript_plots.py`` resolve the ``results/`` and ``figures/``
directories relative to the current working directory, so they are run from the root of the clone.

Installing the project as a package
-----------------------------------

If the modules are to be imported from outside the repository, the project can also be installed
from its ``pyproject.toml``:

.. code-block:: bash

    pip install .      # regular installation
    pip install -e .   # editable installation for development

This installs the top-level modules :mod:`geometric_model`, :mod:`systematic_sampler`,
:mod:`monte_carlo_sampler` and :mod:`optimizer` together with the ``utils`` and
``visualization`` packages. The two runner scripts above are deliberately **not** installed, and
the measurement and simulation data under ``results/`` are not part of the distribution, so
reproducing the publication figures still requires a clone.

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
simulation with the geometric model, optimization and the evaluation of radiometry and raytracing
data — and is described step by step in :doc:`tutorial`.
