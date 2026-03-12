# Direct Irradiation Module Optimization

## Overview

![Geometric Model Overview](doc/images/geometric_model_overview.png)


This project is part of the publication ["Making Photocatalyst Screening Photo-Efficient – Geometric Radiation Field Model Assisted Design of a Direct Irradiation Module for the Multi-Batch Screening Reactor"](https://doi.org/10.5281/zenodo.17867874), DOI: [10.5281/zenodo.17867874](https://doi.org/10.5281/zenodo.17867874). It contains the python code used for the simulations and optimizations described in the publication as well as the corresponding radiometry- and raytracing data. The full dataset including radiometry raw data and eletronic lab journal entries about the conducted simulations is available on Zenodo: [https://zenodo.org/records/18610997](https://zenodo.org/records/18610997).

## Features
- Systematic-sampler and Monte-Carlo-sampler to generate candidate led arrangements for direct irradiation module optimization 
- Geometric radiation field model
- Optimization algorithms to find best trade-offs between mean irradiance and homogeneity
- Raytracing- and radiometry data evaluation
- Data analysis and visualization

## Installation
Follow these steps to set up the project:

1. Download or copy the web address (URL) of this project from GitHub.
2. Open the Command Prompt (Windows) or Terminal (Mac/Linux) on your computer.
3. Type this command to copy the project to your computer:

	```bash
	git clone https://github.com/photonZfeed/DirectIrradiationModule
	```

4. Move into the project folder:

	```bash
	cd direct_irradiation_module_optimization
	```

5. Install the required Python packages:

	```bash
	pip install -r requirements.txt
	```

If you have never used `git` or `pip` before, you can find instructions online for installing [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and [Python with pip](https://realpython.com/installing-python/). The project was developed and tested with Python 3.13.3.

## Usage
Run the usage example script `usage_example.py` to see an example of how to use the modules in this script or run `create_manuscript_plots.py` to generate the plots used in the publication.

```bash
python usage_example.py
```

or

```bash
python create_manuscript_plots.py
```

Refer to the [documentation](https://photonzfeed.github.io/DirectIrradiationModule/)  and `usage_example.py` for detailed usage instructions.

## Documentation
Documentation is available at [https://photonzfeed.github.io/DirectIrradiationModule](https://photonzfeed.github.io/DirectIrradiationModule) and includes guides on how to build the direct irradiation module, how to run the radiometry measurements, how to perform the raytracing simulations and how to use the optimization code including a detailed overview on the API.

## AI Assistance Disclosure

This repository was developed with assistance of AI language models. These tools were used to support code development, documentation creation, and identification of code correctness issues. All AI-generated changes were reviewed and approved by the authors before being committed. 

The following AI tools were used:

- **Claude Code** (command-line interface by Anthropic), powered by the **Claude Sonnet 4.6** model (Anthropic, PBC).
- **GitHub Copilot** (Visual Studio Code extension by GitHub, Inc.), powered by the **OpenAI GPT-4** model (OpenAI, L.L.C.).
