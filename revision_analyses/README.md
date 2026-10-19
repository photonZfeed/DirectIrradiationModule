# Analyses of the second review round

These scripts produce the additional results requested in the second review round: the
sensitivity of the geometric model to the source assumptions, the influence of the irradiance
bound, the convergence of the Monte-Carlo optimization and the uncertainties of the three methods.
They import the code of this repository (`geometric_model`, `monte_carlo_sampler`, `optimizer`,
`utils`) and do not change it.

All scripts are run from the repository root with the environment of `requirements.txt`. By
default they use the repository that contains this folder (`--repo` or `DIM_REPO` override it) and
write to `results/revision_analyses/`, which is not tracked by git (`--out` overrides it).

```
python revision_analyses/r3_1_source_model.py
```

## Scripts

| Script | Purpose | Reviewer point | Output in `results/revision_analyses/` | Supporting Information |
|---|---|---|---|---|
| `r2_common.py` | Shared helpers: import of the repository, exact fast evaluation of many configurations (single-LED field bank, Gram matrix for the standard deviation, half-grid product for the maximum), Pareto filter and hypervolume | – | – | – |
| `validate_fast_eval.py` | Checks the fast evaluation against the stored results of the systematic 8-LED study (optionally `--mc-uv` for the 5 million-sample UV file) | – | printed | – |
| `extract_mc_configs.py` | Reads the LED positions of a published 5 million-sample file (simulation results or sampled set) without unpickling it (more than 5 GB otherwise) | R3-2, R3-3 | `mc_published_seed42_5000000_configs.npz` | – |
| `r3_1_source_model.py` | Finite emitting disk and emitter-plane offset against the point-source model | R3-1 | `r3_1_source_model.csv` | S9, Table S3 |
| `r3_mc_runs.py` | Sampling with `MonteCarloSampler` (sequential, seeded, duplicates removed, i.e. `generate_unique_samples`; with seed 42 the first 5 million configurations are the published sample) and evaluation of the sampled configurations | R3-2, R3-3 | `mc_*_configs.npz`, `mc_*_metrics_{UV,GRN}.npz` | – |
| `check_seed42_published.py` | Compares seed 42 with the published sample (row by row) and with both published simulation files (as a set, and their stored mean, standard deviation and maximum) | R3-3 | printed | S11 |
| `run_all_seeds.sh` | Runs `r3_mc_runs.py` for seeds 42, 1, 2 and 3 with 10 million distinct candidates each; finished stages are skipped | R3-3 | as above | – |
| `r3_2_ebound.py` | Selection of the manuscript repeated for irradiance bounds of 100–300 W m⁻² and without a bound | R3-2 | `r3_2_ebound_sensitivity.csv`, `r3_2_nondominated_fronts.csv` | S10, Table S4 |
| `r3_2_ebound_8led.py` | Highest local irradiance of all 482,931 symmetric 8-LED configurations at 13 cm | R3-2 | `r3_2_ebound_8led.csv` | S10 |
| `r3_3_convergence.py` | Hypervolume, additive ε indicator and selected designs per seed (42, 1, 2, 3) and sample size (0.5–10 million), and the best configuration with its LED positions; with `--local-search` also a local search around the best candidates (not reported in the manuscript) | R3-3 | `r3_3_convergence.csv`, `r3_3_best_configurations.csv`, `r3_3_reference_fronts.npz` (`r3_3_local_search_{UV,GRN}.csv`) | S11 |
| `r3_3_figure.py` | Pareto fronts per seed and best score versus sample size; best configuration of every seed at 5 million candidates; drawn in the style of manuscript Figures 1a, 4 and 5 | R3-3 | `r3_3_convergence.svg/.png`, `r3_3_best_configurations.svg/.png` | Figures S18, S19 |
| `r3_4_uncertainty.py` | Metrics of all three methods, ray-tracing node noise, normalization of the radiometry, sensitivity of the model (height, LED flux tolerance and spread), reflector gains and model–radiometry deviations | R3-4, R3-6, R2-1 | `r3_4_metrics_all_methods.csv`, `r3_4_model_sensitivity.csv`, `r3_4_radiometry_normalisation.csv`, `r3_6_reflector_gains_and_deviations.csv` | S12, Table S5 |
| `r3_4_consistency.py` | Ray-tracing homogeneity on the 0.5 cm grid, absolute scale of the radiometry, rescaling of the 16-LED group, model homogeneity at the scan distance | R3-4 | `r3_4_rt_homogeneity_grid.csv`, `r3_4_radiometry_scale.csv`, `r3_4_radiometry_rescaled_16led_group.csv`, `r3_4_model_H_plane_matched.csv` | S12 |

`r3_1_source_model.py`, `r3_2_ebound_8led.py`, `r3_4_uncertainty.py`, `r3_4_consistency.py` and
`validate_fast_eval.py` need only the data in this repository and take less than a minute each.

## Monte-Carlo studies (R3-2, R3-3)

R3-2 uses the published 5 million-sample UV simulation file, which is too large for the repository and is
available from the Zenodo dataset. Place it in `results/geometric/LST1-01G01-UV01-00/`, then:

```
python revision_analyses/extract_mc_configs.py results/geometric/LST1-01G01-UV01-00/simulation_results_MonteCarlo_1-16_leds_5-15_cm_5000000_samples.npz results/revision_analyses/mc_published_seed42_5000000_configs.npz
python revision_analyses/r3_mc_runs.py evaluate --configs results/revision_analyses/mc_published_seed42_5000000_configs.npz
python revision_analyses/r3_2_ebound.py
```

The green metrics of the published set are recomputed from the same configurations (the rows of the green
simulation file are stored in a different order).

R3-3 needs only the repository: seeds 42, 1, 2 and 3 are drawn with `MonteCarloSampler`, and the first
5 million configurations of seed 42 are the published sample.

```
PY=python bash revision_analyses/run_all_seeds.sh results/revision_analyses
python revision_analyses/r3_3_convergence.py
python revision_analyses/r3_3_figure.py
```

`run_all_seeds.sh` takes about two hours on two cores and writes intermediate files of 135–191 MB each; a
sampler needs about 2 GB of memory at the end of a run. On Windows, run it from Git Bash or run the
`r3_mc_runs.py` commands it lists one by one. The intermediate files of all four seeds, the tables and their
SHA-256 checksums are part of the Zenodo dataset (`results/geometric/monte_carlo_convergence/`); copied to
`results/revision_analyses/`, they let `r3_3_convergence.py` and `r3_3_figure.py` run without sampling.
`check_seed42_published.py` compares seed 42 with the published files after their positions have been read
with `extract_mc_configs.py` (see its docstring for the arguments).
