"""Reviewer 3, comments 4 and 6, and Reviewer 2, comment 1: the numbers behind the answers.

A. Mean irradiance and homogeneity of every configuration for the three methods, computed with the
   repository functions and conventions used for Figures 3, 7, 8 and S17 (validation of the chain).
B. Ray tracing: Monte-Carlo noise of the node values and its effect on H. The published H uses
   all mesh nodes (~0.71 mm pitch); node noise adds variance and lowers H. The noise variance is
   estimated from first differences of neighbouring nodes (var(dE) = 2 sigma^2 for uncorrelated
   noise, fitted as sigma^2 = k * E as expected for ray counting) and cross-checked by block
   averaging.
C. Radiometry: the absolute scale is fixed by assuming that the reference scan contains
   N x P_LED. The fraction of the emitted power that actually falls on the scanned area is
   estimated (model for the reflector-free reference, ray tracing for the reflector reference);
   1/f is the implied bias of all radiometric mean irradiances of that group.
D. Geometric model: sensitivity of E_mean and H to height (+-1 mm), to the LED power tolerance
   stated in the datasheets, and to an LED-to-LED spread of the radiant power.
E. Reflector gains (R3-6) and the model-vs-radiometry deviations quoted in the abstract (R2-1).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from r2_common import GRID_KW, RESOLUTION, common_args, import_repo, repo_path

PREFIX = {"UV": "UV", "GRN": "GRN"}
LED_NAME = {"UV": "LST1-01G01-UV01-00", "GRN": "LST1-01F06-GRN1-00"}
DATASHEET_TOL = {"UV": 0.06, "GRN": 0.07}   # Luminus: +/-6 % flux; Osram: +/-7 % brightness
N_RAYS = 1.6e7
PLANE_AREA_M2 = 0.33 * 0.34


def cases(repo: Path, key: str):
    led = LED_NAME[key]
    out = [("best manual", "best_manual", repo / "results/sampled_configs/best_manual.json", 8, "best_manual_no_reflector")]
    for k in (1, 2, 3):
        out.append((f"{PREFIX[key]}{k}", f"{PREFIX[key]}{k}",
                    repo / f"results/sampled_configs/near_optimal_configs/{led}/{PREFIX[key]}{k}.json", 16,
                    f"{PREFIX[key]}3_reflector"))
    return out


def read_json(path: Path):
    d = json.loads(path.read_text())
    return [tuple(p) for p in d["(x, y) Positions in [cm]"]], float(d["Height in [cm]"])


def rt_grid(df: pd.DataFrame) -> np.ndarray:
    return df.pivot(index="y", columns="x", values="Irradiance").to_numpy()


def rt_noise(E: np.ndarray) -> dict:
    """Node-noise variance model sigma^2 = k E from first differences along x and y."""
    d = np.concatenate([(E[:, 1:] - E[:, :-1]).ravel(), (E[1:, :] - E[:-1, :]).ravel()])
    m = np.concatenate([((E[:, 1:] + E[:, :-1]) / 2).ravel(), ((E[1:, :] + E[:-1, :]) / 2).ravel()])
    ok = np.isfinite(d) & np.isfinite(m) & (m > 0)
    d, m = d[ok], m[ok]
    k = np.sum(d ** 2 / 2) / np.sum(m)                      # least squares through origin of d^2/2 = k m
    vals = E[np.isfinite(E)]
    mean, var = vals.mean(), vals.var()
    var_true = var - k * mean
    out = dict(rt_node_noise_rel_at_mean_pct=100 * np.sqrt(k / mean), rt_H_raw=1 - np.sqrt(var) / mean,
               rt_H_noise_corrected=1 - np.sqrt(max(var_true, 0)) / mean)
    for b in (2, 4, 7):                                       # 7 nodes ~ 0.5 cm
        ny, nx = (E.shape[0] // b) * b, (E.shape[1] // b) * b
        Eb = np.nanmean(E[:ny, :nx].reshape(ny // b, b, nx // b, b), axis=(1, 3))
        out[f"rt_H_block{b}"] = 1 - np.nanstd(Eb) / np.nanmean(Eb)
    return out


def main() -> None:
    ap = common_args(__doc__)
    args = ap.parse_args()
    repo = repo_path(args.repo)
    GeometricModel, Grid, _, leds = import_repo(repo)
    from utils.process_radiometry_data import calculate_correction_factor, process_radiometry_file, read_radiometry_file
    from utils.process_raytracing_data import calc_metrics, process_raytracing_file

    rows, sens, norm = [], [], []
    rng = np.random.default_rng(2026)
    for key in ("UV", "GRN"):
        led = leds[key]
        name = LED_NAME[key]
        model = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        for label, stem, jpath, n_led, ref_stem in cases(repo, key):
            pos, z = read_json(jpath)
            Eg = model.simulate(pos, z)
            cf = calculate_correction_factor(str(repo / f"results/radiometry/{name}/{ref_stem}.csv"), n_led * led.total_power)
            for refl in (False, True):
                suffix = "reflector" if refl else "no_reflector"
                rt = process_raytracing_file(f"{stem}_{suffix}.txt", str(repo / f"results/raytracing/{name}"))
                E_rt, H_rt = calc_metrics(rt["Irradiance"])
                _, E_rad, H_rad = process_radiometry_file(str(repo / f"results/radiometry/{name}/{stem}_{suffix}.csv"), cf, 33, 34)
                row = dict(led=key, config=label, z_cm=z, n_leds=n_led, reflector=refl,
                           model_E=np.nan if refl else Eg.mean(), model_H=np.nan if refl else 1 - Eg.std() / Eg.mean(),
                           rt_E=E_rt, rt_H=H_rt, rad_E=E_rad, rad_H=H_rad, rad_reference=ref_stem)
                row.update(rt_noise(rt_grid(rt)))
                row["rt_E_mean_stat_rel_pct"] = 100 / np.sqrt(N_RAYS * E_rt * PLANE_AREA_M2 / (n_led * led.total_power))
                rows.append(row)
            # D. model sensitivity (reflector-free)
            m0, h0 = Eg.mean(), 1 - Eg.std() / Eg.mean()
            for dz in (-0.1, 0.1):
                E1 = model.simulate(pos, z + dz)
                sens.append(dict(led=key, config=label, z_cm=z, test="height", variant=f"dz = {dz * 10:+.0f} mm",
                                 dE_rel_pct=100 * (E1.mean() / m0 - 1), dH=(1 - E1.std() / E1.mean()) - h0))
            sens.append(dict(led=key, config=label, z_cm=z, test="LED power tolerance (datasheet)",
                             variant=f"+-{100 * DATASHEET_TOL[key]:.0f} % common", dE_rel_pct=100 * DATASHEET_TOL[key], dH=0.0))
            singles = np.stack([model.simulate([p], z) for p in pos])
            for sigma in (0.05, 0.10):
                hs, es = [], []
                for _ in range(2000):
                    a = rng.normal(1.0, sigma, len(pos))
                    E1 = np.tensordot(a, singles, axes=1)
                    hs.append(1 - E1.std() / E1.mean()); es.append(E1.mean())
                hs = np.asarray(hs) - h0
                sens.append(dict(led=key, config=label, z_cm=z, test="LED-to-LED power spread",
                                 variant=f"sigma = {100 * sigma:.0f} %", dE_rel_pct=100 * np.std(np.asarray(es) / m0),
                                 dH=float(np.mean(hs)), dH_p2_5=float(np.percentile(hs, 2.5)), dH_p97_5=float(np.percentile(hs, 97.5))))
        # C. normalisation of the radiometry
        for ref_label, stem, jpath, n_led, has_refl in (("best manual (no reflector)", "best_manual", repo / "results/sampled_configs/best_manual.json", 8, False),
                                                         (f"{PREFIX[key]}3 (reflector)", f"{PREFIX[key]}3", repo / f"results/sampled_configs/near_optimal_configs/{name}/{PREFIX[key]}3.json", 16, True)):
            pos, z = read_json(jpath)
            df = read_radiometry_file(str(repo / f"results/radiometry/{name}/{stem}_{'reflector' if has_refl else 'no_reflector'}.csv"))
            span = float(df.columns.to_numpy()[-1] - df.columns.to_numpy()[0]) + float(df.columns.to_numpy()[1] - df.columns.to_numpy()[0])
            if not has_refl:
                for zs in (z, z + 1.0):
                    big = GeometricModel(Grid(width=span, height=span, step=2.5, side_space=1.5, top_bottom_space=2.0), led, resolution_xy=0.1)
                    Eb = big.simulate(pos, zs)
                    frac = Eb.sum() * (0.001 ** 2) / (n_led * led.total_power)
                    norm.append(dict(led=key, reference=ref_label, scan_side_cm=span, assumed_distance_cm=zs,
                                     method="geometric model, flux on scanned square", captured_fraction=frac,
                                     implied_radiometry_bias_pct=100 * (1 / frac - 1)))
            else:
                rt = process_raytracing_file(f"{stem}_reflector.txt", str(repo / f"results/raytracing/{name}"))
                E_rt, _ = calc_metrics(rt["Irradiance"])
                frac = E_rt * PLANE_AREA_M2 / (n_led * led.total_power)
                norm.append(dict(led=key, reference=ref_label, scan_side_cm=span, assumed_distance_cm=z,
                                 method="ray tracing, flux on 33 x 34 cm floor inside the reflector", captured_fraction=frac,
                                 implied_radiometry_bias_pct=100 * (1 / frac - 1)))
    df = pd.DataFrame(rows)
    # E. reflector gains and deviations
    gains = []
    for (key, label), g in df.groupby(["led", "config"], sort=False):
        a, b = g[~g.reflector].iloc[0], g[g.reflector].iloc[0]
        gains.append(dict(led=key, config=label, z_cm=a.z_cm,
                          gain_rad_pct=100 * (b.rad_E / a.rad_E - 1), gain_rt_pct=100 * (b.rt_E / a.rt_E - 1),
                          H_rad_no=a.rad_H, H_rad_refl=b.rad_H, H_rt_no=a.rt_H, H_rt_refl=b.rt_H,
                          E_rad_no=a.rad_E, E_rad_refl=b.rad_E,
                          model_vs_rad_E_rel_pct=100 * (a.model_E / a.rad_E - 1),
                          model_vs_rad_H_rel_pct=100 * (a.model_H / a.rad_H - 1), model_vs_rad_H_abs=a.model_H - a.rad_H,
                          rt_vs_rad_E_rel_pct=100 * (a.rt_E / a.rad_E - 1), rt_vs_rad_H_abs=a.rt_H - a.rad_H,
                          rt_refl_vs_rad_refl_E_rel_pct=100 * (b.rt_E / b.rad_E - 1), rt_refl_vs_rad_refl_H_abs=b.rt_H - b.rad_H))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "r3_4_metrics_all_methods.csv", index=False, float_format="%.6g")
    pd.DataFrame(sens).to_csv(out / "r3_4_model_sensitivity.csv", index=False, float_format="%.6g")
    pd.DataFrame(norm).to_csv(out / "r3_4_radiometry_normalisation.csv", index=False, float_format="%.6g")
    pd.DataFrame(gains).to_csv(out / "r3_6_reflector_gains_and_deviations.csv", index=False, float_format="%.6g")
    with pd.option_context("display.width", 300, "display.max_columns", 40, "display.max_rows", 100):
        for name_, t in (("A/B metrics", df), ("D sensitivity", pd.DataFrame(sens)), ("C normalisation", pd.DataFrame(norm)), ("E gains", pd.DataFrame(gains))):
            print(f"\n=== {name_}")
            print(t.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
