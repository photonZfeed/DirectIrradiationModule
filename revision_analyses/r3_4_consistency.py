"""Reviewer 3, comment 4: two systematic effects found while quantifying the uncertainties.

1. Ray tracing homogeneity on the common 0.5 cm grid. The published ray tracing H is computed from
   all mesh nodes (0.71 mm pitch), whose Monte-Carlo noise lowers H. Block-averaging onto the
   0.5 cm grid of Supporting Information S6 (``utils.compare_methods.field_raytracing``) removes
   most of it. The same block sizes applied to the noise-free geometric model show how much of the
   change is noise and how much is the averaging itself.
2. Absolute scale of the radiometry. ``calculate_correction_factor`` scales a reference scan so that
   it contains N x P_LED. For the reflector-free 8-LED reference, the geometric model gives the
   fraction of the emitted flux that falls on the scanned square at the scan distance (nominal
   height + 1 cm safety offset). The share of the scan power inside the 33 cm x 34 cm crop is
   compared between model and measurement as a check of the angular tails.
3. The same rescaling for the 16-LED group, using the reflector-free X3 scan as reference
   (``rescaled_group_values``), and the model homogeneity at the scan distance
   (``plane_matched_homogeneity``).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from r2_common import GRID_KW, RESOLUTION, common_args, import_repo, repo_path

LED_NAME = {"UV": "LST1-01G01-UV01-00", "GRN": "LST1-01F06-GRN1-00"}


def H(a):
    a = a[np.isfinite(a)]
    return 1 - a.std() / a.mean()


def main() -> None:
    ap = common_args(__doc__)
    args = ap.parse_args()
    repo = repo_path(args.repo)
    GeometricModel, Grid, _, leds = import_repo(repo)
    from geometric_model import calculate_irradiance
    from utils.compare_methods import common_grid, field_raytracing
    from utils.process_raytracing_data import calc_metrics, process_raytracing_file
    from utils.process_radiometry_data import calculate_correction_factor, process_radiometry_file, read_radiometry_file

    def field(X, Y, pos, z, led):
        X = np.ascontiguousarray(X, dtype=np.float64)
        Y = np.ascontiguousarray(Y, dtype=np.float64)
        return sum(calculate_irradiance(X - p[0], Y - p[1], float(z), led.angle_data, led.intensity_data) for p in pos) * 1e4

    xg, yg = common_grid()
    rows, scale = [], []
    for key, name in LED_NAME.items():
        led = leds[key]
        GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)       # normalises the LED profile
        prefix = "UV" if key == "UV" else "GRN"
        stems = [("best_manual", repo / "results/sampled_configs/best_manual.json")] + [
            (f"{prefix}{k}", repo / f"results/sampled_configs/near_optimal_configs/{name}/{prefix}{k}.json") for k in (1, 2, 3)]
        for stem, jpath in stems:
            d = json.loads(jpath.read_text())
            pos = [tuple(p) for p in d["(x, y) Positions in [cm]"]]
            z = float(d["Height in [cm]"])
            for refl in ("no_reflector", "reflector"):
                df = process_raytracing_file(f"{stem}_{refl}.txt", str(repo / f"results/raytracing/{name}"))
                _, h_nodes = calc_metrics(df["Irradiance"])
                h_grid = H(field_raytracing(name, f"{stem}_{refl}", xg, yg, results_path=str(repo / "results")))
                row = dict(led=key, config=stem, z_cm=z, reflector=refl == "reflector",
                           rt_H_nodes_published=h_nodes, rt_H_common_05cm=h_grid, change=h_grid - h_nodes)
                if refl == "no_reflector":
                    P = df.pivot(index="y", columns="x", values="Irradiance")
                    X, Y = np.meshgrid(P.columns.to_numpy(float), P.index.to_numpy(float))
                    Em = field(X, Y, pos, z, led)
                    Xg, Yg = np.meshgrid(xg, yg)
                    # same block averaging applied to the noise-free model sampled at the RT nodes
                    from utils.compare_methods import regrid_averaged
                    hm_nodes = H(Em)
                    hm_grid = H(regrid_averaged(P.columns.to_numpy(float), P.index.to_numpy(float), Em, xg, yg))
                    row.update(model_H_at_rt_nodes=hm_nodes, model_H_same_block_average=hm_grid,
                               change_explained_by_averaging=hm_grid - hm_nodes,
                               estimated_rt_noise_bias_in_H=(h_nodes - h_grid) - (hm_nodes - hm_grid))
                rows.append(row)
        # absolute scale of the reflector-free reference scan (best manual, 8 LEDs)
        d = json.loads((repo / "results/sampled_configs/best_manual.json").read_text())
        pos = [tuple(p) for p in d["(x, y) Positions in [cm]"]]
        scan = read_radiometry_file(str(repo / f"results/radiometry/{name}/best_manual_no_reflector.csv"))
        c, i = scan.columns.to_numpy(float), scan.index.to_numpy(float)
        A = scan.to_numpy(float)
        Xc, Yc = np.meshgrid(c - c.mean(), i - i.mean())
        crop = (np.abs(Xc) <= 17.0) & (np.abs(Yc) <= 16.5)
        pitch = c[1] - c[0]
        cf = calculate_correction_factor(str(repo / f"results/radiometry/{name}/best_manual_no_reflector.csv"), 8 * led.total_power)
        _, e_rad, _ = process_radiometry_file(str(repo / f"results/radiometry/{name}/best_manual_no_reflector.csv"), cf, 33, 34)
        model_plane = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        for zs in (13.0, 14.0):
            fine = np.arange(-(c.max() - c.min() + pitch) / 2 + 0.025, (c.max() - c.min() + pitch) / 2, 0.05)
            Xf, Yf = np.meshgrid(fine, fine)
            Ef = field(Xf, Yf, pos, zs, led)
            captured = Ef.sum() * (0.0005 ** 2) / (8 * led.total_power)
            Es = field(Xc, Yc, pos, zs, led)
            e_model = model_plane.simulate(pos, zs).mean()
            scale.append(dict(led=key, distance_cm=zs, scan_side_cm=c.max() - c.min() + pitch,
                              crop_share_of_scan_measured=A[crop].sum() / A.sum(),
                              crop_share_of_scan_model=Es[crop].sum() / Es.sum(),
                              captured_fraction_model=captured,
                              rad_E_mean_published=e_rad, rad_E_mean_rescaled=e_rad * captured,
                              model_E_mean_at_this_distance=e_model,
                              model_vs_rescaled_rad_pct=100 * (e_model / (e_rad * captured) - 1)))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    a, b = pd.DataFrame(rows), pd.DataFrame(scale)
    a.to_csv(out / "r3_4_rt_homogeneity_grid.csv", index=False, float_format="%.6g")
    b.to_csv(out / "r3_4_radiometry_scale.csv", index=False, float_format="%.6g")
    with pd.option_context("display.width", 300, "display.max_columns", 40):
        print(a.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
        print(b.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


def rescaled_group_values(repo: Path) -> pd.DataFrame:
    """Alternative normalisation for the 16-LED group: reflector-free X3 scan + model capture fraction.

    The scan of X3 without reflector (nominal 6 cm, scanned at 7 cm) is scaled so that it contains
    the share of 16 x P_LED that the geometric model puts on the scanned square at 7 cm. The same
    factor is then applied to all scans of the group, as in the published convention.
    """
    GeometricModel, Grid, _, leds = import_repo(repo)
    from geometric_model import calculate_irradiance
    from utils.process_radiometry_data import calculate_correction_factor, process_radiometry_file, read_radiometry_file
    rows = []
    for key, name in LED_NAME.items():
        led = leds[key]
        model = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        prefix = "UV" if key == "UV" else "GRN"
        d = json.loads((repo / f"results/sampled_configs/near_optimal_configs/{name}/{prefix}3.json").read_text())
        pos3 = [tuple(p) for p in d["(x, y) Positions in [cm]"]]
        ref = str(repo / f"results/radiometry/{name}/{prefix}3_no_reflector.csv")
        scan = read_radiometry_file(ref)
        c = scan.columns.to_numpy(float)
        side = c.max() - c.min() + (c[1] - c[0])
        fine = np.arange(-side / 2 + 0.025, side / 2, 0.05)
        Xf, Yf = np.meshgrid(fine, fine)
        Ef = sum(calculate_irradiance(np.ascontiguousarray(Xf - p[0]), np.ascontiguousarray(Yf - p[1]), 7.0,
                                      led.angle_data, led.intensity_data) for p in pos3) * 1e4
        captured = Ef.sum() * 0.0005 ** 2 / (16 * led.total_power)
        cf_new = calculate_correction_factor(ref, 16 * led.total_power * captured)
        cf_pub = calculate_correction_factor(str(repo / f"results/radiometry/{name}/{prefix}3_reflector.csv"), 16 * led.total_power)
        for k in (1, 2, 3):
            dk = json.loads((repo / f"results/sampled_configs/near_optimal_configs/{name}/{prefix}{k}.json").read_text())
            posk = [tuple(p) for p in dk["(x, y) Positions in [cm]"]]
            zk = float(dk["Height in [cm]"])
            for refl in ("no_reflector", "reflector"):
                f = str(repo / f"results/radiometry/{name}/{prefix}{k}_{refl}.csv")
                _, e_pub, h = process_radiometry_file(f, cf_pub, 33, 34)
                _, e_new, _ = process_radiometry_file(f, cf_new, 33, 34)
                row = dict(led=key, config=f"{prefix}{k}", z_cm=zk, reflector=refl == "reflector",
                           captured_fraction_X3_no_reflector_at_7cm=captured, scale_new_over_published=cf_new / cf_pub,
                           rad_E_published=e_pub, rad_E_rescaled=e_new, rad_H=h)
                if refl == "no_reflector":
                    row["model_E_at_z_plus_1cm"] = model.simulate(posk, zk + 1).mean()
                    row["model_vs_rescaled_pct"] = 100 * (row["model_E_at_z_plus_1cm"] / e_new - 1)
                rows.append(row)
        scan3 = read_radiometry_file(str(repo / f"results/radiometry/{name}/{prefix}3_reflector.csv"))
        c3 = scan3.columns.to_numpy(float)
        pix = ((c3[1] - c3[0]) / 100) ** 2
        rows.append(dict(led=key, config=f"{prefix}3 reflector: share of 16 x P_LED in the scan after rescaling",
                         rad_E_rescaled=scan3.to_numpy(float).sum() * cf_new * pix / (16 * led.total_power)))
    return pd.DataFrame(rows)


def plane_matched_homogeneity(repo: Path) -> pd.DataFrame:
    """Model homogeneity at the nominal height and at the scan distance (+1 cm) vs radiometry."""
    GeometricModel, Grid, _, leds = import_repo(repo)
    from utils.process_radiometry_data import calculate_correction_factor, process_radiometry_file
    rows = []
    for key, name in LED_NAME.items():
        led = leds[key]
        model = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        prefix = "UV" if key == "UV" else "GRN"
        cases = [("best manual", repo / "results/sampled_configs/best_manual.json", "best_manual")] + [
            (f"{prefix}{k}", repo / f"results/sampled_configs/near_optimal_configs/{name}/{prefix}{k}.json", f"{prefix}{k}")
            for k in (1, 2, 3)]
        for label, jpath, stem in cases:
            d = json.loads(jpath.read_text())
            pos = [tuple(p) for p in d["(x, y) Positions in [cm]"]]
            z = float(d["Height in [cm]"])
            _, _, h_rad = process_radiometry_file(str(repo / f"results/radiometry/{name}/{stem}_no_reflector.csv"), 1.0, 33, 34)
            E0, E1 = model.simulate(pos, z), model.simulate(pos, z + 1)
            h0, h1 = 1 - E0.std() / E0.mean(), 1 - E1.std() / E1.mean()
            rows.append(dict(led=key, config=label, z=z, H_rad=h_rad, H_model_z=h0, H_model_z1=h1,
                             dH_z=h0 - h_rad, dH_z1=h1 - h_rad, rel_z1_pct=100 * (h1 / h_rad - 1)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
    _args = common_args(__doc__).parse_args()
    _repo, _out = repo_path(_args.repo), Path(_args.out)
    for _name, _func in (("r3_4_radiometry_rescaled_16led_group.csv", rescaled_group_values),
                         ("r3_4_model_H_plane_matched.csv", plane_matched_homogeneity)):
        _df = _func(_repo)
        _df.to_csv(_out / _name, index=False, float_format="%.6g")
        with pd.option_context("display.width", 300, "display.max_columns", 40):
            print(_df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
