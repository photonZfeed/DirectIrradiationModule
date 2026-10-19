"""Reviewer 3, comment 1: error of the point-source assumption and of the emitter reference plane.

Two sensitivity calculations for representative configurations, all evaluated on the model's
0.5 cm grid over the 33 cm x 34 cm detection plane:

(a) Finite emitter. Every LED is replaced by a uniformly emitting planar disk of diameter d whose
    surface elements all radiate with the tabulated relative profile I_rel(theta) (the way a Speos
    surface source with an intensity file is set up). d = 2.6 mm (UV) and 2.15 mm (green) are the
    values used in the ray tracing; 5 mm is added as a deliberately pessimistic bound.
(b) Emitter-plane offset. The emitting surface sits above the mounting plane that defines z;
    the field is recomputed for z + dz, dz = -1, -2, -3 mm (emitter closer to the plane).

Outputs: <out>/r3_1_source_model.csv and a printed table.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from r2_common import GRID_KW, RESOLUTION, common_args, import_repo, repo_path

CONFIGS = {
    "UV": [("best manual", "results/sampled_configs/best_manual.json", None),
           ("UV1", "results/sampled_configs/near_optimal_configs/LST1-01G01-UV01-00/UV1.json", None),
           ("UV3", "results/sampled_configs/near_optimal_configs/LST1-01G01-UV01-00/UV3.json", None),
           ("UV3 at 5 cm", "results/sampled_configs/near_optimal_configs/LST1-01G01-UV01-00/UV3.json", 5.0)],
    "GRN": [("best manual", "results/sampled_configs/best_manual.json", None),
            ("GRN3", "results/sampled_configs/near_optimal_configs/LST1-01F06-GRN1-00/GRN3.json", None),
            ("GRN3 at 5 cm", "results/sampled_configs/near_optimal_configs/LST1-01F06-GRN1-00/GRN3.json", 5.0)],
}
SPEOS_DIAMETER_MM = {"UV": 2.6, "GRN": 2.15}


def read_json(path: Path):
    d = json.loads(path.read_text())
    return [tuple(p) for p in d["(x, y) Positions in [cm]"]], float(d["Height in [cm]"])


def disk_nodes(d_cm: float, n_r: int = 12, n_phi: int = 24):
    """Gauss-Legendre (radial) x uniform (azimuthal) quadrature on a disk; weights sum to 1."""
    if d_cm == 0:
        return np.zeros(1), np.zeros(1), np.ones(1)
    a = d_cm / 2
    xg, wg = np.polynomial.legendre.leggauss(n_r)
    r = a * np.sqrt((xg + 1) / 2)                   # area-uniform mapping
    wr = wg / 2
    phi = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    R, P = np.meshgrid(r, phi, indexing="ij")
    W = np.repeat(wr[:, None], n_phi, axis=1) / n_phi
    return (R * np.cos(P)).ravel(), (R * np.sin(P)).ravel(), W.ravel()


def field(X, Y, positions, z, angle_deg, intensity, d_cm=0.0):
    """Irradiance in W m-2 (positions and z in cm, intensity in W sr-1)."""
    ex, ey, ew = disk_nodes(d_cm)
    E = np.zeros_like(X)
    for (px, py) in positions:
        for dx, dy, w in zip(ex, ey, ew):
            xs, ys = X - (px + dx), Y - (py + dy)
            rho2 = xs ** 2 + ys ** 2
            theta = np.arctan(np.sqrt(rho2) / z)
            I = np.interp(np.degrees(theta), angle_deg, intensity)
            E += w * z / (rho2 + z ** 2) ** 1.5 * I
    return E * 1e4


def metrics(E):
    m = E.mean()
    return m, 1 - E.std() / m, E.max()


def main() -> None:
    ap = common_args(__doc__)
    args = ap.parse_args()
    repo = repo_path(args.repo)
    GeometricModel, Grid, _, leds = import_repo(repo)
    rows = []
    for key, cases in CONFIGS.items():
        led = leds[key]
        model = GeometricModel(Grid(**GRID_KW), led, resolution_xy=RESOLUTION)
        X, Y = model.X, model.Y
        ang, inten = led.angle_data, led.intensity_data
        for label, rel, z_override in cases:
            pos, z = read_json(repo / rel)
            z = z_override or z
            E0 = field(X, Y, pos, z, ang, inten)
            check = model.simulate(pos, z)
            assert np.allclose(E0, check, rtol=1e-9, atol=1e-9), "point-source baseline differs from the repository model"
            m0, h0, x0 = metrics(E0)
            variants = [("finite disk", f"d = {SPEOS_DIAMETER_MM[key]} mm (Speos)", dict(d_cm=SPEOS_DIAMETER_MM[key] / 10)),
                        ("finite disk", "d = 5 mm (pessimistic)", dict(d_cm=0.5))]
            for label_v, desc, kw in variants:
                E1 = field(X, Y, pos, z, ang, inten, **kw)
                m1, h1, x1 = metrics(E1)
                rows.append(dict(led=key, config=label, z_cm=z, test=label_v, variant=desc,
                                 E_mean_point=m0, H_point=h0,
                                 dE_mean_rel_pct=100 * (m1 - m0) / m0, dH_abs=h1 - h0,
                                 max_local_dev_pct=100 * np.max(np.abs(E1 - E0)) / m0,
                                 rms_local_dev_pct=100 * np.sqrt(np.mean((E1 - E0) ** 2)) / m0))
            for dz_mm in (-1, -2, -3):
                E1 = field(X, Y, pos, z + dz_mm / 10, ang, inten)
                m1, h1, x1 = metrics(E1)
                rows.append(dict(led=key, config=label, z_cm=z, test="emitter-plane offset", variant=f"dz = {dz_mm} mm",
                                 E_mean_point=m0, H_point=h0,
                                 dE_mean_rel_pct=100 * (m1 - m0) / m0, dH_abs=h1 - h0,
                                 max_local_dev_pct=100 * np.max(np.abs(E1 - E0)) / m0,
                                 rms_local_dev_pct=100 * np.sqrt(np.mean((E1 - E0) ** 2)) / m0))
    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "r3_1_source_model.csv", index=False, float_format="%.6g")
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(df.drop(columns=["E_mean_point", "H_point"]).to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
