#!/usr/bin/env python3
"""
Out-of-sample ANGLE validation for the weight-six side-law model.

Requires in the same folder:
    weight6_torus_full_recovery.py
    recovered_constants_high_order.csv

The recovered constants were obtained from side angles
    x = 0.4, 0.8, 1.2.

This script NEVER uses those angles.  It evaluates
    x = 0.55, 0.95, 1.35

on:
    geometry 1  : pure angle holdout (geometry was in the fit),
    geometries 8--10 : simultaneous geometry + angle holdouts.

For each observation it extracts s6 using
    [1, rho^2, rho^4, rho^6, rho^8]
from the opposite-direction averaged intrinsic geodesic data and compares
against the quadratic-in-r^2 angular model.

Dependencies:
    mpmath, pandas
"""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import mpmath as mp
import pandas as pd

import weight6_torus_full_recovery as w6


TEST_GEOMETRY_IDS = (1, 8, 9, 10)
DEFAULT_ANGLES = ("0.55", "0.95", "1.35")
DEFAULT_RHOS = ("0.015", "0.018", "0.022", "0.027", "0.033", "0.040")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dps", type=int, default=70)
    p.add_argument("--ode-degree", type=int, default=30)
    p.add_argument("--ode-tol-digits", type=int, default=54)
    p.add_argument("--workers", type=int, default=max(1, min(2, os.cpu_count() or 1)))
    p.add_argument(
        "--constants",
        type=Path,
        default=None,
        help=(
            "Path to recovered_constants_high_order.csv. "
            "If omitted, the script searches the script folder and current folder."
        ),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("weight6_angle_holdout_outputs"),
    )
    return p.parse_args()


def mp_ls_s6(rhos, e6s):
    coeff, rms = w6.mp_least_squares(rhos, e6s, (0,2,4,6,8))
    return coeff[0], rms


def load_fit21_constants(path: Path):
    df = pd.read_csv(path, dtype=str)
    vals = {}
    for _, row in df.iterrows():
        vals[row["constant"]] = mp.mpf(row["fit21_high_order"])
    return vals


def model_prediction(constants, geometry, x):
    R, a, th, beta = map(mp.mpf, geometry)
    inv = w6.torus_weight6_invariants(R, a, th, beta)

    q = R + a*mp.cos(th)
    K = mp.cos(th)/(a*q)

    z = mp.sin(x/2)
    rr = mp.cos(x/2)
    t = rr**2

    Pk3 = (
        mp.mpf(31)/15120
        - mp.mpf(11)/280*t
        + mp.mpf(5)/112*t**2
    )

    names = [
        "K DeltaK",
        "K K_ww",
        "|grad K|^2",
        "K_w^2",
        "S_wwww",
        "(tr S)_ww",
        "tr^2 S",
    ]

    poly = mp.mpf("0")
    for I, name in zip(inv, names):
        A = constants[f"A[{name}]"]
        B = constants[f"B[{name}]"]
        C = constants[f"C[{name}]"]
        poly += I*(A + B*t + C*t**2)

    s6_pred = z*(1-t)*(K**3*Pk3 + poly)
    return s6_pred


def main():
    args = parse_args()
    mp.mp.dps = args.dps

    here = Path(__file__).resolve().parent

    if args.constants is not None:
        constants_path = args.constants.expanduser().resolve()
    else:
        candidates = [
            here / "recovered_constants_high_order.csv",
            here.parent / "data" / "recovered_constants_high_order.csv",
            Path.cwd() / "recovered_constants_high_order.csv",
        ]

        # Also tolerate browser/download renaming such as "... (1).csv".
        candidates.extend(sorted(here.glob("recovered_constants_high_order*.csv")))
        candidates.extend(sorted(Path.cwd().glob("recovered_constants_high_order*.csv")))

        constants_path = next(
            (p.resolve() for p in candidates if p.exists()),
            None,
        )

    if constants_path is None or not constants_path.exists():
        raise FileNotFoundError(
            "Could not find recovered_constants_high_order.csv. "
            "Either place it beside this script or run with "
            "--constants C:\\Users\\Wayne\\August\\<filename>.csv"
        )

    print(f"Using constants file: {constants_path}", flush=True)
    constants = load_fit21_constants(constants_path)

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    payloads = []
    for gid in TEST_GEOMETRY_IDS:
        geometry = w6.DESIGN_GEOMETRIES[gid-1]
        for rho in DEFAULT_RHOS:
            payloads.append(
                (
                    gid,
                    geometry,
                    rho,
                    DEFAULT_ANGLES,
                    args.dps,
                    args.ode_degree,
                    args.ode_tol_digits,
                )
            )

    print("WEIGHT-SIX OUT-OF-SAMPLE ANGLE VALIDATION", flush=True)
    print("="*72, flush=True)
    print("Training angles were: 0.4, 0.8, 1.2", flush=True)
    print("Validation angles are: 0.55, 0.95, 1.35", flush=True)
    print("Geometry 1 = angle-only holdout", flush=True)
    print("Geometries 8--10 = geometry + angle holdout", flush=True)
    print()

    rows = []
    with ProcessPoolExecutor(
        max_workers=args.workers,
        max_tasks_per_child=1,
    ) as pool:
        futs = {
            pool.submit(w6.evaluate_geometry_rho, p): (p[0], p[2])
            for p in payloads
        }
        for fut in as_completed(futs):
            gid, rho = futs[fut]
            _, _, outrows = fut.result()
            rows.extend(outrows)
            print(f"completed geometry={gid} rho={rho}", flush=True)

    diagnostics = [r for r in rows if r["direction"] == "EVEN_DIAGNOSTIC"]

    results = []
    for gid in TEST_GEOMETRY_IDS:
        geometry = w6.DESIGN_GEOMETRIES[gid-1]
        for xt in DEFAULT_ANGLES:
            x = mp.mpf(xt)
            sub = [
                r for r in diagnostics
                if int(r["geometry_id"]) == gid
                and abs(mp.mpf(r["x"]) - x) < mp.mpf("1e-50")
            ]
            sub.sort(key=lambda r: mp.mpf(r["rho"]))

            rhos = [mp.mpf(r["rho"]) for r in sub]
            e6s = [mp.mpf(r["E6"]) for r in sub]

            s6_obs, fit_rms = mp_ls_s6(rhos, e6s)
            s6_pred = model_prediction(constants, geometry, x)
            residual = s6_obs - s6_pred
            rel = abs(residual) / max(abs(s6_obs), mp.mpf("1e-80"))

            results.append({
                "geometry_id": gid,
                "role": "ANGLE_ONLY" if gid == 1 else "GEOMETRY_AND_ANGLE",
                "x": mp.nstr(x, 20),
                "s6_observed": mp.nstr(s6_obs, 60),
                "s6_predicted": mp.nstr(s6_pred, 60),
                "residual": mp.nstr(residual, 40),
                "relative_residual": mp.nstr(rel, 20),
                "fit_rms": mp.nstr(fit_rms, 30),
            })

    residuals = [mp.mpf(r["residual"]) for r in results]
    rms = mp.sqrt(mp.fsum(v*v for v in residuals) / len(residuals))
    maxabs = max(abs(v) for v in residuals)

    angle_only = [mp.mpf(r["residual"]) for r in results if r["role"] == "ANGLE_ONLY"]
    double_hold = [mp.mpf(r["residual"]) for r in results if r["role"] == "GEOMETRY_AND_ANGLE"]

    rms_angle = mp.sqrt(mp.fsum(v*v for v in angle_only)/len(angle_only))
    rms_double = mp.sqrt(mp.fsum(v*v for v in double_hold)/len(double_hold))

    csv_path = out / "angle_holdout_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    lines = [
        "WEIGHT-SIX OUT-OF-SAMPLE ANGLE VALIDATION SUMMARY",
        "="*72,
        "training angles:   0.4, 0.8, 1.2",
        "validation angles: 0.55, 0.95, 1.35",
        "",
        f"all 12 RMS residual       = {mp.nstr(rms,30)}",
        f"all 12 max |residual|     = {mp.nstr(maxabs,30)}",
        f"geometry 1 angle-only RMS = {mp.nstr(rms_angle,30)}",
        f"geometries 8--10 RMS      = {mp.nstr(rms_double,30)}",
        "",
    ]

    for row in results:
        lines.extend([
            f"geometry={row['geometry_id']} role={row['role']} x={row['x']}",
            f"  observed  = {row['s6_observed']}",
            f"  predicted = {row['s6_predicted']}",
            f"  residual  = {row['residual']}",
            f"  relative  = {row['relative_residual']}",
            f"  fit RMS   = {row['fit_rms']}",
        ])

    txt_path = out / "angle_holdout_summary.txt"
    txt_path.write_text("\n".join(lines), encoding="utf-8")

    print()
    print("\n".join(lines), flush=True)
    print()
    print("Wrote:", csv_path, flush=True)
    print("Wrote:", txt_path, flush=True)


if __name__ == "__main__":
    main()
