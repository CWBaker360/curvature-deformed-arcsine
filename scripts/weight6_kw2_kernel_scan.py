#!/usr/bin/env python3
"""
Direct high-precision scan of the unresolved weight-six K_w^2 angular kernel.

Background
----------
The 12 genuinely new-angle holdouts show that the quadratic-in-r^2
weight-six angular ansatz fails, but the residual is proportional to K_w^2
alone across all tested torus geometries.

This script therefore fixes:
  * the exact spherical K^3 row, and
  * the 18 other recovered weight-six constants at their simple rational values,

then directly extracts the remaining angular function P_Kw2(r) from

    s6 = z(1-r^2) [
        K^3 P_K3(r^2)
        + sum_{I != K_w^2} I P_I(r^2)
        + K_w^2 P_Kw2(r)
    ].

Two independent torus geometries are scanned at NEW side angles.
The resulting P_Kw2 values should agree between geometries.

The script then tests polynomial models in r of increasing degree and reports
the residual for each degree.  This determines the minimal finite kernel
supported by the data before any rational reconstruction is frozen.

Requires in the same folder:
    weight6_torus_full_recovery.py

Outputs:
    weight6_kw2_scan_outputs/
        kw2_kernel_values.csv
        kw2_kernel_fit_summary.txt
"""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path

import mpmath as mp

import weight6_torus_full_recovery as w6


# Two geometries not simultaneously degenerate in the gradient invariants.
GEOMETRY_IDS = (8, 9)

# All are new angles: none is 0.4, 0.8, or 1.2.
ANGLE_TEXTS = (
    "0.50",
    "0.60",
    "0.70",
    "0.90",
    "1.00",
    "1.10",
    "1.30",
    "1.40",
)

RHO_TEXTS = ("0.015", "0.018", "0.022", "0.027", "0.033", "0.040")


# Candidate invariant order returned by torus_weight6_invariants:
INV_NAMES = (
    "K DeltaK",
    "K K_ww",
    "|grad K|^2",
    "K_w^2",
    "S_wwww",
    "(tr S)_ww",
    "tr^2 S",
)


# Exact rational rows recovered at ~1e-20 from the high-order reprocessing.
# K_w^2 is deliberately omitted.
EXACT_ROWS = {
    "K DeltaK": (
        mp.mpf(0),
        -mp.mpf(1) / 112,
        mp.mpf(127) / 15120,
    ),
    "K K_ww": (
        mp.mpf(29) / 5040,
        mp.mpf(53) / 5040,
        -mp.mpf(53) / 1260,
    ),
    "|grad K|^2": (
        mp.mpf(0),
        -mp.mpf(17) / 10080,
        mp.mpf(71) / 30240,
    ),
    "S_wwww": (
        mp.mpf(1) / 1008,
        mp.mpf(1) / 1260,
        mp.mpf(1) / 315,
    ),
    "(tr S)_ww": (
        mp.mpf(0),
        mp.mpf(1) / 2520,
        -mp.mpf(1) / 420,
    ),
    "tr^2 S": (
        mp.mpf(0),
        -mp.mpf(1) / 5040,
        mp.mpf(1) / 5040,
    ),
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dps", type=int, default=70)
    p.add_argument("--ode-degree", type=int, default=30)
    p.add_argument("--ode-tol-digits", type=int, default=54)
    p.add_argument(
        "--workers",
        type=int,
        default=max(1, min(2, os.cpu_count() or 1)),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("weight6_kw2_scan_outputs"),
    )
    return p.parse_args()


def fit_s6(rhos, e6s):
    coeff, rms = w6.mp_least_squares(
        rhos,
        e6s,
        (0, 2, 4, 6, 8),
    )
    return coeff[0], rms


def poly(A, B, C, t):
    return A + B*t + C*t*t


def extract_kw2_kernel(geometry, x, s6):
    R, a, th, beta = map(mp.mpf, geometry)
    inv = w6.torus_weight6_invariants(R, a, th, beta)

    q = R + a*mp.cos(th)
    K = mp.cos(th) / (a*q)

    z = mp.sin(x/2)
    r = mp.cos(x/2)
    t = r*r

    Pk3 = (
        mp.mpf(31)/15120
        - mp.mpf(11)/280*t
        + mp.mpf(5)/112*t*t
    )

    known = K**3 * Pk3

    for name, I in zip(INV_NAMES, inv):
        if name == "K_w^2":
            continue
        A, B, C = EXACT_ROWS[name]
        known += I * poly(A, B, C, t)

    kw2 = inv[3]
    if abs(kw2) < mp.mpf("1e-30"):
        raise ZeroDivisionError("K_w^2 is too small at this geometry.")

    Pkw2 = (
        s6 / (z*(1-t))
        - known
    ) / kw2

    return r, t, Pkw2, kw2


def polynomial_fit(xs, ys, degree):
    A = mp.matrix(len(xs), degree + 1)
    b = mp.matrix(len(xs), 1)

    for i, (x, y) in enumerate(zip(xs, ys)):
        b[i] = y
        for j in range(degree + 1):
            A[i, j] = x**j

    coeff, _ = mp.qr_solve(A, b)

    residuals = []
    for x, y in zip(xs, ys):
        pred = mp.fsum(coeff[j]*x**j for j in range(degree + 1))
        residuals.append(pred-y)

    rms = mp.sqrt(mp.fsum(v*v for v in residuals)/len(residuals))
    maxabs = max(abs(v) for v in residuals)

    return [coeff[j] for j in range(degree+1)], rms, maxabs


def rational_hint(x, max_den=10_000_000):
    f = Fraction(mp.nstr(x, 50)).limit_denominator(max_den)
    q = mp.mpf(f.numerator)/f.denominator
    return f"{f.numerator}/{f.denominator}", abs(x-q)


def main():
    args = parse_args()
    mp.mp.dps = args.dps

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    payloads = []
    for gid in GEOMETRY_IDS:
        geom = w6.DESIGN_GEOMETRIES[gid-1]
        for rho in RHO_TEXTS:
            payloads.append(
                (
                    gid,
                    geom,
                    rho,
                    ANGLE_TEXTS,
                    args.dps,
                    args.ode_degree,
                    args.ode_tol_digits,
                )
            )

    print("WEIGHT-SIX K_w^2 ANGULAR-KERNEL SCAN", flush=True)
    print("="*72, flush=True)
    print("geometries:", GEOMETRY_IDS, flush=True)
    print("new angles:", ", ".join(ANGLE_TEXTS), flush=True)
    print("all other weight-six rows fixed to exact rational candidates", flush=True)
    print()

    raw_rows = []

    with ProcessPoolExecutor(
        max_workers=args.workers,
        max_tasks_per_child=1,
    ) as pool:
        futures = {
            pool.submit(w6.evaluate_geometry_rho, p): (p[0], p[2])
            for p in payloads
        }
        for fut in as_completed(futures):
            gid, rho = futures[fut]
            _, _, rows = fut.result()
            raw_rows.extend(rows)
            print(f"completed geometry={gid} rho={rho}", flush=True)

    diagnostics = [
        row for row in raw_rows
        if row["direction"] == "EVEN_DIAGNOSTIC"
    ]

    kernel_rows = []

    for gid in GEOMETRY_IDS:
        geom = w6.DESIGN_GEOMETRIES[gid-1]

        for xt in ANGLE_TEXTS:
            x = mp.mpf(xt)

            subset = [
                row for row in diagnostics
                if int(row["geometry_id"]) == gid
                and abs(mp.mpf(row["x"]) - x) < mp.mpf("1e-50")
            ]
            subset.sort(key=lambda row: mp.mpf(row["rho"]))

            rhos = [mp.mpf(row["rho"]) for row in subset]
            e6s = [mp.mpf(row["E6"]) for row in subset]

            s6, fit_rms = fit_s6(rhos, e6s)
            r, t, Pkw2, kw2 = extract_kw2_kernel(geom, x, s6)

            kernel_rows.append(
                {
                    "geometry_id": gid,
                    "x": mp.nstr(x, 20),
                    "r": mp.nstr(r, 60),
                    "t": mp.nstr(t, 60),
                    "s6": mp.nstr(s6, 60),
                    "K_w^2": mp.nstr(kw2, 60),
                    "P_Kw2": mp.nstr(Pkw2, 60),
                    "s6_fit_rms": mp.nstr(fit_rms, 30),
                }
            )

    # Compare the two geometries angle by angle and use their mean.
    averaged = []

    for xt in ANGLE_TEXTS:
        vals = [
            mp.mpf(row["P_Kw2"])
            for row in kernel_rows
            if abs(mp.mpf(row["x"]) - mp.mpf(xt)) < mp.mpf("1e-40")
        ]
        if len(vals) != len(GEOMETRY_IDS):
            raise RuntimeError("Missing geometry result.")

        mean = mp.fsum(vals)/len(vals)
        spread = max(abs(v-mean) for v in vals)

        x = mp.mpf(xt)
        r = mp.cos(x/2)

        averaged.append((r, mean, spread))

    # Fit finite polynomials in r.  Degree 2...7 are informative.
    summaries = []

    for degree in range(2, 8):
        coeff, rms, maxabs = polynomial_fit(
            [row[0] for row in averaged],
            [row[1] for row in averaged],
            degree,
        )

        summaries.append(
            (degree, coeff, rms, maxabs)
        )

    csv_path = out / "kw2_kernel_values.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(kernel_rows[0].keys()))
        writer.writeheader()
        writer.writerows(kernel_rows)

    lines = [
        "WEIGHT-SIX K_w^2 ANGULAR-KERNEL SCAN SUMMARY",
        "="*72,
        "",
        "Geometry agreement at each new angle:",
    ]

    for (r, mean, spread), xt in zip(averaged, ANGLE_TEXTS):
        lines.extend(
            [
                f"x={xt}",
                f"  r             = {mp.nstr(r,30)}",
                f"  P_Kw2 mean    = {mp.nstr(mean,50)}",
                f"  geometry spread = {mp.nstr(spread,20)}",
            ]
        )

    lines.extend(
        [
            "",
            "Polynomial-in-r model comparison:",
        ]
    )

    for degree, coeff, rms, maxabs in summaries:
        lines.extend(
            [
                f"degree {degree}:",
                f"  RMS     = {mp.nstr(rms,25)}",
                f"  max abs = {mp.nstr(maxabs,25)}",
            ]
        )

        if degree >= 4:
            lines.append("  coefficients:")
            for j, c in enumerate(coeff):
                rat, err = rational_hint(c)
                lines.append(
                    f"    r^{j}: {mp.nstr(c,35)}"
                    f"   ~ {rat}  (err {mp.nstr(err,8)})"
                )

    best = summaries[-1]
    lines.extend(
        [
            "",
            "Interpretation:",
            (
                "The minimal degree whose residual is comparable to the "
                "geometry-to-geometry extraction spread is the numerical "
                "candidate for the true finite K_w^2 angular kernel."
            ),
            (
                "Do not freeze rational coefficients until that degree is "
                "stable under an additional angle or radius-window check."
            ),
        ]
    )

    txt_path = out / "kw2_kernel_fit_summary.txt"
    txt_path.write_text("\n".join(lines), encoding="utf-8")

    print()
    print("\n".join(lines), flush=True)
    print()
    print("Wrote:", csv_path, flush=True)
    print("Wrote:", txt_path, flush=True)


if __name__ == "__main__":
    main()
