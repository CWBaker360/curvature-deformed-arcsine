#!/usr/bin/env python3
"""
High-precision pilot for extracting the weight-six normalized side-law
coefficient s6 on a general ring torus.

Why this pilot exists
---------------------
Before running the full 30-observation / 21-constant recovery, we first
verify that the rho^6 coefficient can be extracted stably from intrinsic
geodesic distances at one generic torus geometry.

Key numerical improvements
--------------------------
1. All geometry uses mpmath arbitrary precision.
2. The normalization c'(0) is NOT estimated from tiny angular samples.
   It is computed directly from the scalar Jacobi equation along the
   central radial geodesic:
       j'' + K j = 0,  j(0)=0, j'(0)=1.
   Thus c'(0)=j(rho).
3. Opposite central directions w and -w are averaged.  This cancels all
   odd geometric weights, so after subtracting the known rho^2 and rho^4
   terms,
       E6(rho) = s6 + O(rho^2).
4. Three moderate side angles are tested independently.

Metric
------
    ds^2 = a^2 du^2 + (R + a cos u)^2 dv^2.

Default pilot geometry
----------------------
    R=3, a=1, u0=1.1, beta=0,
which is the first geometry in the full-rank weight-six design.

Dependencies
------------
    mpmath >= 1.3
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Sequence

import mpmath as mp


def mpstr(x: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(x, digits, strip_zeros=False)


def parse_list(text: str) -> list[str]:
    vals = [part.strip() for part in text.split(",") if part.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("Expected comma-separated values.")
    return vals


@dataclass
class SideResult:
    rho: str
    x: str
    direction: str
    side_length: str
    jacobi_cprime: str
    normalized_S: str
    shooting_residual: str
    shooting_iterations: int


class HPTorus:
    def __init__(
        self,
        R: mp.mpf,
        a: mp.mpf,
        u0: mp.mpf,
        beta: mp.mpf,
        dps: int,
        ode_degree: int,
        ode_tol_digits: int,
    ) -> None:
        mp.mp.dps = dps
        self.R = mp.mpf(R)
        self.a = mp.mpf(a)
        self.u0 = mp.mpf(u0)
        self.beta = mp.mpf(beta)
        self.ode_degree = ode_degree
        self.ode_tol = mp.power(10, -ode_tol_digits)
        self.shoot_tol = mp.power(10, -(ode_tol_digits - 8))

    def q(self, u: mp.mpf) -> mp.mpf:
        return self.R + self.a * mp.cos(u)

    def K(self, u: mp.mpf) -> mp.mpf:
        return mp.cos(u) / (self.a * self.q(u))

    def rhs4(self, _s: mp.mpf, y: Sequence[mp.mpf]):
        u, v, up, vp = y
        q = self.q(u)
        return [
            up,
            vp,
            -q * mp.sin(u) / self.a * vp * vp,
            2 * self.a * mp.sin(u) / q * up * vp,
        ]

    def rhs6(self, _s: mp.mpf, y: Sequence[mp.mpf]):
        u, v, up, vp, j, jp = y
        q = self.q(u)
        return [
            up,
            vp,
            -q * mp.sin(u) / self.a * vp * vp,
            2 * self.a * mp.sin(u) / q * up * vp,
            jp,
            -self.K(u) * j,
        ]

    def integrate4(
        self,
        u: mp.mpf,
        v: mp.mpf,
        beta: mp.mpf,
        length: mp.mpf,
    ):
        q0 = self.q(u)
        y0 = [
            u,
            v,
            mp.cos(beta) / self.a,
            mp.sin(beta) / q0,
        ]
        sol = mp.odefun(
            self.rhs4,
            mp.mpf("0"),
            y0,
            tol=self.ode_tol,
            degree=self.ode_degree,
        )
        return sol(length)

    def endpoint(self, beta: mp.mpf, rho: mp.mpf):
        return self.integrate4(self.u0, mp.mpf("0"), beta, rho)

    def central_jacobi(self, beta: mp.mpf, rho: mp.mpf) -> mp.mpf:
        q0 = self.q(self.u0)
        y0 = [
            self.u0,
            mp.mpf("0"),
            mp.cos(beta) / self.a,
            mp.sin(beta) / q0,
            mp.mpf("0"),
            mp.mpf("1"),
        ]
        sol = mp.odefun(
            self.rhs6,
            mp.mpf("0"),
            y0,
            tol=self.ode_tol,
            degree=self.ode_degree,
        )
        y = sol(rho)
        return y[4]

    def intrinsic_distance(self, p1, p2, initial_length: mp.mpf):
        u1, v1 = p1[0], p1[1]
        u2, v2 = p2[0], p2[1]

        du = u2 - u1
        dv = v2 - v1
        q1 = self.q(u1)

        beta0 = mp.atan2(q1 * dv, self.a * du)
        L0 = max(
            initial_length,
            mp.sqrt((self.a * du) ** 2 + (q1 * dv) ** 2),
        )

        iterations = 0

        def residual(beta, logL):
            nonlocal iterations
            iterations += 1
            L = mp.exp(logL)
            end = self.integrate4(u1, v1, beta, L)
            r1 = self.a * (end[0] - u2)
            r2 = self.q(u2) * (end[1] - v2)
            return r1, r2

        beta, logL = mp.findroot(
            residual,
            (beta0, mp.log(L0)),
            solver="mdnewton",
            tol=self.shoot_tol,
            maxsteps=24,
            verify=False,
        )
        L = mp.exp(logL)

        end = self.integrate4(u1, v1, beta, L)
        r1 = self.a * (end[0] - u2)
        r2 = self.q(u2) * (end[1] - v2)
        final_res = mp.sqrt(r1 * r1 + r2 * r2)

        if final_res > self.shoot_tol * 100:
            raise RuntimeError(
                "Shooting residual too large: "
                + mp.nstr(final_res, 12)
            )

        return L, final_res, iterations

    def normalized_side(
        self,
        rho: mp.mpf,
        x: mp.mpf,
        beta: mp.mpf,
    ):
        p1 = self.endpoint(beta - x / 2, rho)
        p2 = self.endpoint(beta + x / 2, rho)

        euclid_guess = 2 * rho * mp.sin(x / 2)
        side, residual, iterations = self.intrinsic_distance(
            p1,
            p2,
            initial_length=euclid_guess,
        )

        cprime = self.central_jacobi(beta, rho)
        S = side / (2 * cprime)
        return side, cprime, S, residual, iterations


# ---------------------------------------------------------------------
# Base-point invariants needed for the known even orders rho^2,rho^4.
# These are evaluated directly from the torus metric.
# ---------------------------------------------------------------------

def torus_even_invariants(
    R: mp.mpf,
    a: mp.mpf,
    u: mp.mpf,
    beta: mp.mpf,
):
    q = R + a * mp.cos(u)
    K = mp.cos(u) / (a * q)

    # K_u and K_uu.
    Ku = -R * mp.sin(u) / (a * q**2)
    Kuu = (
        -R * mp.cos(u) / (a * q**2)
        - 2 * R * mp.sin(u)**2 / q**3
    )

    # Hessian in the orthonormal frame (e_u,e_v).
    # e_u = a^{-1} partial_u, e_v = q^{-1} partial_v.
    Huu = Kuu / a**2

    # Hess_{vv} K = - Gamma^u_vv K_u / q^2,
    # Gamma^u_vv = q sin(u)/a.
    Hvv = -(mp.sin(u) / (a * q)) * Ku

    cb = mp.cos(beta)
    sb = mp.sin(beta)
    Kww = cb**2 * Huu + sb**2 * Hvv

    DeltaK = Huu + Hvv

    return K, DeltaK, Kww


def normalized_known_even_coeffs(
    K: mp.mpf,
    DeltaK: mp.mpf,
    Kww: mp.mpf,
    x: mp.mpf,
):
    """
    Return s2(x), s4(x) in
        S_rho(x/2) = sin(x/2) + s2 rho^2 + s4 rho^4 + ...
    from the frozen normalized side law.
    """
    u = x / 2
    z = mp.sin(u)
    r = mp.cos(u)

    f2 = -K * r**2 / 6
    f20 = -K / 6

    B4 = r**2 / 120 * (
        (1 - 9*z**2) * K**2
        - 3*r**2 * Kww
        - z**2 * (DeltaK - Kww)
    )
    B40 = (K**2 - 3*Kww) / 120

    # Coefficients of F/F0 through rho^4:
    g2 = f2 - f20
    g4 = B4 - B40 + f20**2 - f2*f20

    return z*g2, z*g4


def mp_least_squares(
    xs: Sequence[mp.mpf],
    ys: Sequence[mp.mpf],
    powers: Sequence[int],
):
    A = mp.matrix(len(xs), len(powers))
    b = mp.matrix(len(xs), 1)

    for i, (x, y) in enumerate(zip(xs, ys)):
        b[i] = y
        for j, power in enumerate(powers):
            A[i, j] = x**power

    coeff = mp.lu_solve(A.T*A, A.T*b)

    residuals = []
    for i, x in enumerate(xs):
        model = mp.fsum(
            coeff[j] * x**powers[j]
            for j in range(len(powers))
        )
        residuals.append(ys[i] - model)

    rms = mp.sqrt(
        mp.fsum(v*v for v in residuals) / len(residuals)
    )
    return [coeff[j] for j in range(len(powers))], rms


def evaluate_rho(payload):
    (
        rho_text,
        R_text,
        a_text,
        u0_text,
        beta_text,
        x_texts,
        dps,
        ode_degree,
        ode_tol_digits,
    ) = payload

    mp.mp.dps = dps

    R = mp.mpf(R_text)
    a = mp.mpf(a_text)
    u0 = mp.mpf(u0_text)
    beta = mp.mpf(beta_text)
    rho = mp.mpf(rho_text)
    xs = [mp.mpf(x) for x in x_texts]

    geo = HPTorus(
        R, a, u0, beta,
        dps=dps,
        ode_degree=ode_degree,
        ode_tol_digits=ode_tol_digits,
    )

    K, DeltaK, Kww = torus_even_invariants(
        R, a, u0, beta
    )

    rows = []

    for x in xs:
        values = {}

        for label, direction_beta in (
            ("plus", beta),
            ("minus", beta + mp.pi),
        ):
            side, cprime, S, res, iters = geo.normalized_side(
                rho,
                x,
                direction_beta,
            )
            values[label] = S

            rows.append(
                asdict(
                    SideResult(
                        rho=mpstr(rho),
                        x=mpstr(x),
                        direction=label,
                        side_length=mpstr(side),
                        jacobi_cprime=mpstr(cprime),
                        normalized_S=mpstr(S),
                        shooting_residual=mpstr(res),
                        shooting_iterations=iters,
                    )
                )
            )

        Seven = (values["plus"] + values["minus"]) / 2

        z = mp.sin(x / 2)
        s2, s4 = normalized_known_even_coeffs(
            K, DeltaK, Kww, x
        )

        E6 = (
            Seven
            - z
            - s2*rho**2
            - s4*rho**4
        ) / rho**6

        rows.append(
            {
                "rho": mpstr(rho),
                "x": mpstr(x),
                "direction": "EVEN_DIAGNOSTIC",
                "side_length": "",
                "jacobi_cprime": "",
                "normalized_S": mpstr(Seven),
                "shooting_residual": "",
                "shooting_iterations": "",
                "K": mpstr(K),
                "DeltaK": mpstr(DeltaK),
                "Kww": mpstr(Kww),
                "s2": mpstr(s2),
                "s4": mpstr(s4),
                "E6": mpstr(E6),
            }
        )

    return rho_text, rows


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument("--R", default="3.0")
    p.add_argument("--minor-radius", dest="a", default="1.0")
    p.add_argument("--u0", default="1.1")
    p.add_argument("--beta", default="0.0")

    p.add_argument(
        "--rhos",
        type=parse_list,
        default=parse_list(
            "0.015,0.018,0.022,0.027,0.033,0.040"
        ),
    )
    p.add_argument(
        "--xs",
        type=parse_list,
        default=parse_list("0.4,0.8,1.2"),
    )

    p.add_argument("--dps", type=int, default=60)
    p.add_argument("--ode-degree", type=int, default=28)
    p.add_argument("--ode-tol-digits", type=int, default=46)
    p.add_argument(
        "--workers",
        type=int,
        default=max(1, min(2, os.cpu_count() or 1)),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("weight6_pilot_outputs"),
    )
    return p.parse_args()


def main():
    args = parse_args()

    if args.dps < 50:
        raise ValueError("Use at least 50 decimal digits.")

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    payloads = [
        (
            rho,
            args.R,
            args.a,
            args.u0,
            args.beta,
            args.xs,
            args.dps,
            args.ode_degree,
            args.ode_tol_digits,
        )
        for rho in args.rhos
    ]

    print("WEIGHT-SIX TORUS s6 PILOT", flush=True)
    print("="*64, flush=True)
    print(
        f"R={args.R}, a={args.a}, u0={args.u0}, beta={args.beta}",
        flush=True,
    )
    print(
        f"dps={args.dps}, radii={len(args.rhos)}, "
        f"side angles={len(args.xs)}, workers={args.workers}",
        flush=True,
    )
    print(
        "Opposite-direction averaging enabled: "
        "E6(rho)=s6+O(rho^2).",
        flush=True,
    )
    print()

    all_rows = []

    with ProcessPoolExecutor(
        max_workers=args.workers,
        max_tasks_per_child=1,
    ) as pool:
        futures = {
            pool.submit(evaluate_rho, payload): payload[0]
            for payload in payloads
        }

        for future in as_completed(futures):
            rho = futures[future]
            _, rows = future.result()
            all_rows.extend(rows)
            print(f"completed rho={rho}", flush=True)

    # Sort for stable output.
    def sort_key(row):
        try:
            rr = float(row["rho"])
        except Exception:
            rr = 0.0
        try:
            xx = float(row["x"])
        except Exception:
            xx = 0.0
        order = {
            "plus": 0,
            "minus": 1,
            "EVEN_DIAGNOSTIC": 2,
        }.get(row["direction"], 9)
        return rr, xx, order

    all_rows.sort(key=sort_key)

    # Extract E6 values by x and extrapolate with even powers.
    diagnostics = [
        row for row in all_rows
        if row["direction"] == "EVEN_DIAGNOSTIC"
    ]

    mp.mp.dps = args.dps
    summary = {
        "geometry": {
            "R": args.R,
            "minor_radius": args.a,
            "u0": args.u0,
            "beta": args.beta,
        },
        "dps": args.dps,
        "ode_degree": args.ode_degree,
        "ode_tol_digits": args.ode_tol_digits,
        "fits": {},
    }

    for x_text in args.xs:
        subset = [
            row for row in diagnostics
            if abs(mp.mpf(row["x"]) - mp.mpf(x_text))
            < mp.mpf("1e-40")
        ]
        subset.sort(key=lambda row: mp.mpf(row["rho"]))

        rhos = [mp.mpf(row["rho"]) for row in subset]
        E6s = [mp.mpf(row["E6"]) for row in subset]

        fit3, rms3 = mp_least_squares(
            rhos, E6s, (0, 2, 4)
        )
        fit2, rms2 = mp_least_squares(
            rhos, E6s, (0, 2)
        )

        # Nested smallest-radius fit as a stability diagnostic.
        nsmall = max(4, len(rhos) - 1)
        fit_small, rms_small = mp_least_squares(
            rhos[:nsmall],
            E6s[:nsmall],
            (0, 2, 4),
        )

        q0 = fit3[0]
        q0_small = fit_small[0]
        spread = abs(q0 - q0_small)

        summary["fits"][x_text] = {
            "s6_extrapolated_0_2_4": mpstr(q0),
            "s6_extrapolated_0_2": mpstr(fit2[0]),
            "s6_nested_small_radius": mpstr(q0_small),
            "nested_spread": mpstr(spread),
            "fit_rms_0_2_4": mpstr(rms3),
            "fit_rms_0_2": mpstr(rms2),
            "fit_rms_nested": mpstr(rms_small),
            "E6_values": [
                {
                    "rho": mpstr(rr),
                    "E6": mpstr(ee),
                }
                for rr, ee in zip(rhos, E6s)
            ],
        }

    # Write CSV with union of keys.
    csv_path = out / "pilot_raw_and_diagnostics.csv"
    keys = []
    for row in all_rows:
        for key in row:
            if key not in keys:
                keys.append(key)

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=keys,
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    json_path = out / "pilot_summary.json"
    json_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    txt_lines = [
        "WEIGHT-SIX TORUS s6 PILOT SUMMARY",
        "="*64,
        (
            f"R={args.R}, a={args.a}, "
            f"u0={args.u0}, beta={args.beta}"
        ),
        (
            f"dps={args.dps}, ode_degree={args.ode_degree}, "
            f"ode_tol_digits={args.ode_tol_digits}"
        ),
        "",
        "Opposite-direction average:",
        "  E6(rho)=s6 + O(rho^2)",
        "",
    ]

    for x_text in args.xs:
        fit = summary["fits"][x_text]
        txt_lines.extend(
            [
                f"x={x_text}",
                (
                    "  s6 [1,rho^2,rho^4] = "
                    f"{fit['s6_extrapolated_0_2_4']}"
                ),
                (
                    "  s6 [1,rho^2]       = "
                    f"{fit['s6_extrapolated_0_2']}"
                ),
                (
                    "  nested spread      = "
                    f"{fit['nested_spread']}"
                ),
                (
                    "  fit RMS            = "
                    f"{fit['fit_rms_0_2_4']}"
                ),
                "  E6 values:",
            ]
        )
        for item in fit["E6_values"]:
            txt_lines.append(
                f"    rho={item['rho'][:12]:<12} "
                f"E6={item['E6']}"
            )
        txt_lines.append("")

    txt_path = out / "pilot_summary.txt"
    txt_path.write_text(
        "\n".join(txt_lines),
        encoding="utf-8",
    )

    print()
    print("\n".join(txt_lines), flush=True)
    print("Wrote:", csv_path, flush=True)
    print("Wrote:", json_path, flush=True)
    print("Wrote:", txt_path, flush=True)


if __name__ == "__main__":
    main()
