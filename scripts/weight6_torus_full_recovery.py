#!/usr/bin/env python3
"""
High-precision recovery of the 21 unknown weight-six side-law constants.

Model
-----
For the normalized half-side law

    S_rho(u) = sin(u) + rho^2 s2 + rho^3 s3 + ...,

opposite-direction averaging removes odd geometric weights.  After
subtracting the known rho^2 and rho^4 terms,

    E6(rho) = s6 + O(rho^2).

At geometric weight six, use the candidate angular form

    s6 = z(1-r^2) sum_I I [A_I + B_I r^2 + C_I r^4],

where z=sin(x/2), r=cos(x/2), and the K^3 row is already fixed by the
exact sphere:

    P_K3(t) = 31/15120 - 11/280 t + 5/112 t^2,
    t=r^2.

The seven remaining candidate invariants are

    K DeltaK,
    K K_ww,
    |grad K|^2,
    K_w^2,
    S_wwww,
    (tr S)_ww,
    tr^2 S,

with S = Sym nabla^4 K.

Experimental design
-------------------
Ten torus geometries x three side angles = 30 observations.

Geometries 1--7 are used for the 21 x 21 recovery system.
Geometries 8--10 are NEVER used to determine the constants and supply
nine genuine holdout equations.

An independent all-30 least-squares solution is also computed and
compared with the 21-equation recovery.

Numerics
--------
- mpmath arbitrary precision
- Jacobi-equation normalization c'(0)=j(rho)
- opposite-direction averaging
- intrinsic two-point geodesic shooting
- even extrapolation [1, rho^2, rho^4]

Dependencies
------------
mpmath

Outputs
-------
weight6_recovery_outputs/
    raw_geodesic_data.csv
    s6_observations.csv
    recovered_constants.csv
    recovery_summary.json
    recovery_summary.txt
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import mpmath as mp


# ---------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------

def mpstr(x: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(x, digits, strip_zeros=False)


def parse_list(text: str) -> list[str]:
    vals = [part.strip() for part in text.split(",") if part.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("Expected comma-separated values.")
    return vals


def rational_candidate(x: mp.mpf, max_denominator: int = 100_000_000):
    s = mp.nstr(x, 50)
    frac = Fraction(s).limit_denominator(max_denominator)
    q = mp.mpf(frac.numerator) / frac.denominator
    err = abs(x - q)
    return frac, err


# ---------------------------------------------------------------------
# Intrinsic torus geodesic engine
# ---------------------------------------------------------------------

@dataclass
class RawSide:
    geometry_id: int
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
        return sol(rho)[4]

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

    def normalized_side_with_cprime(
        self,
        rho: mp.mpf,
        x: mp.mpf,
        beta: mp.mpf,
        cprime: mp.mpf,
    ):
        p1 = self.endpoint(beta - x / 2, rho)
        p2 = self.endpoint(beta + x / 2, rho)

        euclid_guess = 2 * rho * mp.sin(x / 2)
        side, residual, iterations = self.intrinsic_distance(
            p1,
            p2,
            initial_length=euclid_guess,
        )

        S = side / (2 * cprime)
        return side, S, residual, iterations


# ---------------------------------------------------------------------
# Torus invariant formulas
# ---------------------------------------------------------------------

def torus_even_invariants(R, a, u, beta):
    q = R + a * mp.cos(u)
    K = mp.cos(u) / (a * q)

    Ku = -R * mp.sin(u) / (a * q**2)
    Kuu = (
        -R * mp.cos(u) / (a * q**2)
        - 2 * R * mp.sin(u)**2 / q**3
    )

    Huu = Kuu / a**2
    Hvv = -(mp.sin(u) / (a * q)) * Ku

    cb = mp.cos(beta)
    sb = mp.sin(beta)

    Kww = cb**2 * Huu + sb**2 * Hvv
    DeltaK = Huu + Hvv

    return K, DeltaK, Kww


def torus_weight6_invariants(R, a, th, beta):
    """
    Closed formulas obtained from the metric and the fully symmetrized
    fourth covariant derivative S=Sym nabla^4 K.

    Return order:
      K DeltaK,
      K K_ww,
      |grad K|^2,
      K_w^2,
      S_wwww,
      (tr S)_ww,
      tr^2 S.
    """
    st = mp.sin(th)
    ct = mp.cos(th)
    sb = mp.sin(beta)
    cb = mp.cos(beta)
    q = R + a*ct

    KD = -R*(R*ct + a)*ct / (a**4 * q**4)

    KKww = (
        R
        * (
            a*sb**2*st**2
            + (-R*ct + a*ct**2 - 2*a)*cb**2
        )
        * ct
        / (a**4 * q**4)
    )

    gradK2 = R**2 * st**2 / (a**4 * q**4)
    Kw2 = R**2 * st**2 * cb**2 / (a**4 * q**4)

    Swwww = (
        R
        * (
            a**2*(-4*R*ct + 5*a*ct**2 - 9*a)
            * sb**4 * st**2
            + a*(
                11*R**2*ct**2
                - 7*R**2
                - 43*R*a*ct**3
                + 51*R*a*ct
                + 18*a**2*st**4
                - 50*a**2*ct**2
                + 54*a**2
            ) * sb**2 * cb**2
            + (
                R**3*ct
                - 11*R**2*a*ct**2
                + 8*R**2*a
                + 11*R*a**2*ct**3
                - 20*R*a**2*ct
                - a**3*st**4
                + 18*a**3*ct**2
                - 23*a**3
            ) * cb**4
        )
        / (a**5 * q**5)
    )

    trSww = (
        -R
        * (
            a*(
                -11*R**2*ct**2
                + 7*R**2
                + 19*R*a*ct**3
                - 27*R*a*ct
                + 12*a**2*st**4
                + 26*a**2*ct**2
                - 30*a**2
            ) * sb**2
            - (
                6*R**3*ct
                - 55*R**2*a*ct**2
                + 41*R**2*a
                + 23*R*a**2*ct**3
                - 69*R*a**2*ct
                + 12*a**3*st**4
                + 58*a**3*ct**2
                - 84*a**3
            ) * cb**2
        )
        / (6*a**5 * q**5)
    )

    tr2S = (
        R
        * (
            3*R**3*ct
            - 22*R**2*a*ct**2
            + 17*R**2*a
            + 2*R*a**2*ct**3
            - 21*R*a**2*ct
            + 16*a**3*ct**2
            - 27*a**3
        )
        / (3*a**5 * q**5)
    )

    return [KD, KKww, gradK2, Kw2, Swwww, trSww, tr2S]


def normalized_known_even_coeffs(K, DeltaK, Kww, x):
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

    g2 = f2 - f20
    g4 = B4 - B40 + f20**2 - f2*f20

    return z*g2, z*g4


# ---------------------------------------------------------------------
# Extrapolation and recovery
# ---------------------------------------------------------------------

def mp_least_squares(xs, ys, powers):
    A = mp.matrix(len(xs), len(powers))
    b = mp.matrix(len(xs), 1)

    for i, (x, y) in enumerate(zip(xs, ys)):
        b[i] = y
        for j, p in enumerate(powers):
            A[i, j] = x**p

    coeff, _ = mp.qr_solve(A, b)

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


def solve_column_scaled_square(M, b):
    nrow, ncol = M.rows, M.cols
    if nrow != ncol:
        raise ValueError("Expected square matrix.")

    norms = []
    Ms = mp.matrix(nrow, ncol)

    for j in range(ncol):
        norm = mp.sqrt(mp.fsum(M[i,j]**2 for i in range(nrow)))
        if norm == 0:
            raise ValueError("Zero design column.")
        norms.append(norm)
        for i in range(nrow):
            Ms[i,j] = M[i,j] / norm

    y = mp.lu_solve(Ms, b)
    x = mp.matrix(ncol, 1)
    for j in range(ncol):
        x[j] = y[j] / norms[j]
    return x


def solve_column_scaled_ls(M, b):
    nrow, ncol = M.rows, M.cols
    norms = []
    Ms = mp.matrix(nrow, ncol)

    for j in range(ncol):
        norm = mp.sqrt(mp.fsum(M[i,j]**2 for i in range(nrow)))
        if norm == 0:
            raise ValueError("Zero design column.")
        norms.append(norm)
        for i in range(nrow):
            Ms[i,j] = M[i,j] / norm

    y, _ = mp.qr_solve(Ms, b)
    x = mp.matrix(ncol, 1)
    for j in range(ncol):
        x[j] = y[j] / norms[j]
    return x


def matrix_predict(M, x):
    return M*x


def residual_stats(M, x, b):
    pred = matrix_predict(M, x)
    res = [pred[i] - b[i] for i in range(M.rows)]
    rms = mp.sqrt(mp.fsum(v*v for v in res) / len(res))
    maxabs = max(abs(v) for v in res)
    return res, rms, maxabs


# ---------------------------------------------------------------------
# Parallel intrinsic-data calculation
# ---------------------------------------------------------------------

DESIGN_GEOMETRIES = [
    ("3.0", "1.0", "1.1", "0.0"),
    ("3.0", "0.7", "2.0", "1.4"),
    ("5.0", "0.7", "2.5", "0.0"),
    ("3.0", "0.7", "0.7", "1.4"),
    ("3.0", "0.7", "2.0", "1.05"),
    ("3.0", "0.7", "0.3", "0.0"),
    ("3.0", "0.7", "2.5", "1.4"),
    ("4.0", "0.7", "0.3", "0.35"),
    ("3.0", "1.0", "2.0", "0.35"),
    ("4.0", "1.0", "1.5", "0.0"),
]


def evaluate_geometry_rho(payload):
    (
        geometry_id,
        geometry,
        rho_text,
        x_texts,
        dps,
        ode_degree,
        ode_tol_digits,
    ) = payload

    mp.mp.dps = dps

    R = mp.mpf(geometry[0])
    a = mp.mpf(geometry[1])
    u0 = mp.mpf(geometry[2])
    beta = mp.mpf(geometry[3])
    rho = mp.mpf(rho_text)
    xs = [mp.mpf(v) for v in x_texts]

    geo = HPTorus(
        R, a, u0, beta,
        dps=dps,
        ode_degree=ode_degree,
        ode_tol_digits=ode_tol_digits,
    )

    K, DeltaK, Kww = torus_even_invariants(R, a, u0, beta)

    # Jacobi normalizations depend only on the central direction and rho.
    cprime_plus = geo.central_jacobi(beta, rho)
    cprime_minus = geo.central_jacobi(beta + mp.pi, rho)

    rows = []

    for x in xs:
        side_p, Sp, res_p, it_p = geo.normalized_side_with_cprime(
            rho, x, beta, cprime_plus
        )
        side_m, Sm, res_m, it_m = geo.normalized_side_with_cprime(
            rho, x, beta + mp.pi, cprime_minus
        )

        rows.append(
            asdict(
                RawSide(
                    geometry_id=geometry_id,
                    rho=mpstr(rho),
                    x=mpstr(x),
                    direction="plus",
                    side_length=mpstr(side_p),
                    jacobi_cprime=mpstr(cprime_plus),
                    normalized_S=mpstr(Sp),
                    shooting_residual=mpstr(res_p),
                    shooting_iterations=it_p,
                )
            )
        )
        rows.append(
            asdict(
                RawSide(
                    geometry_id=geometry_id,
                    rho=mpstr(rho),
                    x=mpstr(x),
                    direction="minus",
                    side_length=mpstr(side_m),
                    jacobi_cprime=mpstr(cprime_minus),
                    normalized_S=mpstr(Sm),
                    shooting_residual=mpstr(res_m),
                    shooting_iterations=it_m,
                )
            )
        )

        Seven = (Sp + Sm) / 2
        z = mp.sin(x/2)
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
                "geometry_id": geometry_id,
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

    return geometry_id, rho_text, rows


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--rhos",
        type=parse_list,
        default=parse_list("0.015,0.018,0.022,0.027,0.033,0.040"),
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
        "--max-rational-denominator",
        type=int,
        default=100_000_000,
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("weight6_recovery_outputs"),
    )
    return p.parse_args()


def main():
    args = parse_args()
    if args.dps < 50:
        raise ValueError("Use at least 50 decimal digits.")

    mp.mp.dps = args.dps

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    print("WEIGHT-SIX HIGH-PRECISION TORUS RECOVERY", flush=True)
    print("="*72, flush=True)
    print(
        f"geometries={len(DESIGN_GEOMETRIES)}, "
        f"radii={len(args.rhos)}, angles={len(args.xs)}, "
        f"dps={args.dps}, workers={args.workers}",
        flush=True,
    )
    print(
        "Fit geometries: 1--7 (21 equations). "
        "Holdout geometries: 8--10 (9 equations).",
        flush=True,
    )
    print()

    payloads = []
    for gid, geometry in enumerate(DESIGN_GEOMETRIES, start=1):
        for rho in args.rhos:
            payloads.append(
                (
                    gid,
                    geometry,
                    rho,
                    args.xs,
                    args.dps,
                    args.ode_degree,
                    args.ode_tol_digits,
                )
            )

    all_rows = []

    with ProcessPoolExecutor(
        max_workers=args.workers,
        max_tasks_per_child=1,
    ) as pool:
        futures = {
            pool.submit(evaluate_geometry_rho, payload):
            (payload[0], payload[2])
            for payload in payloads
        }

        for future in as_completed(futures):
            gid, rho = futures[future]
            _, _, rows = future.result()
            all_rows.extend(rows)
            print(
                f"completed geometry={gid:2d} rho={rho}",
                flush=True,
            )

    # Stable sort.
    def sort_key(row):
        return (
            int(row["geometry_id"]),
            mp.mpf(row["rho"]),
            mp.mpf(row["x"]),
            {"plus":0, "minus":1, "EVEN_DIAGNOSTIC":2}.get(
                row["direction"], 9
            ),
        )

    all_rows.sort(key=sort_key)

    diagnostics = [
        row for row in all_rows
        if row["direction"] == "EVEN_DIAGNOSTIC"
    ]

    # Recover one s6 value per geometry x angle.
    observations = []

    for gid, geometry in enumerate(DESIGN_GEOMETRIES, start=1):
        R = mp.mpf(geometry[0])
        a = mp.mpf(geometry[1])
        th = mp.mpf(geometry[2])
        beta = mp.mpf(geometry[3])

        inv6 = torus_weight6_invariants(R, a, th, beta)

        # K for the already-known K^3 row.
        q = R + a*mp.cos(th)
        K = mp.cos(th)/(a*q)

        for x_text in args.xs:
            x = mp.mpf(x_text)

            subset = [
                row for row in diagnostics
                if int(row["geometry_id"]) == gid
                and abs(mp.mpf(row["x"]) - x) < mp.mpf("1e-40")
            ]
            subset.sort(key=lambda row: mp.mpf(row["rho"]))

            rhos = [mp.mpf(row["rho"]) for row in subset]
            E6s = [mp.mpf(row["E6"]) for row in subset]

            fit024, rms024 = mp_least_squares(
                rhos, E6s, (0,2,4)
            )
            fit02, rms02 = mp_least_squares(
                rhos, E6s, (0,2)
            )

            nsmall = max(4, len(rhos)-1)
            fit_small, rms_small = mp_least_squares(
                rhos[:nsmall],
                E6s[:nsmall],
                (0,2,4),
            )

            s6 = fit024[0]
            spread = abs(s6 - fit_small[0])

            z = mp.sin(x/2)
            rang = mp.cos(x/2)
            t = rang**2

            denom = z*(1-t)
            if abs(denom) < mp.mpf("1e-20"):
                raise ValueError("Angular normalization too small.")

            Pk3 = (
                mp.mpf(31)/15120
                - mp.mpf(11)/280*t
                + mp.mpf(5)/112*t**2
            )

            y = s6/denom - K**3*Pk3

            row = (
                inv6
                + [v*t for v in inv6]
                + [v*t**2 for v in inv6]
            )

            observations.append(
                {
                    "geometry_id": gid,
                    "R": mpstr(R),
                    "minor_radius": mpstr(a),
                    "theta": mpstr(th),
                    "beta": mpstr(beta),
                    "x": mpstr(x),
                    "s6": mpstr(s6),
                    "s6_fit_02": mpstr(fit02[0]),
                    "nested_small": mpstr(fit_small[0]),
                    "nested_spread": mpstr(spread),
                    "fit_rms_024": mpstr(rms024),
                    "fit_rms_02": mpstr(rms02),
                    "fit_rms_nested": mpstr(rms_small),
                    "normalized_rhs": mpstr(y),
                    "_row": row,
                    "_rhs": y,
                }
            )

    if len(observations) != 30:
        raise RuntimeError(
            f"Expected 30 observations, found {len(observations)}."
        )

    # Build matrices.
    M = mp.matrix(30, 21)
    b = mp.matrix(30, 1)

    for i, obs in enumerate(observations):
        for j, value in enumerate(obs["_row"]):
            M[i,j] = value
        b[i] = obs["_rhs"]

    fit_indices = [
        i for i, obs in enumerate(observations)
        if int(obs["geometry_id"]) <= 7
    ]
    hold_indices = [
        i for i, obs in enumerate(observations)
        if int(obs["geometry_id"]) >= 8
    ]

    if len(fit_indices) != 21 or len(hold_indices) != 9:
        raise RuntimeError("Unexpected fit/holdout partition.")

    Mfit = mp.matrix(21,21)
    bfit = mp.matrix(21,1)
    for ii, i in enumerate(fit_indices):
        for j in range(21):
            Mfit[ii,j] = M[i,j]
        bfit[ii] = b[i]

    Mhold = mp.matrix(9,21)
    bhold = mp.matrix(9,1)
    for ii, i in enumerate(hold_indices):
        for j in range(21):
            Mhold[ii,j] = M[i,j]
        bhold[ii] = b[i]

    x_fit = solve_column_scaled_square(Mfit, bfit)
    x_all = solve_column_scaled_ls(M, b)

    fit_res, fit_rms, fit_max = residual_stats(
        Mfit, x_fit, bfit
    )
    hold_res, hold_rms, hold_max = residual_stats(
        Mhold, x_fit, bhold
    )
    all_res_fit, all_rms_fit, all_max_fit = residual_stats(
        M, x_fit, b
    )
    all_res_ls, all_rms_ls, all_max_ls = residual_stats(
        M, x_all, b
    )

    coeff_difference = [
        x_fit[j] - x_all[j]
        for j in range(21)
    ]
    coeff_diff_max = max(abs(v) for v in coeff_difference)

    invariant_names = [
        "K DeltaK",
        "K K_ww",
        "|grad K|^2",
        "K_w^2",
        "S_wwww",
        "(tr S)_ww",
        "tr^2 S",
    ]

    constant_rows = []

    for block, prefix in enumerate(("A","B","C")):
        for j, inv_name in enumerate(invariant_names):
            idx = block*7 + j
            value = x_fit[idx]
            value_all = x_all[idx]
            frac, frac_err = rational_candidate(
                value,
                args.max_rational_denominator,
            )
            constant_rows.append(
                {
                    "constant": f"{prefix}[{inv_name}]",
                    "fit21_value": mpstr(value),
                    "all30_ls_value": mpstr(value_all),
                    "fit_minus_all30": mpstr(value-value_all),
                    "rational_candidate":
                        f"{frac.numerator}/{frac.denominator}",
                    "rational_error": mpstr(frac_err),
                }
            )

    # Add prediction residual to each observation.
    pred_fit = M*x_fit
    pred_all = M*x_all

    for i, obs in enumerate(observations):
        obs["prediction_fit21"] = mpstr(pred_fit[i])
        obs["residual_fit21"] = mpstr(pred_fit[i]-b[i])
        obs["prediction_all30"] = mpstr(pred_all[i])
        obs["residual_all30"] = mpstr(pred_all[i]-b[i])
        obs["role"] = (
            "FIT"
            if int(obs["geometry_id"]) <= 7
            else "HOLDOUT"
        )
        del obs["_row"]
        del obs["_rhs"]

    # Raw-data CSV.
    raw_path = out / "raw_geodesic_data.csv"
    raw_keys = []
    for row in all_rows:
        for key in row:
            if key not in raw_keys:
                raw_keys.append(key)

    with raw_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=raw_keys,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(all_rows)

    # Observation CSV.
    obs_path = out / "s6_observations.csv"
    obs_keys = list(observations[0].keys())
    with obs_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=obs_keys)
        writer.writeheader()
        writer.writerows(observations)

    # Constants CSV.
    const_path = out / "recovered_constants.csv"
    with const_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(constant_rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(constant_rows)

    summary = {
        "numerics": {
            "dps": args.dps,
            "ode_degree": args.ode_degree,
            "ode_tol_digits": args.ode_tol_digits,
            "rhos": args.rhos,
            "xs": args.xs,
        },
        "partition": {
            "fit_geometries": "1-7",
            "holdout_geometries": "8-10",
            "fit_equations": 21,
            "holdout_equations": 9,
        },
        "residuals": {
            "fit21_fit_rms": mpstr(fit_rms),
            "fit21_fit_maxabs": mpstr(fit_max),
            "fit21_holdout_rms": mpstr(hold_rms),
            "fit21_holdout_maxabs": mpstr(hold_max),
            "fit21_all30_rms": mpstr(all_rms_fit),
            "fit21_all30_maxabs": mpstr(all_max_fit),
            "all30_ls_rms": mpstr(all_rms_ls),
            "all30_ls_maxabs": mpstr(all_max_ls),
            "max_coefficient_difference_fit21_vs_all30":
                mpstr(coeff_diff_max),
        },
        "constants": constant_rows,
    }

    json_path = out / "recovery_summary.json"
    json_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    lines = [
        "WEIGHT-SIX HIGH-PRECISION TORUS RECOVERY SUMMARY",
        "="*72,
        (
            f"dps={args.dps}, ode_degree={args.ode_degree}, "
            f"ode_tol_digits={args.ode_tol_digits}"
        ),
        "",
        "Recovery partition:",
        "  geometries 1--7 : FIT (21 equations)",
        "  geometries 8--10: HOLDOUT (9 equations)",
        "",
        "Residual diagnostics:",
        f"  fit-system RMS       = {mpstr(fit_rms)}",
        f"  fit-system max abs   = {mpstr(fit_max)}",
        f"  HOLDOUT RMS          = {mpstr(hold_rms)}",
        f"  HOLDOUT max abs      = {mpstr(hold_max)}",
        f"  all-30 RMS (fit21)   = {mpstr(all_rms_fit)}",
        f"  all-30 RMS (LS)      = {mpstr(all_rms_ls)}",
        (
            "  max coeff |fit21-all30| = "
            f"{mpstr(coeff_diff_max)}"
        ),
        "",
        "Recovered constants:",
    ]

    for row in constant_rows:
        lines.extend(
            [
                f"  {row['constant']}",
                f"    fit21 = {row['fit21_value']}",
                f"    all30 = {row['all30_ls_value']}",
                (
                    "    rational candidate = "
                    f"{row['rational_candidate']}"
                ),
                (
                    "    rational error     = "
                    f"{row['rational_error']}"
                ),
            ]
        )

    lines.extend(
        [
            "",
            "Interpretation rule:",
            (
                "A convincing recovery requires small holdout residuals, "
                "close fit21/all30 constants, and stable rational "
                "reconstruction well above the numerical uncertainty."
            ),
        ]
    )

    txt_path = out / "recovery_summary.txt"
    txt_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print()
    print("\n".join(lines), flush=True)
    print()
    print("Wrote:", raw_path, flush=True)
    print("Wrote:", obs_path, flush=True)
    print("Wrote:", const_path, flush=True)
    print("Wrote:", json_path, flush=True)
    print("Wrote:", txt_path, flush=True)


if __name__ == "__main__":
    main()
