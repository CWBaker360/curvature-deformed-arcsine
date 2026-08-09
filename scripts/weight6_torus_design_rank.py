#!/usr/bin/env python3
"""
Weight-six torus design-rank audit.

This script checks whether a family of ring tori can identify the 21
remaining angular constants in the candidate weight-six side-law ansatz

    s6(z) = z(1-rang^2) sum_I I [A_I + B_I rang^2 + C_I rang^4],

after the K^3 row has already been fixed by the exact sphere.

The seven remaining candidate invariants are

    K DeltaK,
    K K_ww,
    |grad K|^2,
    K_w^2,
    S_wwww,
    (tr S)_ww,
    tr^2 S,

with S = Sym nabla^4 K.

The metric is the general ring torus

    ds^2 = a^2 dtheta^2 + (R + a cos theta)^2 dphi^2.

For each torus geometry and bisector direction w, the script constructs
the curvature derivatives directly from the metric, builds the 21-column
design row, and verifies that the proposed 30-observation design has rank 21.

Dependencies: sympy, numpy
"""

import itertools
import math
import numpy as np
import sympy as sp


# ---------------------------------------------------------------------
# Symbolic torus geometry
# ---------------------------------------------------------------------

theta, phi, R, a, beta = sp.symbols(
    "theta phi R a beta", positive=True, real=True
)
coords = (theta, phi)
q = R + a * sp.cos(theta)

g = sp.Matrix([
    [a**2, 0],
    [0, q**2],
])
ginv = sp.simplify(g.inv())

Gamma = [[[sp.Integer(0) for _ in range(2)] for _ in range(2)] for _ in range(2)]
for k in range(2):
    for i in range(2):
        for j in range(2):
            Gamma[k][i][j] = sp.simplify(
                sp.Rational(1, 2)
                * sum(
                    ginv[k, ell]
                    * (
                        sp.diff(g[ell, j], coords[i])
                        + sp.diff(g[ell, i], coords[j])
                        - sp.diff(g[i, j], coords[ell])
                    )
                    for ell in range(2)
                )
            )

K = sp.simplify(sp.cos(theta) / (a * q))


def covariant_derivative_covariant(T, rank):
    """
    If T is a dict keyed by rank indices, return nabla T with the new
    derivative index placed first.
    """
    out = {}
    for idx in itertools.product(range(2), repeat=rank + 1):
        i = idx[0]
        rest = idx[1:]

        expr = sp.diff(T[rest], coords[i])

        for slot in range(rank):
            j = rest[slot]
            for m in range(2):
                replaced = list(rest)
                replaced[slot] = m
                expr -= Gamma[m][i][j] * T[tuple(replaced)]

        out[idx] = sp.simplify(expr)

    return out


T1 = {(i,): sp.diff(K, coords[i]) for i in range(2)}
T2 = covariant_derivative_covariant(T1, 1)
T3 = covariant_derivative_covariant(T2, 2)
T4 = covariant_derivative_covariant(T3, 3)

# Fully symmetrized fourth derivative.
S4 = {}
for idx in itertools.product(range(2), repeat=4):
    perms = set(itertools.permutations(idx))
    S4[idx] = sp.simplify(
        sum(T4[p] for p in perms) / len(perms)
    )

# Unit direction in the local orthonormal frame:
# w = cos(beta) e_theta + sin(beta) e_phi.
w = (
    sp.cos(beta) / a,
    sp.sin(beta) / q,
)

Kw = sp.simplify(
    sum(T1[(i,)] * w[i] for i in range(2))
)

Kww = sp.simplify(
    sum(
        T2[(i, j)] * w[i] * w[j]
        for i in range(2)
        for j in range(2)
    )
)

DeltaK = sp.simplify(
    sum(
        ginv[i, j] * T2[(i, j)]
        for i in range(2)
        for j in range(2)
    )
)

gradK2 = sp.simplify(
    sum(
        ginv[i, j] * T1[(i,)] * T1[(j,)]
        for i in range(2)
        for j in range(2)
    )
)

Swwww = sp.simplify(
    sum(
        S4[(i, j, k, ell)]
        * w[i] * w[j] * w[k] * w[ell]
        for i in range(2)
        for j in range(2)
        for k in range(2)
        for ell in range(2)
    )
)

trS = {}
for k in range(2):
    for ell in range(2):
        trS[(k, ell)] = sp.simplify(
            sum(
                ginv[i, j] * S4[(i, j, k, ell)]
                for i in range(2)
                for j in range(2)
            )
        )

trSww = sp.simplify(
    sum(
        trS[(k, ell)] * w[k] * w[ell]
        for k in range(2)
        for ell in range(2)
    )
)

tr2S = sp.simplify(
    sum(
        ginv[k, ell] * trS[(k, ell)]
        for k in range(2)
        for ell in range(2)
    )
)

invariant_names = [
    "K DeltaK",
    "K K_ww",
    "|grad K|^2",
    "K_w^2",
    "S_wwww",
    "(tr S)_ww",
    "tr^2 S",
]

invariant_exprs = [
    sp.simplify(K * DeltaK),
    sp.simplify(K * Kww),
    gradK2,
    sp.simplify(Kw**2),
    Swwww,
    trSww,
    tr2S,
]

inv_fun = sp.lambdify(
    (R, a, theta, beta),
    invariant_exprs,
    modules="numpy",
)


# ---------------------------------------------------------------------
# Design matrix
# ---------------------------------------------------------------------

# Side angles x; rang = cos(x/2).
side_angles = (0.4, 0.8, 1.2)

# Ten geometry-direction configurations.
# Tuple = (major radius R, minor radius a, latitude theta, bisector beta).
design_geometries = [
    (3.0, 1.0, 1.1, 0.0),
    (3.0, 0.7, 2.0, 1.4),
    (5.0, 0.7, 2.5, 0.0),
    (3.0, 0.7, 0.7, 1.4),
    (3.0, 0.7, 2.0, 1.05),
    (3.0, 0.7, 0.3, 0.0),
    (3.0, 0.7, 2.5, 1.4),
    (4.0, 0.7, 0.3, 0.35),
    (3.0, 1.0, 2.0, 0.35),
    (4.0, 1.0, 1.5, 0.0),
]


def design_row(Rv, av, thetav, betav, x):
    vals = np.asarray(
        inv_fun(Rv, av, thetav, betav),
        dtype=np.float64,
    ).reshape(7)

    rang = math.cos(x / 2.0)
    t = rang * rang

    # Unknown ordering:
    # A_1,...,A_7, B_1,...,B_7, C_1,...,C_7.
    return np.concatenate(
        [
            vals,
            vals * t,
            vals * t * t,
        ]
    )


rows = []
labels = []

for geometry in design_geometries:
    for x in side_angles:
        rows.append(design_row(*geometry, x))
        labels.append((*geometry, x))

M = np.asarray(rows, dtype=np.float64)

rank_raw = np.linalg.matrix_rank(M)

# Column normalization makes the condition number meaningful despite
# different dimensional scales of the invariants.
column_norms = np.linalg.norm(M, axis=0)
if np.any(column_norms == 0.0):
    raise RuntimeError("Zero design column encountered.")

Mscaled = M / column_norms
rank_scaled = np.linalg.matrix_rank(Mscaled)
condition_scaled = np.linalg.cond(Mscaled)

print("WEIGHT-SIX TORUS DESIGN-RANK AUDIT")
print("=" * 64)
print(f"observations                 : {M.shape[0]}")
print(f"unknown constants            : {M.shape[1]}")
print(f"raw matrix rank              : {rank_raw}")
print(f"column-normalized matrix rank: {rank_scaled}")
print(f"scaled condition number      : {condition_scaled:.12g}")
print()

if rank_scaled != 21:
    raise SystemExit("FAIL: proposed torus design is not full rank.")

print("FULL RANK 21: PASS")
print("Redundant equations:", M.shape[0] - 21)
print()

print("Invariant rows:")
for name in invariant_names:
    print(" ", name)

print()
print("Recommended geometry groups")
print("(R, minor_radius, theta, beta)")
for i, geometry in enumerate(design_geometries, start=1):
    print(f"{i:2d}: {geometry}")

print()
print("Side angles for each geometry:", side_angles)
print()
print("Interpretation:")
print("30 measured s6 values determine 21 remaining constants,")
print("leaving 9 redundant equations for an internal consistency test.")
