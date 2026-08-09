#!/usr/bin/env python3
"""
Exact symbolic audit for:
    A Curvature-Deformed Arcsine Coordinate for Cubic Geometric Refinement

Checks:
  1. Reconstruct the normalized signed side law from the frozen
     fifth-order curvature formulas.
  2. Invert it perturbatively through rho^5 and verify Psi_2,...,Psi_5.
  3. Expand in z and verify gamma_3, gamma_5, gamma_7.
  4. Verify the exact spherical deformed-arcsine coefficients.
  5. Verify the spherical tripling map through z^7.

All algebra is exact SymPy rational arithmetic.
"""

import sympy as sp

z, rho = sp.symbols("z rho")
K, Kw, D, Hww, Dw, Twww = sp.symbols("K Kw D Hww Dw Twww")
u = sp.symbols("u")
s = sp.sin(u)
c = sp.cos(u)
r = sp.sqrt(1 - z**2)

# ---------------------------------------------------------------------
# Frozen signed-side law:
# c_rho,w(x) = 2 rho sin(x/2) [1 + f2 rho^2 + f3 rho^3
#                              + B4 rho^4 + B5 rho^5 + O(rho^6)]
# with u=x/2.
# ---------------------------------------------------------------------

f2 = -K*c**2 / 6
f3 = -Kw*c**3 / 12

B4 = c**2 / 120 * (
    (1 - 9*s**2)*K**2
    - 3*c**2*Hww
    - s**2*(D - Hww)
)

Tsum = 3*Dw - 3*Twww + 2*K*Kw

B5 = c**3 / 1080 * (
    (9*c**2 - 92*s**2)*K*Kw
    - 6*c**2*Twww
    - 2*s**2*Tsum
)

F = 1 + f2*rho**2 + f3*rho**3 + B4*rho**4 + B5*rho**5
F0 = sp.expand(F.subs(u, 0))

# Normalized side law:
# H(2u) = 2 sin(u) F(u)/F(0).
G = sp.series(F / F0, rho, 0, 6).removeO().expand()

# ---------------------------------------------------------------------
# Perturbative inverse:
# Psi(z) = u = asin(z) + rho^2 U2 + ... + rho^5 U5.
# ---------------------------------------------------------------------

U2, U3, U4, U5 = sp.symbols("U2 U3 U4 U5")
delta = U2*rho**2 + U3*rho**3 + U4*rho**4 + U5*rho**5

sin_delta = sp.series(sp.sin(delta), rho, 0, 6).removeO()
cos_delta = sp.series(sp.cos(delta), rho, 0, 6).removeO()

# sin(asin z + delta), cos(asin z + delta)
sin_u = z*cos_delta + r*sin_delta
cos_u = r*cos_delta - z*sin_delta

S, C = sp.symbols("S C")
G_SC = G.subs({sp.sin(u): S, sp.cos(u): C})
G_shift = sp.series(G_SC.subs({S: sin_u, C: cos_u}), rho, 0, 6).removeO().expand()
encoded = sp.series(sin_u * G_shift, rho, 0, 6).removeO().expand()

solutions = {}
for n, U in zip(range(2, 6), (U2, U3, U4, U5)):
    equation = sp.expand(encoded - z).coeff(rho, n).subs(solutions)
    solutions[U] = sp.factor(sp.solve(sp.Eq(equation, 0), U)[0])

Psi2_derived = solutions[U2]
Psi3_derived = solutions[U3]
Psi4_derived = solutions[U4]
Psi5_derived = solutions[U5]

Psi2 = -K*z**3/(6*r)

Psi3 = -Kw*z*(1-r**3)/(12*r)

Psi4 = z**3/(360*r**3) * (
    K**2*(2*r**4 + 13*r**2 + 5)
    + 3*D*r**4
    - 3*Hww*r**2*(4*r**2 + 3)
)

Psi5 = z/(360*r**3) * (
    2*Dw*r**5*z**2
    + K*Kw*(1-r) * (
        5*r**6 + 5*r**5 + 8*r**4 + 18*r**3
        + 8*r**2 + 5*r + 5
    )
    - 2*Twww*r**2*(1-r) * (
        2*r**4 + 2*r**3 + r**2 + r + 1
    )
)

checks = {}

checks["Psi2"] = sp.simplify(Psi2_derived - Psi2)
checks["Psi3"] = sp.simplify(Psi3_derived - Psi3)
checks["Psi4"] = sp.simplify(Psi4_derived - Psi4)
checks["Psi5"] = sp.simplify(Psi5_derived - Psi5)

# ---------------------------------------------------------------------
# Taylor coefficients in z.
# ---------------------------------------------------------------------

Psi_total = (
    sp.asin(z)
    + rho**2*Psi2
    + rho**3*Psi3
    + rho**4*Psi4
    + rho**5*Psi5
)

Psi_z = sp.series(Psi_total, z, 0, 9).removeO().expand()

gamma3_derived = sp.factor(Psi_z.coeff(z, 3))
gamma5_derived = sp.factor(Psi_z.coeff(z, 5))
gamma7_derived = sp.factor(Psi_z.coeff(z, 7))

gamma3 = (
    sp.Rational(1, 6)
    - sp.Rational(1, 6)*K*rho**2
    - sp.Rational(1, 8)*Kw*rho**3
    + (
        sp.Rational(1, 18)*K**2
        + sp.Rational(1, 120)*D
        - sp.Rational(7, 120)*Hww
    )*rho**4
    + (
        sp.Rational(3, 40)*K*Kw
        + sp.Rational(1, 180)*Dw
        - sp.Rational(7, 360)*Twww
    )*rho**5
)

gamma5 = (
    sp.Rational(3, 40)
    - sp.Rational(1, 12)*K*rho**2
    - sp.Rational(1, 32)*Kw*rho**3
    + (
        sp.Rational(13, 360)*K**2
        - sp.Rational(1, 240)*D
        + sp.Rational(1, 240)*Hww
    )*rho**4
    + (
        sp.Rational(3, 160)*K*Kw
        - sp.Rational(1, 180)*Dw
        + sp.Rational(13, 1440)*Twww
    )*rho**5
)

gamma7 = (
    sp.Rational(5, 112)
    - sp.Rational(1, 16)*K*rho**2
    - sp.Rational(5, 192)*Kw*rho**3
    + (
        sp.Rational(7, 180)*K**2
        - sp.Rational(1, 960)*D
        - sp.Rational(1, 192)*Hww
    )*rho**4
    + (
        sp.Rational(5, 144)*K*Kw
        - sp.Rational(1, 576)*Twww
    )*rho**5
)

checks["gamma3"] = sp.simplify(gamma3_derived - gamma3)
checks["gamma5"] = sp.simplify(gamma5_derived - gamma5)
checks["gamma7"] = sp.simplify(gamma7_derived - gamma7)

# ---------------------------------------------------------------------
# Exact spherical benchmark.
# ---------------------------------------------------------------------

a = sp.symbols("alpha", real=True)
sa = sp.sin(a)
ca2 = sp.cos(a)**2

Psi_sphere = sp.asin(sp.sin(z*sa)/sa)
sphere_series = sp.series(Psi_sphere, z, 0, 9).removeO().expand()

sphere_g3 = sp.factor(sphere_series.coeff(z, 3))
sphere_g5 = sp.factor(sphere_series.coeff(z, 5))
sphere_g7 = sp.factor(sphere_series.coeff(z, 7))

sphere_g3_expected = ca2 / 6
sphere_g5_expected = ca2*(ca2 + 8) / 120
sphere_g7_expected = ca2*(ca2**2 + 88*ca2 + 136) / 5040

checks["sphere_gamma3"] = sp.trigsimp(sphere_g3 - sphere_g3_expected)
checks["sphere_gamma5"] = sp.trigsimp(sphere_g5 - sphere_g5_expected)
checks["sphere_gamma7"] = sp.trigsimp(sphere_g7 - sphere_g7_expected)

# Exact spherical tripling map:
# Psi^{-1}(y) = asin(sin(alpha)*sin(y))/sin(alpha).
T_sphere = sp.asin(sa*sp.sin(3*Psi_sphere))/sa
T_series = sp.series(T_sphere, z, 0, 9).removeO().expand()

T_expected = (
    3*z
    - 4*ca2*z**3
    + 16*ca2*(ca2 - 1)*z**5
    - sp.Rational(4, 15)*ca2*(ca2 - 1)*(311*ca2 - 221)*z**7
)

checks["sphere_tripling"] = sp.trigsimp(
    sp.expand(T_series - T_expected)
)

# ---------------------------------------------------------------------
# Report.
# ---------------------------------------------------------------------

failed = {name: value for name, value in checks.items() if sp.simplify(value) != 0}

print("Exact symbolic audit: curvature-deformed arcsine")
print("=" * 55)
for name, value in checks.items():
    status = "PASS" if sp.simplify(value) == 0 else "FAIL"
    print(f"{name:20s} {status}")

print()
if failed:
    print("FAILED CHECKS:")
    for name, value in failed.items():
        print(name, ":", value)
    raise SystemExit(1)

print("ALL CHECKS PASS")
print()
print("Recovered coefficients:")
print("gamma_3 =", gamma3)
print("gamma_5 =", gamma5)
print("gamma_7 =", gamma7)


# ---------------------------------------------------------------------
# General curved tripling map through z^7.
# If Psi(z)=z+g3 z^3+g5 z^5+g7 z^7+..., then
# T=Psi^{-1}(3 Psi) has coefficients below.
# ---------------------------------------------------------------------

g3, g5, g7 = sp.symbols("g3 g5 g7")
tau3 = -24*g3
tau5 = 24*(27*g3**2 - 10*g5)
tau7 = -24*(945*g3**3 - 675*g3*g5 + 91*g7)

tau3_rho = sp.series(tau3.subs(g3, gamma3), rho, 0, 6).removeO().expand()
tau5_rho = sp.series(tau5.subs({g3: gamma3, g5: gamma5}), rho, 0, 6).removeO().expand()
tau7_rho = sp.series(tau7.subs({g3: gamma3, g5: gamma5, g7: gamma7}), rho, 0, 6).removeO().expand()

tau3_expected = (
    -4
    + 4*K*rho**2
    + 3*Kw*rho**3
    - (20*K**2 + 3*D - 21*Hww)*rho**4/15
    - (27*K*Kw + 2*Dw - 7*Twww)*rho**5/15
)

tau5_expected = (
    -16*K*rho**2
    - sp.Rational(39, 2)*Kw*rho**3
    + sp.Rational(2, 15)*(160*K**2 + 21*D - 102*Hww)*rho**4
    + (1161*K*Kw + 76*Dw - 191*Twww)*rho**5/30
)

tau7_expected = (
    24*K*rho**2
    + sp.Rational(455, 8)*Kw*rho**3
    - (1724*K**2 + 219*D - 930*Hww)*rho**4/15
    - (7085*K*Kw + 450*Dw - 991*Twww)*rho**5/24
)

tripling_checks = {
    "tau3_general": sp.simplify(tau3_rho - tau3_expected),
    "tau5_general": sp.simplify(tau5_rho - tau5_expected),
    "tau7_general": sp.simplify(tau7_rho - tau7_expected),
}

print()
print("General curved tripling coefficients")
print("=" * 55)
for name, value in tripling_checks.items():
    status = "PASS" if sp.simplify(value) == 0 else "FAIL"
    print(f"{name:20s} {status}")

if any(sp.simplify(v) != 0 for v in tripling_checks.values()):
    raise SystemExit(1)

print()
print("tau_3 =", sp.collect(tau3_expected, rho))
print("tau_5 =", sp.collect(tau5_expected, rho))
print("tau_7 =", sp.collect(tau7_expected, rho))
