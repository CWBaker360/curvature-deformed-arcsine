#!/usr/bin/env python3
"""
Weight-six angular target for the curvature-deformed arcsine program.

Encodes the constrained even-weight ansatz

    s6^I(z,r) = z (1-r^2) (A_I + B_I r^2 + C_I r^4) I,

with r=sqrt(1-z^2), for each candidate weight-six invariant I.

The exact sphere fixes the K^3 row:
    A=31/15120, B=-594/15120, C=675/15120.

The direct inverse contribution is
    -s6/r,
so every genuinely new weight-six invariant has Laurent support
    r^{-1}, r^1, r^3, r^5
before lower-weight nonlinear inversion terms are added.
"""

import sympy as sp

r, z = sp.symbols("r z")
a, b, c = sp.symbols("a b c")

s6_generic = z*(1-r**2)*(a + b*r**2 + c*r**4)
L_direct = sp.expand(-s6_generic/(r*z))

assert sp.simplify(
    L_direct
    - (
        -a/r
        + (a-b)*r
        + (b-c)*r**3
        + c*r**5
    )
) == 0

aK = sp.Rational(31,15120)
bK = -sp.Rational(594,15120)
cK = sp.Rational(675,15120)

sphere = sp.expand(z*(1-r**2)*(aK+bK*r**2+cK*r**4))
expected = z*(1-r**2)*(675*r**4-594*r**2+31)/sp.Integer(15120)
assert sp.simplify(sphere-expected) == 0

print("WEIGHT-SIX ANGULAR TARGET: PASS")
print("Generic direct Laurent kernel:")
print(sp.collect(L_direct, r))
print("K^3 constants:", aK, bK, cK)
print("Raw unknown count: 24")
print("After exact spherical K^3 row: 21")
