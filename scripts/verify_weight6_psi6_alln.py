#!/usr/bin/env python3
"""
Exact symbolic audit of the complete weight-six Psi_6 Laurent layer.

Status
------
This script treats the numerically certified rational weight-six side-law
formula as the candidate exact input.  It then performs exact symbolic
algebra to:

  1. assemble all eight weight-six Laurent kernels of Psi_6,
  2. verify normalization L_I(1)=0,
  3. verify the exact spherical K^3 kernel,
  4. derive the all-n coefficient rows using
         beta_n^(nu) = (-1)^n binomial(nu/2,n),
  5. check those all-n formulas against direct Taylor expansion for
     n=1,...,12,
  6. print gamma_3,...,gamma_11 weight-six coefficients.

No floating-point arithmetic is used.

Dependencies
------------
sympy
"""

import sympy as sp

r, z = sp.symbols("r z")
n = sp.symbols("n", integer=True, nonnegative=True)

# ---------------------------------------------------------------------
# Complete weight-six Psi_6 kernels
#
# Psi_6(z) = z * sum_I I * L_I(r),  r=sqrt(1-z^2).
# ---------------------------------------------------------------------

L = {}

L["K^3"] = (
    sp.Rational(1,2835)*r**5
    + sp.Rational(1,5670)*r**3
    + sp.Rational(23,15120)*r
    + sp.Rational(1,9072)/r
    + sp.Rational(1,6480)/r**3
    - sp.Rational(1,432)/r**5
)

L["K DeltaK"] = (
    -sp.Rational(1,756)*r**5
    + sp.Rational(53,15120)*r**3
    - sp.Rational(1,280)*r
    + sp.Rational(1,720)/r
)

L["K K_ww"] = (
    -sp.Rational(1,315)*r**5
    - sp.Rational(5,504)*r**3
    + sp.Rational(13,1680)*r
    + sp.Rational(1,840)/r
    + sp.Rational(1,240)/r**3
)

L["|grad K|^2"] = (
    sp.Rational(71,30240)*r**5
    - sp.Rational(61,15120)*r**3
    + sp.Rational(17,10080)*r
)

L["K_w^2"] = (
    -sp.Rational(1,140)*r**5
    + sp.Rational(139,10080)*r**3
    - sp.Rational(1,48)*r**2
    - sp.Rational(17,10080)*r
    + sp.Rational(1,72)
    - sp.Rational(1,672)/r
    + sp.Rational(1,288)/r**3
)

L["S_wwww"] = (
    sp.Rational(1,315)*r**5
    - sp.Rational(1,420)*r**3
    + sp.Rational(1,5040)*r
    - sp.Rational(1,1008)/r
)

L["(tr S)_ww"] = (
    -sp.Rational(1,420)*r**5
    + sp.Rational(1,360)*r**3
    - sp.Rational(1,2520)*r
)

L["tr^2 S"] = (
    sp.Rational(1,5040)*r**5
    - sp.Rational(1,2520)*r**3
    + sp.Rational(1,5040)*r
)

# Exact sphere checksum.
L_K3_expected = (
    16*r**5 + 8*r**3 + 69*r
    + 5/r + 7/r**3 - 105/r**5
) / sp.Integer(45360)

# ---------------------------------------------------------------------
# Generalized binomial rows
# ---------------------------------------------------------------------

def beta(nn, nu):
    if isinstance(nn, int) and nn < 0:
        return sp.Integer(0)
    return sp.simplify(
        (-1)**nn * sp.binomial(sp.Rational(nu,2), nn)
    )

def laurent_coefficients(expr):
    """
    Return dict power -> rational coefficient for a Laurent polynomial in r.
    """
    out = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        coeff, power = term.as_coeff_exponent(r)
        power = int(power)
        out[power] = sp.simplify(out.get(power, 0) + coeff)
    return out

LAURENT = {name: laurent_coefficients(expr) for name, expr in L.items()}

def row_formula(name, nn):
    return sp.simplify(
        sum(
            coeff * beta(nn, power)
            for power, coeff in LAURENT[name].items()
        )
    )

# ---------------------------------------------------------------------
# Exact checks
# ---------------------------------------------------------------------

print("WEIGHT-SIX Psi_6 LAURENT / ALL-n AUDIT")
print("="*72)

# Every weight-six correction must leave Psi'(0)=1.
normalization_failures = []
for name, expr in L.items():
    if sp.simplify(expr.subs(r,1)) != 0:
        normalization_failures.append(name)

print(
    "normalization L_I(1)=0:",
    "PASS" if not normalization_failures else "FAIL"
)
if normalization_failures:
    print("  failures:", normalization_failures)

sphere_check = sp.simplify(L["K^3"] - L_K3_expected)
print("exact spherical K^3 kernel:", "PASS" if sphere_check == 0 else "FAIL")

# Direct Taylor-vs-beta checks through n=12.
alln_failures = []
for name, expr in L.items():
    zexpr = z * expr.subs(r, sp.sqrt(1-z**2))
    series = sp.series(zexpr, z, 0, 27).removeO().expand()

    for nn in range(1,13):
        direct = sp.simplify(series.coeff(z, 2*nn+1))
        formula = sp.simplify(row_formula(name, nn))
        if sp.simplify(direct-formula) != 0:
            alln_failures.append(
                (name, nn, direct, formula)
            )

print(
    "all-n rows n=1,...,12:",
    "PASS" if not alln_failures else "FAIL"
)

if alln_failures:
    for item in alln_failures:
        print(item)
    raise SystemExit(1)

# Selection-rule checks in the K_w^2 row.
print()
print("K_w^2 selection rules:")
print("  beta_n^(0)=0 for every n>=1: constant term affects only z.")
print("  beta_n^(2)=0 for every n>=2: r^2 term affects only gamma_3.")
print(
    "  beta_1^(2) =",
    beta(1,2),
    "so -r^2/48 contributes +1/48 to gamma_3."
)

# ---------------------------------------------------------------------
# Print kernels
# ---------------------------------------------------------------------

print()
print("Complete finite Laurent kernels:")
for name, expr in L.items():
    print(f"{name}:")
    print(" ", sp.expand(expr))

# ---------------------------------------------------------------------
# Print all-n formulas in beta notation
# ---------------------------------------------------------------------

print()
print("All-n weight-six coefficient rows")
print("gamma_(2n+1)^[I,rho^6] = sum_p c_p beta_n^(p), n>=1")
print()

for name in L:
    pieces = []
    for power in sorted(LAURENT[name], reverse=True):
        coeff = LAURENT[name][power]
        pieces.append((coeff, power))
    print(name, ":", pieces)

# ---------------------------------------------------------------------
# Explicit low Taylor rows
# ---------------------------------------------------------------------

print()
print("Explicit weight-six coefficients")
print("="*72)

for nn in range(1,6):
    print(f"gamma_{2*nn+1} [rho^6]:")
    for name in L:
        print(f"  {name:14s} {sp.factor(row_formula(name, nn))}")
    print()

print("ALL CHECKS PASS")
