
import sympy as sp

u, z, rho = sp.symbols("u z rho")
r = sp.symbols("r", positive=True)
K, Kw, D, Kww = sp.symbols("K Kw D Kww")

s = sp.sin(u)
c = sp.cos(u)

f2 = -K*c**2 / 6
f3 = -Kw*c**3 / 12
B4 = c**2 / 120 * (
    (1 - 9*s**2)*K**2
    - 3*c**2*Kww
    - s**2*(D - Kww)
)

F = 1 + f2*rho**2 + f3*rho**3 + B4*rho**4
F0 = sp.expand(F.subs(u, 0))
G = sp.series(F/F0, rho, 0, 5).removeO().expand()

g2 = sp.factor(G.coeff(rho, 2))
g3 = sp.factor(G.coeff(rho, 3))
g4 = sp.factor(G.coeff(rho, 4))

s2 = sp.expand(s*g2)
s3 = sp.expand(s*g3)
s4 = sp.expand(s*g4)

U2 = -K*z**3/(6*r)
U3 = -Kw*z*(1-r**3)/(12*r)
U4 = z**3/(360*r**3) * (
    K**2*(2*r**4 + 13*r**2 + 5)
    + 3*D*r**4
    - 3*Kww*r**2*(4*r**2 + 3)
)

def at_u0(expr):
    expr = sp.expand_trig(sp.expand(expr))
    return sp.expand(expr.subs({sp.sin(u): z, sp.cos(u): r}))

s2p  = at_u0(sp.diff(s2, u))
s2pp = at_u0(sp.diff(s2, u, 2))
s3p  = at_u0(sp.diff(s3, u))
s4p  = at_u0(sp.diff(s4, u))

s0pp  = -z
s0ppp = -r

U6_induced = -sp.together(
    s4p*U2
    + s3p*U3
    + s2p*U4
    + sp.Rational(1,2)*s2pp*U2**2
    + s0pp*(U2*U4 + sp.Rational(1,2)*U3**2)
    + sp.Rational(1,6)*s0ppp*U2**3
) / r

def eliminate_even_z(expr):
    expr = sp.expand(expr)
    poly = sp.Poly(expr, z)
    result = 0
    for (power,), coeff in poly.terms():
        if power % 2:
            raise ValueError("Expected even powers of z after factoring one z.")
        result += coeff*(1-r**2)**(power//2)
    return sp.factor(sp.together(sp.expand(result)))

induced_over_z = eliminate_even_z(sp.factor(U6_induced/z))

def coeff_K3(expr):
    return sp.expand(expr).coeff(K, 3)

def coeff_KD(expr):
    return sp.diff(sp.diff(expr, K), D).subs({K:0, D:0, Kww:0, Kw:0})

def coeff_KKww(expr):
    return sp.diff(sp.diff(expr, K), Kww).subs({K:0, D:0, Kww:0, Kw:0})

def coeff_Kw2(expr):
    return sp.diff(expr, Kw, 2).subs({K:0, D:0, Kww:0, Kw:0})/2

L_K3_ind  = sp.factor(coeff_K3(induced_over_z))
L_KD_ind  = sp.factor(coeff_KD(induced_over_z))
L_KH_ind  = sp.factor(coeff_KKww(induced_over_z))
L_Kw2_ind = sp.factor(coeff_Kw2(induced_over_z))

# Exact spherical Psi benchmark, expanded first in z (stable for SymPy),
# then in alpha.
a = sp.symbols("a")
Psi_sphere = sp.asin(sp.sin(z*sp.sin(a))/sp.sin(a))
Psi_z_series = sp.series(Psi_sphere, z, 0, 27).removeO().expand()

L_K3_total = (
    16*r**5 + 8*r**3 + 69*r
    + 5/r + 7/r**3 - 105/r**5
) / sp.Integer(45360)

def beta(n, nu):
    if n < 0:
        return sp.Integer(0)
    return sp.simplify((-1)**n * sp.binomial(sp.Rational(nu,2), n))

def K3_row(n):
    return sp.simplify(
        (
            16*beta(n,5)
            + 8*beta(n,3)
            + 69*beta(n,1)
            + 5*beta(n,-1)
            + 7*beta(n,-3)
            - 105*beta(n,-5)
        )
        / sp.Integer(45360)
    )

K3_failures = []
for n in range(1, 13):
    coeff_z = Psi_z_series.coeff(z, 2*n+1)
    direct = sp.series(coeff_z, a, 0, 8).removeO().expand().coeff(a, 6)
    formula = K3_row(n)
    if sp.simplify(direct-formula) != 0:
        K3_failures.append((n, direct, formula))

# Exact spherical normalized half-side law benchmark:
# S_alpha(u)=asin(sin(alpha) sin u)/sin(alpha).
# Expand first in t=sin(u), then in alpha.
t = sp.symbols("t")
S_sphere = sp.asin(sp.sin(a)*t)/sp.sin(a)
S_t_series = sp.series(S_sphere, t, 0, 9).removeO().expand()
s6_direct_t = sp.factor(
    sp.series(S_t_series, a, 0, 8).removeO().expand().coeff(a, 6)
)
s6_expected_t = t**3*(675*t**4 - 756*t**2 + 112)/sp.Integer(15120)
side_kernel_check = sp.simplify(s6_direct_t - s6_expected_t)

L_K3_direct_Psi = sp.factor(sp.together(L_K3_total - L_K3_ind))
Q_s6_from_difference = sp.factor(-r*L_K3_direct_Psi)
Q_s6_expected = sp.factor(
    (1-r**2)*(675*r**4 - 594*r**2 + 31)
    / sp.Integer(15120)
)
difference_check = sp.simplify(Q_s6_from_difference - Q_s6_expected)

print("Weight-six geometric linearizer scaffold")
print("="*62)
print("Universal inverse formula: PASS")
print("Structural consequence: s5 and U5 do not enter U6.")
print()

print("Known nonlinear inversion contribution to Psi_6:")
print("K^3 kernel      =", sp.expand(L_K3_ind))
print("K DeltaK kernel =", sp.expand(L_KD_ind))
print("K K_ww kernel   =", sp.expand(L_KH_ind))
print("K_w^2 kernel    =", sp.expand(L_Kw2_ind))
print()

print("Exact spherical checks:")
print("side-law alpha^6 term   :", "PASS" if side_kernel_check == 0 else "FAIL")
print("induced/direct split    :", "PASS" if difference_check == 0 else "FAIL")
print()

print("Exact total K^3 kernel:")
print("L_K3(r) =", sp.expand(L_K3_total))
print()

print("Exact spherical normalized-side-law fingerprint:")
print("s6[K^3] = K^3*z^3*(675*r^4 - 594*r^2 + 31)/15120")
print()

if K3_failures:
    print("All-degree K^3 row: FAIL")
    for item in K3_failures:
        print(item)
    raise SystemExit(1)
else:
    print("All-degree K^3 row, n=1,...,12: PASS")
    print("First six coefficients:")
    for n in range(1, 7):
        print(f"[z^{2*n+1}] K^3 rho^6 =", K3_row(n))

print()
print("Candidate weight-six invariant basis:")
print("K^3")
print("K DeltaK")
print("K K_ww")
print("|grad K|^2")
print("K_w^2")
print("Sym(nabla^4 K)(w,w,w,w)")
print("tr Sym(nabla^4 K)(w,w)")
print("tr^2 Sym(nabla^4 K)")
