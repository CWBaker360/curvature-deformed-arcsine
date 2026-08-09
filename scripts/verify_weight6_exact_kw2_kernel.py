#!/usr/bin/env python3
"""
Certificate for the corrected weight-six K_w^2 angular kernel.

Requires in the same folder:
    kw2_kernel_values.csv
    angle_holdout_results.csv
    weight6_torus_full_recovery.py

No geodesic solves are performed.

Checks
------
1. Remove the normalization-forced term
       (1-r^3)/(144(1-r^2))
   from the extracted P_{K_w^2}(r).

2. Verify that the remainder is exactly modeled numerically by
       -1/504 - r^2/3360 - 317 r^4/10080.

3. Insert the corrected exact K_w^2 kernel together with the other
   18 exact rational rows and the exact spherical K^3 row.

4. Re-predict the twelve earlier out-of-sample angle measurements
   (x=0.55,0.95,1.35 on geometries 1,8,9,10).

Expected scale from the current data:
    scan corrected-kernel RMS ~ 1e-18 geometry extraction scale or better
    12-angle holdout RMS      ~ 3.5e-24
"""

from pathlib import Path
import importlib.util
import sys
import pandas as pd
import mpmath as mp

mp.mp.dps = 80
HERE = Path(__file__).resolve().parent

DATA = HERE.parent / "data"
SCAN = HERE / "kw2_kernel_values.csv"
HOLD = HERE / "angle_holdout_results.csv"
if not SCAN.exists():
    SCAN = DATA / "kw2_kernel_values.csv"
if not HOLD.exists():
    HOLD = DATA / "angle_holdout_results.csv"
HELPER = HERE / "weight6_torus_full_recovery.py"

for p in (SCAN, HOLD, HELPER):
    if not p.exists():
        raise FileNotFoundError(f"Missing required file: {p}")

spec = importlib.util.spec_from_file_location("w6helper", str(HELPER))
w6 = importlib.util.module_from_spec(spec)
sys.modules["w6helper"] = w6
spec.loader.exec_module(w6)

scan = pd.read_csv(SCAN, dtype=str)
hold = pd.read_csv(HOLD, dtype=str)

EXACT_ROWS = {
    "K DeltaK": (
        mp.mpf(0),
        -mp.mpf(1)/112,
        mp.mpf(127)/15120,
    ),
    "K K_ww": (
        mp.mpf(29)/5040,
        mp.mpf(53)/5040,
        -mp.mpf(53)/1260,
    ),
    "|grad K|^2": (
        mp.mpf(0),
        -mp.mpf(17)/10080,
        mp.mpf(71)/30240,
    ),
    "S_wwww": (
        mp.mpf(1)/1008,
        mp.mpf(1)/1260,
        mp.mpf(1)/315,
    ),
    "(tr S)_ww": (
        mp.mpf(0),
        mp.mpf(1)/2520,
        -mp.mpf(1)/420,
    ),
    "tr^2 S": (
        mp.mpf(0),
        -mp.mpf(1)/5040,
        mp.mpf(1)/5040,
    ),
}

INV_NAMES = (
    "K DeltaK",
    "K K_ww",
    "|grad K|^2",
    "K_w^2",
    "S_wwww",
    "(tr S)_ww",
    "tr^2 S",
)

def q_kw2(r):
    t = r*r
    return (
        (1-r**3)/(mp.mpf(144)*(1-r**2))
        - mp.mpf(1)/504
        - t/mp.mpf(3360)
        - mp.mpf(317)*t*t/mp.mpf(10080)
    )

def s6_kw2_polynomial_bracket(r):
    return (
        mp.mpf(5)/1008
        + mp.mpf(17)*r**2/10080
        - r**3/mp.mpf(144)
        - mp.mpf(157)*r**4/5040
        + mp.mpf(317)*r**6/10080
    )

# ------------------------------------------------------------------
# Scan check
# ------------------------------------------------------------------
scan_residuals = []

for _, row in scan.iterrows():
    r = mp.mpf(row["r"])
    observed = mp.mpf(row["P_Kw2"])
    predicted = q_kw2(r)
    scan_residuals.append(observed-predicted)

scan_rms = mp.sqrt(
    mp.fsum(v*v for v in scan_residuals)/len(scan_residuals)
)
scan_max = max(abs(v) for v in scan_residuals)

# ------------------------------------------------------------------
# Full exact s6 model
# ------------------------------------------------------------------
def predict_s6(gid, x):
    geom = w6.DESIGN_GEOMETRIES[int(gid)-1]
    R, a, th, beta = map(mp.mpf, geom)

    inv = w6.torus_weight6_invariants(R, a, th, beta)

    q = R + a*mp.cos(th)
    K = mp.cos(th)/(a*q)

    z = mp.sin(x/2)
    r = mp.cos(x/2)
    t = r*r

    Pk3 = (
        mp.mpf(31)/15120
        - mp.mpf(11)*t/280
        + mp.mpf(5)*t*t/112
    )

    total = K**3 * Pk3

    for I, name in zip(inv, INV_NAMES):
        if name == "K_w^2":
            total += I*q_kw2(r)
        else:
            A, B, C = EXACT_ROWS[name]
            total += I*(A+B*t+C*t*t)

    return z*(1-t)*total

hold_residuals = []
angle_only = []
double_hold = []

for _, row in hold.iterrows():
    gid = int(row["geometry_id"])
    x = mp.mpf(row["x"])
    obs = mp.mpf(row["s6_observed"])
    pred = predict_s6(gid, x)
    res = obs-pred

    hold_residuals.append(res)
    if gid == 1:
        angle_only.append(res)
    else:
        double_hold.append(res)

hold_rms = mp.sqrt(
    mp.fsum(v*v for v in hold_residuals)/len(hold_residuals)
)
hold_max = max(abs(v) for v in hold_residuals)

angle_rms = mp.sqrt(
    mp.fsum(v*v for v in angle_only)/len(angle_only)
)
double_rms = mp.sqrt(
    mp.fsum(v*v for v in double_hold)/len(double_hold)
)

print("WEIGHT-SIX EXACT K_w^2 KERNEL CERTIFICATE")
print("="*72)
print()
print("Corrected factored kernel:")
print("Q_Kw2(r) = (1-r^3)/(144(1-r^2)) - 1/504 - r^2/3360 - 317 r^4/10080")
print()
print("Equivalent normalized-side-law polynomial:")
print("s6[K_w^2]/(z K_w^2) =")
print("  5/1008 + 17 r^2/10080 - r^3/144 - 157 r^4/5040 + 317 r^6/10080")
print()
print("Eight-angle scan:")
print("  RMS residual     =", mp.nstr(scan_rms, 30))
print("  max |residual|   =", mp.nstr(scan_max, 30))
print()
print("Twelve independent angle holdouts:")
print("  all-12 RMS       =", mp.nstr(hold_rms, 30))
print("  all-12 max abs   =", mp.nstr(hold_max, 30))
print("  geometry-1 RMS   =", mp.nstr(angle_rms, 30))
print("  geometries 8-10  =", mp.nstr(double_rms, 30))
print()
print("PASS if these residuals remain at the extraction/roundoff scale.")
