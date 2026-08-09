#!/usr/bin/env python3
"""
Reprocess weight-six raw torus data with a higher even extrapolant.

Requires in the same folder:
  raw_geodesic_data.csv
  weight6_torus_full_recovery.py

Uses:
  E6(rho) = s6 + c2 rho^2 + c4 rho^4 + c6 rho^6 + c8 rho^8

This avoids rerunning the expensive intrinsic geodesic solver.
"""
from pathlib import Path
import importlib.util, sys, pandas as pd, mpmath as mp
from fractions import Fraction

mp.mp.dps = 80
HERE = Path(__file__).resolve().parent
RAW = HERE / "raw_geodesic_data.csv"
if not RAW.exists():
    RAW = HERE.parent / "data" / "raw_geodesic_data.csv"
HELPER = HERE / "weight6_torus_full_recovery.py"

spec = importlib.util.spec_from_file_location("w6helper", str(HELPER))
w6 = importlib.util.module_from_spec(spec)
sys.modules["w6helper"] = w6
spec.loader.exec_module(w6)

raw = pd.read_csv(RAW, dtype=str)
angles = ["0.4","0.8","1.2"]
powers = (0,2,4,6,8)

def get_s6(gid, x_text):
    sub = raw[(raw.geometry_id.astype(int)==gid) & (raw.direction=="EVEN_DIAGNOSTIC")].copy()
    sub = sub[(sub.x.astype(float)-float(x_text)).abs() < 1e-12]
    sub = sub.sort_values("rho", key=lambda s:s.astype(float))
    rhos = [mp.mpf(v) for v in sub.rho]
    e6 = [mp.mpf(v) for v in sub.E6]
    coeff, rms = w6.mp_least_squares(rhos, e6, powers)
    return coeff[0], rms

rows=[]; rhs=[]
for gid,geom in enumerate(w6.DESIGN_GEOMETRIES,start=1):
    R,a,th,beta = map(mp.mpf,geom)
    inv = w6.torus_weight6_invariants(R,a,th,beta)
    q = R+a*mp.cos(th)
    K = mp.cos(th)/(a*q)
    for xt in angles:
        x=mp.mpf(xt); z=mp.sin(x/2); r=mp.cos(x/2); t=r*r
        s6,_ = get_s6(gid,xt)
        Pk3 = mp.mpf(31)/15120 - mp.mpf(11)/280*t + mp.mpf(5)/112*t*t
        y = s6/(z*(1-t)) - K**3*Pk3
        rows.append(inv+[v*t for v in inv]+[v*t*t for v in inv])
        rhs.append(y)

M=mp.matrix(30,21); b=mp.matrix(30,1)
for i in range(30):
    for j in range(21): M[i,j]=rows[i][j]
    b[i]=rhs[i]

Mfit=mp.matrix(21,21); bfit=mp.matrix(21,1)
Mhold=mp.matrix(9,21); bhold=mp.matrix(9,1)
for i in range(21):
    for j in range(21): Mfit[i,j]=M[i,j]
    bfit[i]=b[i]
for ii,i in enumerate(range(21,30)):
    for j in range(21): Mhold[ii,j]=M[i,j]
    bhold[ii]=b[i]

xf=w6.solve_column_scaled_square(Mfit,bfit)
xa=w6.solve_column_scaled_ls(M,b)
_,_,_=w6.residual_stats(Mfit,xf,bfit)
_,hrms,hmax=w6.residual_stats(Mhold,xf,bhold)
_,arms,amax=w6.residual_stats(M,xa,b)
diff=max(abs(xf[j]-xa[j]) for j in range(21))

print("WEIGHT-SIX HIGH-ORDER REPROCESSING")
print("="*64)
print("holdout RMS       =", mp.nstr(hrms,30))
print("holdout max abs   =", mp.nstr(hmax,30))
print("all-30 LS RMS     =", mp.nstr(arms,30))
print("max |fit21-all30| =", mp.nstr(diff,30))
