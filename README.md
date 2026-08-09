# A Curvature-Deformed Arcsine Coordinate for Cubic Geometric Refinement

**Author:** Wayne Baker  
**Version:** 1.0.0  
**Release date:** August 8, 2026  
**Status:** Frozen research manuscript and reproducibility archive

## Manuscript

The release manuscript is:

- `paper/baker_curvature_deformed_arcsine.tex`
- `paper/baker_curvature_deformed_arcsine.pdf`

The manuscript develops the normalized intrinsic geodesic side law as an exact local linearizing coordinate for geometric encode--scale--decode refinement. The associated half-angle inverse is a curvature-deformed arcsine. The paper contains:

- exact geometric Koenigs/Schroeder conjugacy and local semigroup structure;
- an exact spherical curvature-deformed arcsine;
- an analytic variable-curvature expansion through geometric weight five;
- a curved cubic coefficient-lifting theorem;
- a high-precision rational reconstruction of the general weight-six side law;
- exact weight-six finite Laurent kernels conditional on that reconstructed side law;
- all Taylor coefficients at weight six in generalized-binomial form; and
- a finite Laurent-kernel conjecture for the full hierarchy.

## Mathematical status boundary

The paper is **analytic through geometric weight five**.

At geometric weight six, the normalized side-law coefficient is a **high-precision rational reconstruction** obtained from intrinsic torus geodesics and tested on unused geometries and unused side angles. It is not presented as interval arithmetic and is not a substitute for the still-missing analytic degree-eight two-point geodesic-distance derivation.

Conditional on the reconstructed weight-six side law, the Laurent kernels and all-`n` Taylor formulas are exact symbolic consequences.

The exact spherical model independently fixes the full `K^3 rho^6` sector.

## Reproducibility

Python dependencies are listed in `requirements.txt`.

From the extracted release root:

```powershell
python -m pip install -r requirements.txt
```

A fast retained-data verification sequence is:

```powershell
python scripts/verify_curvature_deformed_arcsine.py
python scripts/verify_weight6_scaffold.py
python scripts/weight6_angular_target.py
python scripts/weight6_torus_design_rank.py
python scripts/weight6_reprocess_high_order.py
python scripts/verify_weight6_exact_kw2_kernel.py
python scripts/verify_weight6_psi6_alln.py
```

The high-precision intrinsic-geodesic runs are intentionally separate because they are substantially more expensive:

```powershell
python scripts/weight6_torus_s6_pilot.py
python scripts/weight6_torus_full_recovery.py
python scripts/weight6_angle_holdout_validation_fixed.py
python scripts/weight6_kw2_kernel_scan.py
```

See `VERIFICATION_QUICK_RUN.txt` for the retained-data verification record and `REFEREE_AUDIT.md` for the final manuscript/package audit.

## Key numerical checkpoints

- Torus design: 30 observations, 21 unknowns, rank 21.
- Higher-order unused-geometry holdout RMS: approximately `9.6326e-23`.
- Initial restricted new-angle ansatz: rejected, RMS approximately `1.94e-8`.
- Corrected `K_w^2` model: twelve independent angle-holdout RMS approximately `3.4600e-24`.
- Exact Laurent/all-`n` symbolic audit: all checks pass through `n=12`.

## Why `UNRESOLVED_K_W_SQUARED` appears in an intermediate CSV

`data/recovered_constants_high_order.csv` deliberately preserves an intermediate recovery stage. The three entries marked `UNRESOLVED_K_W_SQUARED` belong to the original overly restrictive quadratic angular fit. They are **not unresolved parameters in the final manuscript**.

New-angle testing falsified that restricted ansatz. The discrepancy localized to the `K_w^2` sector and was explained analytically by normalization of the known weight-three side-law term. The corrected kernel is verified by `scripts/verify_weight6_exact_kw2_kernel.py` and is the one used in the final weight-six formula.

The historical intermediate file is retained because the falsification-and-recovery path is part of the reproducibility record.

## Directory layout

```text
paper/      Frozen LaTeX source and compiled PDF
scripts/    Symbolic and high-precision verification/recovery programs
data/       Retained raw data, extrapolated observations, and validation outputs
```

Additional release files:

- `CITATION.cff` - citation metadata
- `LICENSE.md` - dual-license notice
- `RELEASE_NOTES.md` - version 1.0.0 release notes
- `SUBMISSION_METADATA.md` - target-neutral journal submission metadata
- `RELEASE_VERIFICATION.txt` - fresh quick verification run from the v1.0.0 release tree
- `MANIFEST.txt` - distributed file list
- `SHA256SUMS.txt` - SHA-256 hashes for release files
- `VERSION` - release version

## License

The manuscript and documentation are released under **Creative Commons Attribution 4.0 International (CC BY 4.0)**. Python source code is released under the **MIT License**. See `LICENSE.md`.

## Outstanding analytic problem

The principal unresolved proof obligation is an independent analytic derivation of the degree-eight two-point geodesic-distance term sufficient to prove the reconstructed general weight-six side law without numerical recovery.
