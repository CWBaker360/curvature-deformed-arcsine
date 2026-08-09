# Release Notes - Version 1.0.0

**Release date:** August 8, 2026

This is the first frozen release of *A Curvature-Deformed Arcsine Coordinate for Cubic Geometric Refinement* and its reproducibility archive.

## Frozen mathematical scope

- Exact normalized-side-law linearization and local multiplicative semigroup.
- Exact spherical curvature-deformed arcsine.
- Analytic variable-curvature expansion through geometric weight five.
- Curved cubic coefficient-lifting theorem.
- High-precision rational reconstruction of the weight-six normalized side law.
- Exact symbolic weight-six Laurent kernels conditional on the reconstructed side law.
- Exact generalized-binomial all-`n` formulas at weight six.
- Finite Laurent-kernel conjecture.

## Final validation state

- Manuscript compiles cleanly with no undefined references and no overfull/underfull boxes in the audited build.
- Exact spherical `K^3 rho^6` checksum passes.
- Retained-data unused-geometry holdout RMS: approximately `9.6326e-23`.
- Corrected independent angle-holdout RMS: approximately `3.4600e-24`.
- Exact Laurent/all-`n` audit: all checks pass through `n=12`.
- Release archive includes raw geodesic data, intermediate failed angular model, corrected `K_w^2` analysis, scripts, and verification records.

## Deliberately retained failed intermediate model

The archive preserves the original restricted quadratic angular ansatz and the intermediate `UNRESOLVED_K_W_SQUARED` labels. New-angle validation falsified that restriction. The resulting discrepancy localized to `K_w^2`, was explained by the normalization-forced `(1-r^3)` contribution, and led to the corrected final kernel. This record is retained to document the actual discovery and falsification path.

## Outstanding proof obligation

The general weight-six side law remains a high-precision reconstruction pending an independent analytic degree-eight two-point geodesic-distance derivation.
