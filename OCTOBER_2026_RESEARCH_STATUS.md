# October 2026 research-status cross-reference

**Date:** 21 October 2026  
**Applies to:** Related directional-curvature verification campaign; **not** a new proof of this repository's weight-six reconstruction.

## Release integrity

The August 8, 2026 v1.0.0 manuscript and supporting numerical reconstruction remain frozen. This informational supplement does not alter any `paper/`, `scripts/`, or `data/` artifact and does not replace the historical `SHA256SUMS.txt` manifest. Because this supplementary file postdates that manifest, it is not covered by those historical checksums.

## Related work and its verified scope

The October 20, 2026 verification closeout concerns **Directional Curvature Expansion of the Baker Cubic Coefficient on Riemannian Surfaces Through Fifth Order**, a related but distinct manuscript. Its focused v1.6 audit reports:

- A local even real-analytic chord expansion and cubic chord/Jacobi identity under *real-analytic metric*, nonvanishing angular Jacobi field, and convex normal-neighborhood hypotheses. These local analytic results do not automatically apply to every smooth metric.
- Exact-rational five-node projection and Richardson cancellation formulas, conditional on the applicable chord regularity.
- High-precision torus numerical consistency evidence, **not** a certified interval enclosure for the resulting curvature estimator.
- Separate Whitney derivative/jet compatibility checks, which do not certify geodesic quadrature, a global minimizing joining branch, or the directional-curvature tensor coefficients.

The October 20 campaign closeout explicitly retains **`release_certified = false`** for the related main directional-curvature manuscript. The independent degree-six/seven tensor transcription and full fifth-order coefficient audit remain open.

## What this does *not* establish

The above work does **not** discharge the outstanding obligation described in this repository's original README and referee audit: an independent analytic degree-eight two-point geodesic-distance derivation proving the **general weight-six side law**. That part of this repository still has the status **high-precision rational reconstruction**, with conditional exact symbolic consequences.

The October directional-curvature verification packages also have different scopes and should not be presented as replacements for this repository's frozen manuscript or original validation data.

## Source record

These scope statements reflect the October 20, 2026 internal research records:

- `baker_directional_curvature_final_audit_v1_6.md`
- `baker_torus_analytic_lemmas_v1_5.pdf`
- `baker_curvature_verification_closeout_2026_10_20.md`
- `Baker_Directional_Curvature_Referee_Brief_2026_10_20.pdf`

Those files have **not** been added to this repository by this informational update; consult the corresponding directional-curvature verification archive when it is released.

## Research status

This repository remains a frozen, reproducible **v1.0.0** research release with a subsequent documentation-only cross-reference. No new theorem, certified weight-six result, or independent rederivation is claimed.
