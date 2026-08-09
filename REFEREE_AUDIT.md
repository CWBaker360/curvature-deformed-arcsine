# Referee-Style Internal Audit

## Verdict

No blocking mathematical inconsistency was found in the current manuscript. The paper is suitable to freeze as a first submission candidate **provided the weight-six status boundary remains explicit**: analytic through geometric weight five; high-precision rational reconstruction at weight six; exact symbolic consequences conditional on that reconstructed side law.

## Substantive revisions made in the audited draft

1. **Curvature convention made explicit.** The geometric four-slot convention and the opposite-sign Brewin source array are stated before the side-law calculation, preventing the sign ambiguity already known to matter in the fifth-order derivation.

2. **Lower-order derivation made standalone.** The manuscript now records the analytic normalized side-law inputs through weight five and explains the perturbative inversion that produces \(\Psi_2,\ldots,\Psi_5\), rather than merely importing those formulas from the earlier curvature paper.

3. **Weight-six status language tightened.** The manuscript uses “high-precision reconstruction” / “computational validation,” not language implying an interval-arithmetic or analytic certificate.

4. **Candidate-basis scope clarified.** The eight weight-six invariants are described as the candidate family used by this side law and recovery design, not as a new abstract completeness theorem for natural invariants.

5. **Exact spherical weight-six checksum included.** The full \(K^3\rho^6\) Laurent kernel agrees identically with the coefficient obtained from the exact spherical deformed-arcsine formula.

6. **Falsification chain retained.** The failed restricted new-angle ansatz is documented because it localized the missing structure to \(K_w^2\); the normalization correction then predicts independent angles at the numerical extraction scale.

7. **Standalone reproducibility repaired.** Scripts that consume retained CSV data now locate the bundle's `data/` directory correctly when run from the package root.

8. **MSC corrected.** The numerical-analysis classification is 65B05 (extrapolation to the limit / deferred corrections), replacing 65D30.

## Remaining proof obligation

The only major unresolved analytic step is an independent degree-eight two-point geodesic-distance / normal-coordinate derivation of the reconstructed weight-six normalized side law. Completing that derivation would promote the weight-six input proposition from high-precision reconstruction to an analytic theorem.

## Scope / novelty boundary

The manuscript does not claim novelty for Schröder/Koenigs linearization, Riemann normal-coordinate expansions, or natural-tensor invariant organization. Its specific claims concern the normalized intrinsic side-law realization, the exact spherical deformed-arcsine model, the explicit variable-curvature kernels, the curved cubic coefficient-lifting mechanism, and the observed finite Laurent-kernel structure for this refinement law.
