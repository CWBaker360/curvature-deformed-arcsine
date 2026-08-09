# Submission Metadata - Version 1.0.0

## Full title

**A Curvature-Deformed Arcsine Coordinate for Cubic Geometric Refinement: Exact Spherical Conjugacy, Riemannian Expansion Through Fifth Order, and a High-Precision Weight-Six Reconstruction**

## Short title

**A Curvature-Deformed Arcsine Coordinate**

## Author

Wayne Baker  
Independent researcher

## Manuscript date

August 8, 2026

## Version

1.0.0

## Abstract

A normalized geodesic side law on a Riemannian surface defines, by inversion, a local linearizing coordinate for geometric encode--scale--decode refinement. In the Euclidean limit this coordinate is the ordinary arcsine function. This paper develops the resulting curvature-deformed arcsine

\[
\Psi_{\rho,w}(z)=\frac12 H_{\rho,w}^{-1}(2z),
\]

where \(H_{\rho,w}\) is the normalized signed side law for equal-radius geodesic endpoints based at \(p\) and centered about a unit tangent direction \(w\).

The construction gives an exact Schroeder/Koenigs conjugacy

\[
\Psi_{\rho,w}\!\left(\mathcal{T}_{\rho,w,q}(z)\right)
=q\,\Psi_{\rho,w}(z),
\]

and hence an exact local multiplicative semigroup. On a sphere the deformed arcsine is available in closed form,

\[
\Psi_\alpha(z)
=
\arcsin\!\left(
\frac{\sin(z\sin\alpha)}{\sin\alpha}
\right),
\]

which interpolates between \(\arcsin z\) in the flat limit and the identity map at the hemispherical value \(\alpha=\pi/2\).

For a general Riemannian surface the variable-curvature expansion is derived analytically through geometric weight five. The coefficient functions are finite Laurent polynomials in \(r=\sqrt{1-z^2}\), and therefore every geometric-weight row of the Taylor coefficients is a finite combination of generalized binomial sequences. At weight six, where the available geodesic-distance expansion no longer directly supplies the required degree-eight side-law term, a high-precision intrinsic torus computation is used to recover a rational candidate formula in eight natural curvature invariants. Independent geometry and angle holdouts, including a failed restricted ansatz and its subsequent correction, provide a stringent computational validation. Conditional on that recovered side law, exact symbolic inversion yields the complete weight-six Laurent kernels and all Taylor coefficients at that weight.

The results motivate a finite Laurent-kernel conjecture for the full curvature hierarchy and extend the cubic \(N\)-series coefficient-lifting mechanism from the ordinary arcsine to its Riemannian deformation.

## 2020 Mathematics Subject Classification

- 53B20 - Local Riemannian geometry
- 53A55 - Differential invariants (local theory)
- 39B12 - Iteration theory, iterative and composite equations
- 65B05 - Extrapolation to the limit, deferred corrections

## Keywords

Riemannian normal coordinates; geodesic distance; Schroeder equation; Koenigs linearization; Gaussian curvature; asymptotic expansion; geometric refinement.

## Research-status statement

The manuscript is analytic through geometric weight five. The general weight-six normalized side-law coefficient is a high-precision rational reconstruction supported by intrinsic torus computations, unused-geometry holdouts, unused-angle tests, and an exact spherical `K^3 rho^6` checksum. It is not presented as an interval-arithmetic proof. Conditional on the reconstructed side law, the complete weight-six Laurent kernels and all-`n` formulas are exact symbolic consequences.

The outstanding analytic proof obligation is an independent derivation of the degree-eight two-point geodesic-distance term sufficient to prove the reconstructed general weight-six side law directly.

## Data and code availability statement

All scripts and retained numerical data used for the weight-six reconstruction and verification are distributed with the version 1.0.0 reproducibility archive. The archive includes raw geodesic data, intermediate and corrected angular-model outputs, exact symbolic audits, SHA-256 checksums, and a verification record.

## Reproducibility note

The release deliberately retains the failed initial restricted angular ansatz and the intermediate `UNRESOLVED_K_W_SQUARED` labels. These are historical intermediate results, not unresolved parameters in the final manuscript. The new-angle failure localized to the `K_w^2` sector and led to the normalization-forced correction used in the final formula.

## Conflict of interest

The author declares no conflict of interest.

## Funding

No external funding is declared in this release.

## License

Manuscript, documentation, and distributed data: CC BY 4.0.  
Python source code: MIT License.
