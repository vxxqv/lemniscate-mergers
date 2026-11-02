# Adversarial self-review, 2 November 2026

Primary verdict: **proved as written**, for the four stated analytic theorems
and their stated corollaries, after the restrictions recorded below. This is a
self-review by the same assistant, not an independent mathematical audit or a
formal proof certificate. Exact reviewed file hashes are in `release_hashes.json`.

## Dependency structure

1. Factorized cluster normalization -> uniformly convergent log series ->
   moment and derivative tail estimates -> modulus bracket.
2. Modulus bracket -> sublevel inclusions and minimax comparison -> equality
   of labelled partitions when the two bounding partitions coincide.
   Boundary exclusion is a separate required leaf for global cluster claims.
3. Complex Taylor remainder -> strict Rouche comparison of F' with its affine
   Taylor polynomial -> unique simple critical point -> completing the complex
   square -> explicit stationary-value disk and modulus interval.
4. Exact quadratic identity -> minimum separation below the saddle and a
   uniform tube above it -> lattice covering radius and Voronoi adjacency ->
   root-to-root graph connectivity for all rotations and translations.

These arguments use standard Rouche, complex Taylor, implicit-function, and
elementary compactness/Voronoi facts. No main theorem from an inaccessible
specialist paper is used as an unchecked premise.

## Obligation matrix

| Obligation | Verdict | Check |
| --- | --- | --- |
| Logarithm branch and convergence | Passed | Each epsilon R/|d_l| < 1; normalized log is zero at w=0. |
| Signs of inverse-moment coefficients | Passed | log(1+x) and its derivative reconstructed term by term. |
| K=0 and p=1 | Passed | Empty finite sum, geometric tail still valid. |
| Constant external factor | Passed after explicit convention | N>=1 in main setting; G=1 handled separately. |
| Exact versus nearly vanishing moments | Passed | Order corollary requires exact cancellation; broken-ring experiment shows the first-order return. |
| Restricted versus full-plane paths | Passed | Explicit strict boundary barrier, not merely a local norm estimate. |
| Event ordering versus root labels | Passed | Partition sandwich and reflection example retain the distinction. |
| Critical-point uniqueness | Passed | Strict Rouche inequality counts derivative zeros with multiplicity one. |
| Complex quadratic correction | Passed | f1 squared, not |f1| squared; completion checked with translated nonreal quadratic. |
| Third-derivative bound | Passed | Finite Taylor coefficient majorant over entire closed complex disk. |
| Frozen-value comparison | Passed | Error O(epsilon^(2p)) is against F_e(gamma), not Q(gamma). |
| Digital-graph theorem scope | Passed after explicit restriction | Exact vertex classification, root pair only, eta<=1/16, h/W<1/sqrt(2), sufficient finite window. |
| Voronoi vertex crossing | Passed | Insert intermediate nearest vertex to retain 4-neighbor path. |
| Lower-level endpoint inclusion | Passed | 2sqrt(eta)+eta <= 1-eta on stated eta range. |
| Numerical finite coverage | Passed | 64 critical centers, 32 brackets, 5376 raster cases, no exclusions. |
| Floating raster decisions as certification | Out of scope | Not claimed; independent theorem supplies exact geometric statement. |
| Universal novelty | Not addressed | Full-text overlap and expert review still pending. |

## Decisive recomputations

- The stationary value of F(z)=3+(z-2i)^2 at arbitrary center z0=1 is exactly
  3, with M3=0. The ball implementation reproduces it exactly.
- A cubic with a multiple derivative zero at the selected center is rejected
  when f2=0; no simple-saddle assertion is inferred.
- Factorized and expanded evaluations agree through 85 decimal places in
  eight cross-checks using a separate mpmath evaluation path.
- All 2304 fine-grid trials protected by the exact theorem were checked
  individually. This consistency check is finite evidence, not its proof.
- The upper tube identity was evaluated at 63 rational-coordinate samples.
  Its proof uses the exact polynomial identity for the whole stated interval.

## Corrections during review

- Removed an extra brace in one displayed formula before compilation.
- Stated N>=1 explicitly and separated the constant-factor case.
- Retained the aligned-field epsilon=1/2, gamma=-1 enclosure failure.
- Replaced the first raster window implementation, whose coarse translated
  origin could crop a root, with outward-rounded indices covering a fixed
  window. Run main-002 is superseded; main-003 is authoritative.
- Clarified that the predictor is the stationary value of the quadratic
  Taylor polynomial, not evaluation after a Newton step.

## Remaining gaps to a publishable claim

Complete the exact overlap comparison with Farouki-Han, the polynomial-fireworks
chapter, and validated Newton/Taylor models. Obtain independent expert review.
A fair adaptive-method performance benchmark has not been run; the draft makes
no superiority claim. Extending the quadratic lattice theorem requires uniform
control of Morse-coordinate distortion and a global connecting path.
