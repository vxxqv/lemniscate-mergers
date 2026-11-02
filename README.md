# Quantitative lemniscate mergers near clustered zeros

Reproducible code and data, version 0.1.0.

**Author:** Vivaan Patni, Independent Researcher, Dubai, UAE  
**ORCID:** https://orcid.org/0009-0005-1859-5107  
**Contact:** contact.vpatni@gmail.com

## Scientific scope

This package evaluates retained inverse-moment bounds for the effect of distant
zeros on a small polynomial-root cluster, local critical-value enclosures from a
quadratic Taylor model, and grid resolution near an exact quadratic merger.
The associated manuscript is a working draft. Proofs have received assistant
self-review only; independent review and comprehensive historical originality
verification remain pending. The package does not claim an advantage over a
faithful implementation of an established adaptive contouring method.

## Reproduce

Use Python 3.12 and install the pinned dependencies:

```bash
python3 -m pip install -r requirements.txt
python3 src/experiments.py --out computations/reproduction-001
python3 tests/check_mathematics.py
python3 src/figures.py
```

The retained final scientific run is `computations/main-003`. Use a new output
directory when reproducing it. Figures and the consistency checks deliberately
read that retained run; change the `DATA` constant in `src/figures.py` explicitly
to inspect a new run. The scientific experiment uses deterministic inputs.

## Recorded results

- 64 critical-center cases, with 63 successful strict localization/value bounds.
- 32 hierarchy cases: 26 moment-aware certificates versus 12 distance-only ones.
- 5376 raster trials: 4069 correct, 619 wrong and 688 roots unresolved.
- All 2304 trials with h/W <= 1/2 are correct. A separate exact-model proof gives
  a sufficient condition valid for all lattice translations and orientations.
- The unresolved local enclosure is the aligned field at epsilon=1/2, gamma=-1.

The exact inputs enter the Arb/acb path as integers or rational strings at
384-bit precision. A separate 95-digit mpmath solve provides numerical reference
values. The raster diagnostic uses float64; it is not a pixel certificate.

## Data and provenance

`computations/main-003/` contains CSV results, full ball inequality records,
software versions and a hashed execution manifest. Earlier runs are preserved:
`main-001` failed at an output serialization interface; `main-002` is superseded
because its coarse translated raster window could crop a root. All conclusions
and figures use `main-003`.

`src/` contains the analytic implementation, deterministic experiment and figure
code. `tests/` contains orthogonal checks. `figures/` contains editable vector
exports and PNG previews. `records/` documents the target, proof self-review,
source checks and limitations. No reuse license has been selected by the author. The manuscript PDF is
distributed separately.

## Citation and archive status

Version 0.1.0 is archived at Zenodo: https://doi.org/10.5281/zenodo.23889014.
The matching GitHub release is https://github.com/vxxqv/lemniscate-mergers/releases/tag/v0.1.0.
Creator metadata is supplied in `CITATION.cff` and `.zenodo.json`. The archived
ZIP preserves the files as deposited; this README and citation metadata on the
main branch were updated after the DOI was registered.

## AI assistance and review

OpenAI ChatGPT/Codex substantially assisted derivations, code, experiments,
figures, source checking and drafting. Generated plots use explicit formulas
and recorded numerical results. No independent mathematical reviewer is claimed.
