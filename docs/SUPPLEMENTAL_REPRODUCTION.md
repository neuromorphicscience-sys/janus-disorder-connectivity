# Supplemental reproduction

The supplemental tables use archived graph ensembles and deterministic reanalysis. No new base-network realizations are included.

## Source tables and figures

    python scripts/validate_archive.py
    python scripts/summarize_supplemental.py
    python scripts/reproduce_supplemental.py
    python -m pytest

The summary script recomputes condition medians, IQRs, variances, bootstrap median intervals, bracketed crossover diagnostics, intervention summaries, and deletion summaries from the included graph-level tables. It compares these with the archived source summaries. Outputs go to `outputs/supplemental_summaries`. The figure script writes PDF, SVG, and PNG files to `figures/generated/supplemental`.

| Figure | Evidence | Source subdirectory under data/supplemental |
|---|---|---|
| S1 | Connection budgets q=4,5,6,7,8,10,12 | 04_parameter_robustness |
| S2 | Native distribution parameters and spatial correlations | 04_parameter_robustness |
| S3 | All 48 canonical intervention pairs | 02_ell_sensitivity |
| S4 | Paired correlation lengths 1,2,4,8 R | 02_ell_sensitivity |
| S5 | Finite-size restoration medians and interpolated crossings | 05_finite_size_diagnostics |
| S6 | Low-order largest-component mass and geometry | 05_finite_size_diagnostics |
| S7 | Internal participation, concentration, counts and retention | 06_si |
| S8 | All eligible boundary deletions and representative placement | 03_boundary_deletion_ensemble |

The single-column source is `supplemental/supplemental.tex`. It uses the reference PDFs in `figures/reference/supplemental`. Compile it from its own directory with a standard LaTeX installation. The accompanying Letter supplies the literature reference list.

## Frozen graph replay

The optional replay script reconstructs only graph identities in the supplied manifests:

    python scripts/replay_final_hardening.py ell --workers 4
    python scripts/replay_final_hardening.py deletion --workers 4

A selected frozen identity can be replayed separately:

    python scripts/replay_final_hardening.py ell --graph-id S1C_W0.790_r00 --workers 1
    python scripts/replay_final_hardening.py deletion --graph-id gaussian_iid_L128_q7_W0.815_r34 --workers 1

Outputs go to `outputs/final_hardening`. Each base graph must reproduce its CSR hash before analysis. The canonical 4 R intervention must reproduce the historical reassigned graph hash and SCC fraction. Deletion orders are classified in the complete graph, checked against archived full-graph low-order masses, and held fixed when inducing the three cores. No outgoing edge is replaced after deletion. Archived 4 R core values are also checked where available.

The replay field is the canonical 256 by 256 grid with unit spacing, periodic Gaussian filtering of width ell/sqrt(2), offset (64,64), and bilinear interpolation. It is intentionally specified explicitly in this script; a generic field helper with another canvas-to-domain mapping does not reproduce these frozen pairs. The same white-noise seed is reused across ell values for each historical realization. There is one attempt per base and length. Full-edge marginal distances, rate change, outdegree equality, multiset equality, and changed-edge fraction are saved before surrogate SCC calculation. Rejected attempts retain fidelity diagnostics and have no surrogate connectivity outcome in the reported comparison.

## Populations and definitions

The robustness inventory contains 29,148 distinct hashes. Independent disorder families use their own parameters: Gaussian W, uniform half-width A, von Mises concentration kappa, and Laplace scale b. Correlated scans use latent mixing alpha at fixed Gaussian W=0.82. Numeric ell_x and ell_y columns give the actual covariance scales from the generating kernel; historical family names may contain a different nominal length. In particular, the old rank_iso_2R family has actual ell=4 R.

The full-graph finite-size inventory contains 4,634 distinct hashes at four sizes. It includes overlapping archived campaigns after deduplication. The 784-member low-order mass cohort and 1,528-member internal-recurrence cohort retain their own definitions and sample inventories. Counts from overlapping analyses should not be added as independent realizations.

The boundary-deletion population contains every eligible frozen L=128, q=7, sigma=24 Gaussian graph at W=.805,.815,.86,.90: respectively 96,48,24,24. Each has widths 2,4,8 R. Graph-level values retain the complete observed spread. Main-figure representatives were selected from the earlier mass cohort by proximity of full-graph S2 to its condition median, with graph identifier as tie-break; selection did not optimize a deletion outcome. `REPRESENTATIVE_PLACEMENT.csv` gives the earlier selection rule and midrank percentiles in the enlarged deletion ensemble.

`SOURCE_TABLE_IDENTITIES.csv` records the names and content digests of archived source tables. Each compact row has a graph hash and, where applicable, a source identity. The included data manifest verifies the actual public files. Original private directory paths and transient graph caches are omitted.
