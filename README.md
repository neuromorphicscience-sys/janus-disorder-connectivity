# Disorder-Restored Strong Connectivity in Degree-Limited Spatial Networks

This is the public reproducibility repository for the manuscript “Disorder-Restored Strong Connectivity in Degree-Limited Spatial Networks.” It contains compact source data, frozen manifests, analysis code, and figure-reproduction scripts supporting the results reported in the manuscript. The model is a fixed-N spatial directed graph in which every source retains at most q neighbors ranked by projection on a quenched local direction.

## Submission version — 2026-10-10

For review of the current manuscript, use the fixed reproducibility snapshot
[`v1.2.0-submission-20261010`](https://github.com/neuromorphicscience-sys/janus-disorder-connectivity/tree/v1.2.0-submission-20261010), or its
[ZIP archive](https://github.com/neuromorphicscience-sys/janus-disorder-connectivity/archive/refs/tags/v1.2.0-submission-20261010.zip).
This snapshot contains the data, code, classifier proof and validation cases, and
all 14 reference figure PDFs corresponding to the submission. Do not use a
moving branch or `releases/latest` to identify the submission data.
[Version and integrity instructions](docs/SUBMISSION_VERSION.md) describe how
to verify the download. The older `v1.0-submission` is a historical snapshot.

## Repository scope

The repository provides:

- fixed-N point and quenched-orientation generation;
- projection-ranked top-q graph construction;
- strongly connected component observables;
- the exact matrix-free feedback-order hierarchy through order two;
- fixed-multiset spatial reassignment and fidelity observables;
- surface and internal recurrence observables;
- fixed-N local-alignment and clean-boundary predictions;
- source tables and scripts for manuscript Figs. 1–4 and Supplemental Figs. S1–S10;
- frozen correlation-length intervention and boundary-deletion replays;
- parameter controls and descriptive finite-size restoration diagnostics.

## Main physical result

Quenched orientational disorder restores a giant strongly connected component while individual sources remain strongly direction selective. Reassigning the same orientation multiset over a fixed point cloud suppresses strong connectivity by orders of magnitude, identifying spatial assignment as causally important for global reachability. The intervention also modifies joint structural statistics, including indegree heterogeneity and edge reciprocity, and therefore does not isolate a unique mediator of connectivity. Near restoration, the largest low-order recurrent components form a direction-selected surface sector, while substantial low-order recurrence remains fragmented in the interior. Macroscopic mutual reachability therefore depends on the collective organization of recurrence, rather than its mere presence.

## Repository structure

- src/janus_connectivity: model, graph, SCC, feedback-order, intervention, recurrence, and theory code.
- scripts: validation and Fig. 1–4 reproduction entry points.
- data: compact source tables and frozen identity manifests.
- theory: standalone fixed-N and clean-boundary calculations.
- tests: lightweight exact and invariant tests.
- figures/reference: manuscript figure PDFs used as visual references.
- docs: model definitions, provenance, and reproduction details.

## Installation

Python 3.10 or newer is recommended.

    python -m venv .venv
    source .venv/bin/activate
    python -m pip install -r requirements.txt
    python -m pip install .

On Windows, activate the environment with .venv\Scripts\activate.

## Quick start

    python scripts/validate_archive.py
    python -m pytest
    python scripts/reproduce_fig1.py
    python scripts/reproduce_fig2.py
    python scripts/reproduce_fig3.py
    python scripts/reproduce_fig4.py

To reproduce all figures:

    python scripts/reproduce_all_figures.py

Generated PDF and PNG files are written to figures/generated.

## Reproducing the manuscript figures

| Figure | Main source data | Script | Expected PDF |
|---|---|---|---|
| Fig. 1 | data/figure1 | scripts/reproduce_fig1.py | figures/generated/fig1.pdf |
| Fig. 2 | data/figure2 | scripts/reproduce_fig2.py | figures/generated/fig2.pdf |
| Fig. 3 | data/figure3 | scripts/reproduce_fig3.py | figures/generated/fig3.pdf |
| Fig. 4 | data/figure4 | scripts/reproduce_fig4.py | figures/generated/fig4.pdf |

The reference PDFs in figures/reference are the manuscript versions. Small rendering differences can arise from local font and Matplotlib versions; the plotted values and definitions are fixed by the supplied source tables.

## Reproducing the Supplemental Material

    python scripts/analyze_final_si.py --data data/supplemental_final --output outputs/supplemental_statistics
    python scripts/validate_si_classifier.py --data data/supplemental_final --output outputs/supplemental_classifier
    python scripts/validate_full_graph_classifier.py
    python scripts/reproduce_supplemental.py

These commands recompute descriptive statistics, independently validate the five archived
toy graphs, compare every feedback edge on exhaustive four-vertex topologies and complete
small-to-medium model graphs with two unpruned oracles, and reproduce Figs. S1–S10 as vector PDF/SVG and 600 dpi PNG.
Exact sample counts, rejected intervention attempts, source-table mappings, and optional
frozen-graph replay are documented in [Supplemental reproduction](docs/SUPPLEMENTAL_REPRODUCTION.md).
The single-column source is [Supplemental_Material.tex](supplemental/Supplemental_Material.tex).
No new base-network realization is needed to reproduce these figures. The separate
full-graph algorithm check generates only its fixed validation cases; it does not change
the frozen physical ensembles. See [classifier proof and validation](docs/CLASSIFIER_CORRECTNESS.md).

## Frozen data and provenance

The plotted ensembles and interventions were frozen before repository preparation. data/manifests/frozen_identifiers.csv records available seeds and graph hashes; ensemble_counts.csv records plotted sample counts; source_data_sha256.json verifies every supplied data file. The source tables are derived from the archived production ensembles and were not selected by visual appearance.

Large production graph arrays and transient caches are omitted. The compact source tables reproduce the reported figures, while the public algorithms and tests expose the exact graph definitions and feedback-order classifier. See docs/DATA_PROVENANCE.md.

## Computational notes

Figure reproduction from supplied source data is lightweight. A full production replay at the largest manuscript sizes requires substantially more memory and CPU time. The order-two classifier is sparse and matrix free: it does not materialize a feedback meta-graph, an M by M reachability matrix, or a transitive-closure matrix.

## Citation

The manuscript has not been assigned a publication DOI. Citation metadata for this submission snapshot is provided in CITATION.cff.

## License

Source code is licensed under the BSD 3-Clause License. Data and documentation are licensed under CC BY 4.0 unless otherwise noted. See `LICENSE` and `LICENSE-DATA`.

Third-party components, including embedded fonts, remain subject to their original licenses and are not relicensed by this repository.

## Final matched bridge and submission diagnostics

The final source update adds the 48-pair intervention–filtration bridge, classifier coverage audit, shared-stream crossover intervals, clean depth profile and complete Main Fig.4(c) deletion ensemble. See [definitions and source data](data/final_bridge/README.md) and [reproduction notes](docs/FINAL_BRIDGE_REPRODUCTION.md). The reassigned strict cores have no nontrivial low-order recurrence, making their conditional concentration undefined; the intervention and natural-disorder filtration are complementary evidence.
