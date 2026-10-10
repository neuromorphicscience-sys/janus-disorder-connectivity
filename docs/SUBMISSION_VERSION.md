# Submission reproducibility snapshot — 2026-10-10

Version: `v1.2.0-submission-20261010`. This version corresponds to the manuscript **Disorder-Restored Strong Connectivity in Degree-Limited Spatial Networks**, with the clean-limit benchmark wording, intervention mediator qualification, 11-page Supplemental Material, and classifier-bound proof/full-graph validation.

Use [this tagged snapshot](https://github.com/neuromorphicscience-sys/janus-disorder-connectivity/tree/v1.2.0-submission-20261010) or [download its ZIP](https://github.com/neuromorphicscience-sys/janus-disorder-connectivity/archive/refs/tags/v1.2.0-submission-20261010.zip). The submission uses the data and algorithms frozen here. Subsequent development on `main` is a separate version. The old `v1.0-submission` tag predates the present submission; retain it only for historical provenance. Do not move or reuse submission tags. Corrections require a new version.

## Scope and integrity

The snapshot supplies compact source tables for main Figs. 1–4 and supplemental Figs. S1–S10, 14 exact submission reference PDFs, deterministic model/analysis scripts, the classifier proof and independent full-graph validation, archived labels, sample inventories, and licenses. Large transient production arrays and build caches are excluded. Data and algorithms are unchanged from the previously validated classifier revision; the four main reference PDFs are synchronized to the submission's vector figures with transparent legends.

`data/manifests/source_data_sha256.json` checks the supplied data, selected analysis sources, and all 14 reference figures. `data/manifests/submission_payload_sha256.json` additionally freezes all model code, scripts, tests, theory, dependency declarations and license files. Its keys are repository-relative paths, and values are SHA-256 digests of the distributed file bytes. The manifest excludes itself and version-navigation documents so it does not contain a circular self-hash.

After extracting the archive, run:

```sh
python scripts/validate_archive.py
python scripts/validate_full_graph_classifier.py
python -m pytest
```

The full-graph command generates only the fixed algorithm-validation cases. Figure reproduction uses the supplied source tables and does not require replaying the large physical ensembles. Generated outputs belong under `outputs/` and `figures/generated/` and do not alter the submission snapshot.

The manuscript's data links identify the snapshot by its full commit SHA; the tag is a readable version label. Use the SHA-addressed archive if a tag and a commit ever disagree. The authors' final TeX/PDF files and their checksums are supplied separately from this data/code snapshot.
