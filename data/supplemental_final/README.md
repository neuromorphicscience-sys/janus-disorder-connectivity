# Supplemental source tables

Figures S1--S9 read the graph-level CSV tables in this directory.
The current supplement also includes Figure S10, whose source tables are in
`data/final_bridge`; `scripts/reproduce_supplemental.py` reproduces S1--S10.
Green/circle and blue/square curves denote G1 and G2; full-graph results use magenta.
All ensemble summaries use medians and IQRs. Counts describe overlapping cohorts.

| Figure | Source files |
|---|---|
| S1 | source_figS1.csv |
| S2 | source_figS2.csv |
| S3 | source_figS3.csv; finite_size_archived_summary.csv |
| S4 | source_figS4.csv |
| S5 | source_figS5.csv; ell_base_manifest.csv |
| S6 | source_figS6.csv; classifier_toy_nodes.csv; classifier_toy_edges.csv |
| S7 | source_figS7.csv; source_figS7_scaling.csv; source_figS7_direction.csv; direction_base_manifest.csv |
| S8 | source_figS8.csv |
| S9 | source_figS9.csv; deletion_base_manifest.csv |
| Table S1 | source_tableS1.csv; source_tableS1_print.csv |
| Table S2 | source_tableS2.csv; source_figS5.csv |
| Table S3 | source_tableS3.csv; source_figS8.csv |

source_provenance.json maps extracted tables to original archived content hashes.
Seeds are decimal identity strings, not floating-point values. Missing rejected surrogate
SCC outcomes are intentional: those attempts failed fidelity before connectivity evaluation.
Covariance lengths in family controls are the actual smoothing-kernel lengths.
The finite-size cohort has 4--150 graphs per sampled (L,W); the internal cohort has 8--96.

Plotting: scripts/reproduce_final_si.py. Statistics: scripts/analyze_final_si.py.
Independent unchanged-toy validation: scripts/validate_si_classifier.py.
