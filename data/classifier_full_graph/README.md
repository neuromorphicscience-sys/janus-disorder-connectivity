# Complete-graph classifier checks

These deterministic algorithm-validation results supplement the five archived
toys and 80 induced-subgraph checks. They are separate from the frozen physical
ensembles plotted in the manuscript.

Run `python scripts/validate_full_graph_classifier.py` from the repository root.
See `docs/CLASSIFIER_CORRECTNESS.md` for the proof, fixed suite and two independent
reference algorithms. Results are written to `outputs/full_graph_classifier`.

- `FIXED_CASES.json`: the 26 complete-model cases, parameters and seeds.
- `FOUR_VERTEX_CASES.csv`: every one of the 4096 directed four-vertex topologies.
- `MODEL_CASES.csv`: complete-model counts, graph hashes and agreement status.
- `MODEL_EDGE_COMPARISONS.csv`: every complete-model directed edge and its labels.
- `ARCHIVED_SUBGRAPH_RECHECK.csv`: all 80 archived subset labels, also checked after a negative coordinate translation.
- `SUMMARY.json`: aggregate counts, checks and production classifier source hash.

Label 0 denotes a downward edge, 1/2 an upward edge of exact order one/two, and
3 an upward edge of higher or infinite return order. The script compares both
references to production; the edge CSV records the state-space reference labels.
