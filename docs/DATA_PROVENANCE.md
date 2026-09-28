# Data provenance

The public source tables were extracted from the frozen ensembles used by the manuscript. They retain condition labels, sample counts, paired identifiers, available random seeds, graph hashes, and the summary statistics plotted in Figs. 1–4.

The release preparation performed column selection, path normalization, and checksum generation only. It generated no new scientific realizations and did not change any plotted value or selection rule.

The repository excludes full production sparse graphs, per-vertex component masks, temporary NPZ collections, benchmark logs, manuscript-development history, exploratory feature searches, and internal review material. These omissions keep the repository small while preserving the exact plotted data and the algorithms needed to inspect every reported definition.

data/manifests/source_data_sha256.json records hashes of all included source-data files. frozen_identifiers.csv records available seeds and graph hashes. figure_map.json maps each manuscript figure to its input files, script, and output.
