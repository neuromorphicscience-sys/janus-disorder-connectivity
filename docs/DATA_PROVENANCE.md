# Data provenance

The public source tables were extracted from the frozen ensembles used by the manuscript. They retain condition labels, sample counts, paired identifiers, available random seeds, graph hashes, and the summary statistics plotted in Figs. 1–4.

The release preparation performed column selection, path normalization, and checksum generation only. It generated no new scientific realizations and did not change any plotted value or selection rule.

The repository excludes full production sparse graphs, per-vertex component masks, temporary NPZ collections, benchmark logs, manuscript-development history, exploratory feature searches, and internal review material. These omissions keep the repository small while preserving the exact plotted data and the algorithms needed to inspect every reported definition.

data/manifests/source_data_sha256.json records hashes of all included source-data files. frozen_identifiers.csv records available seeds and graph hashes. figure_map.json maps each manuscript figure to its input files, script, and output.

## Figure 1(a) vector schematic

The current schematic is drawn from explicit coordinates and directed edges in
`data/figure1/schematic_vector_definition.json` by `scripts/_fig1_schematic.py`.
It is a conceptual illustration, rather than a sampled spatial-network realization.
The local three-edge selection follows the drawn projection ranks; the aligned
branches have strictly downstream links; all eleven highlighted restored vertices
are mutually reachable, with maximum outdegree three.

The supplied reference layout was reconstructed as vector geometry. The earlier
`schematic_source.png` is retained as a historical reference; current figure
reproduction does not read or embed it. PDF and SVG exports contain vector paths
and editable text throughout, including the schematic.
