# Source data

The figure directories contain compact tables used directly by the reproduction scripts.

- figure1: restoration curves, frozen bulk marginal checks, fixed-N predictions, and the schematic source image.
- figure2: 48 paired fixed-multiset interventions and the displacement and edge-length controls.
- figure3: feedback-order scans, finite-size summaries, depth and direction interventions, occupancy profiles, and internal recurrence summaries.
- figure4: clean feedback counts, boundary depth theory, deletion checks, and the parameter-free prediction record.
- manifests: available frozen seeds and hashes, ensemble counts, figure-to-input mapping, and SHA-256 checksums.

These are source tables derived from frozen production ensembles. Large per-realization sparse graphs, component masks, and temporary production caches are not included in ordinary Git history. The supplied tables are sufficient to reproduce Figs. 1–4. Model realizations can be regenerated from recorded seeds where both seeds are available; graph hashes identify frozen realizations used for independent checks.

## Supplemental source data

`data/supplemental` supplies the 48-base paired correlation-length scan, the complete 192-base deletion ensemble at three widths, the archived parameter controls, the combined finite-size diagnostics, the 784-member low-order mass cohort, and the 1,528-member internal-recurrence cohort. Graph hashes distinguish realizations; manifests provide exact seeds for the new reanalyses. See [Supplemental reproduction](../docs/SUPPLEMENTAL_REPRODUCTION.md) for populations, definitions, and scripts.
