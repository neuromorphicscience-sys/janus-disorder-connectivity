# Reproducibility audit

Audit date: 2026-09-28

## Environment

- Platform: Linux 6.6.87.2-microsoft-standard-WSL2, x86_64
- Python: 3.10.12
- NumPy: 2.2.6
- SciPy: 1.15.3
- pandas: 2.3.3
- Matplotlib: 3.10.9
- Numba: 0.65.1
- Pillow: 12.2.0
- pytest: 9.0.3

## Installation and tests

- Isolated virtual-environment package install: PASS
- Archive and source-data validation: PASS (22 hashed source files)
- Lightweight test suite: PASS (15 tests)
- Clean committed-snapshot clone install: PASS
- Windows Git clean clone with `core.autocrlf=true`: PASS under the repository LF policy

## Figure reproduction

All four scripts completed successfully from the supplied compact source tables:

- `scripts/reproduce_fig1.py` -> `figures/generated/fig1.pdf`, `fig1.png`
- `scripts/reproduce_fig2.py` -> `figures/generated/fig2.pdf`, `fig2.png`
- `scripts/reproduce_fig3.py` -> `figures/generated/fig3.pdf`, `fig3.png`
- `scripts/reproduce_fig4.py` -> `figures/generated/fig4.pdf`, `fig4.png`

The PDFs use embedded TrueType fonts. Generated files are intentionally ignored by Git; manuscript reference PDFs are retained in `figures/reference`.

## Integrity checks

- Source-data hashes: PASS against `data/manifests/source_data_sha256.json`
- Absolute Windows path audit: PASS
- Credential and secret-pattern audit: PASS
- Internal-review term audit: PASS
- Missing public input audit: PASS

## Repository size

- Largest tracked file: `data/figure1/schematic_source.png` (1,212,684 bytes)
- Total tracked content: 3401915 bytes

Large production graph arrays, transient caches, internal reviews, and exploratory material are outside this compact public release.
