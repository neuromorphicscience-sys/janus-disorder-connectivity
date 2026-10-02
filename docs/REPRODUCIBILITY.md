# Reproducibility workflow

## Lightweight validation

    python scripts/validate_archive.py
    python -m pytest

The validator checks source-data hashes, required schemas, reference outputs, and imports. The test suite checks graph construction against a brute-force implementation, SCC invariants, feedback orders one and two, an order-two edge whose partner is order one, graph nesting, internal recurrence, deletion monotonicity, the fixed multiset, and analytical predictions.

## Figures

Run any figure separately or all four:

    python scripts/reproduce_fig1.py
    python scripts/reproduce_fig2.py
    python scripts/reproduce_fig3.py
    python scripts/reproduce_fig4.py
    python scripts/reproduce_all_figures.py

Outputs are written to figures/generated. Reference manuscript PDFs are in figures/reference.

## Full realization replay

The source package exposes deterministic model generation and sparse analysis algorithms. Full replay of the largest ensembles is intentionally separate from figure reproduction because it is computationally expensive. The public release supplies compact derived source tables, ensemble counts, and available seed/hash identities rather than large production graph caches.

## Supplemental statistics and figures

    python scripts/summarize_supplemental.py
    python scripts/reproduce_supplemental.py

See [Supplemental reproduction](SUPPLEMENTAL_REPRODUCTION.md) for optional deterministic graph replay and the actual sample inventories. The canonical field recipe and wrapped-angle fidelity definition are given there explicitly.
