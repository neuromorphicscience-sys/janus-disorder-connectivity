"""Evaluate the exact fixed-N bulk selection kernel."""
from __future__ import annotations
from janus_connectivity.boundary_theory import (
    circular_cap_area,
    fixed_n_bulk_predictions,
)

__all__ = ["circular_cap_area", "fixed_n_bulk_predictions"]

if __name__ == "__main__":
    import json
    rows = fixed_n_bulk_predictions(393216, 128.0, 7, 2.0)
    print(json.dumps(rows, indent=2))
