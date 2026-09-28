"""Reproduce the local-alignment and global-drift predictions."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from janus_connectivity.boundary_theory import fixed_n_bulk_predictions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("figures/generated/local_alignment_predictions.csv"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        fixed_n_bulk_predictions(393216, 128.0, 7, 2.0)
    ).to_csv(args.output, index=False)
    print(args.output)


if __name__ == "__main__":
    main()
