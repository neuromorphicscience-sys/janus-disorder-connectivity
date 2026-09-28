"""Reproduce the clean boundary-source depth profile and total-count law."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from janus_connectivity.boundary_theory import (
    clean_boundary_depth_profile,
    clean_feedback_high_density,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("figures/generated/clean_boundary_profile.csv"))
    args = parser.parse_args()
    depths = np.linspace(0.0, 0.22, 161)
    values = clean_boundary_depth_profile(depths, 393216, 128.0, 7, 2.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {"d_over_R": depths, "mean_Z_up_fixedN_x_averaged": values}
    ).to_csv(args.output, index=False)
    print(f"high-density total: {clean_feedback_high_density(128.0, 7, 2.0):.12f}")
    print(args.output)


if __name__ == "__main__":
    main()
