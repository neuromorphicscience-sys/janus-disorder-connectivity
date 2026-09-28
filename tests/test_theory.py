import json
from pathlib import Path
import numpy as np
import pandas as pd
from janus_connectivity.boundary_theory import (
    fixed_n_bulk_predictions,
    clean_boundary_depth_profile,
    clean_feedback_high_density,
)

ROOT = Path(__file__).resolve().parents[1]


def test_fixed_n_alignment_matches_archived_prediction():
    expected = pd.read_csv(ROOT / "data/figure1/bulk_fixedN_predictions.csv")
    actual = fixed_n_bulk_predictions(393216, 128.0, 7, 2.0)
    assert abs(actual[0]["local_alignment_fixedN"] - expected.local_alignment_fixedN.iloc[0]) < 2e-13


def test_clean_boundary_profile_matches_source_table():
    expected = pd.read_csv(ROOT / "data/figure4/boundary_depth_prediction.csv")
    select = np.array([0, 40, 80, 120, 160])
    depth = expected.d_over_R.to_numpy()[select]
    actual = clean_boundary_depth_profile(depth, 393216, 128.0, 7, 2.0)
    assert np.allclose(actual, expected.mean_Z_up_fixedN_x_averaged.to_numpy()[select], rtol=0, atol=2e-12)


def test_high_density_boundary_law():
    expected = json.loads((ROOT / "data/figure4/clean_feedback_prediction.json").read_text())
    assert abs(clean_feedback_high_density(128.0, 7, 2.0) - expected["high_density_square_prediction"]) < 1e-12
