"""Validate included source data and public imports."""
from __future__ import annotations
from pathlib import Path
import hashlib
import importlib
import json
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


manifest = json.loads((ROOT / "data/manifests/source_data_sha256.json").read_text())
for relative, expected in manifest.items():
    path = ROOT / relative
    if not path.is_file():
        raise SystemExit(f"missing source data: {relative}")
    observed = sha256(path)
    if observed != expected:
        raise SystemExit(f"hash mismatch: {relative}")

required = {
    "data/figure1/restoration_size_summary.csv": {"L", "W", "n", "median__S"},
    "data/figure2/paired_reassignment.csv": {"pair_id", "original_S", "reassigned_S"},
    "data/figure3/internal_recurrence_summary.csv": {"W", "k", "quantity", "median"},
    "data/figure4/boundary_deletion_checks.csv": {"W", "cut_kind", "G_fraction_remaining", "G2_fraction_remaining"},
}
for relative, columns in required.items():
    observed = set(pd.read_csv(ROOT / relative, nrows=2).columns)
    if not columns.issubset(observed):
        raise SystemExit(f"schema mismatch: {relative}")

if len(pd.read_csv(ROOT / "data/figure2/paired_reassignment.csv")) != 48:
    raise SystemExit("expected 48 paired reassignments")
if len(pd.read_csv(ROOT / "data/figure4/clean_feedback_counts.csv")) != 24:
    raise SystemExit("expected 24 clean feedback-count realizations")
for number in range(1, 5):
    if not (ROOT / f"figures/reference/fig{number}.pdf").is_file():
        raise SystemExit(f"missing reference Figure {number}")

for name in [
    "janus_connectivity.model",
    "janus_connectivity.graph_build",
    "janus_connectivity.scc",
    "janus_connectivity.feedback_order",
    "janus_connectivity.recurrence",
    "janus_connectivity.intervention",
    "janus_connectivity.boundary_theory",
]:
    importlib.import_module(name)

files = [p for p in ROOT.rglob("*") if p.is_file() and ".git" not in p.parts]
largest = max(files, key=lambda p: p.stat().st_size)
report = {
    "status": "PASS",
    "source_files_checked": len(manifest),
    "paired_interventions": 48,
    "clean_feedback_realizations": 24,
    "reference_figures": 4,
    "largest_file": str(largest.relative_to(ROOT)),
    "largest_file_bytes": largest.stat().st_size,
    "repository_bytes": sum(p.stat().st_size for p in files),
}
print(json.dumps(report, indent=2))
