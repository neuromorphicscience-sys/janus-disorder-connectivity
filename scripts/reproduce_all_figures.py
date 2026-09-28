"""Reproduce all four manuscript figures."""
from __future__ import annotations
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
for number in range(1, 5):
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / f"reproduce_fig{number}.py")],
        cwd=ROOT,
        check=True,
    )
