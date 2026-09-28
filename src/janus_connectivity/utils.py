"""Small deterministic utilities."""
from __future__ import annotations
from pathlib import Path
import hashlib
import json


def sha256(path) -> str:
    path = Path(path)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
