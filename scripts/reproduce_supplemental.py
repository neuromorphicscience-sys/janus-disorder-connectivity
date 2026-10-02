"""Reproduce the current Supplemental Figs. S1--S9 from archived source tables."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
subprocess.run([sys.executable,str(ROOT/'scripts/reproduce_final_si.py'),
               '--data',str(ROOT/'data/supplemental_final'),
               '--output',str(ROOT/'figures/generated/supplemental')],check=True)
