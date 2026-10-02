# Reproducing the Supplemental Material

From the repository root, with the project and requirements installed:

    python scripts/validate_archive.py
    python scripts/analyze_final_si.py --data data/supplemental_final --output outputs/supplemental_statistics
    python scripts/validate_si_classifier.py --data data/supplemental_final --output outputs/supplemental_classifier
    python scripts/reproduce_supplemental.py

The plotting command creates nine vector PDFs, editable SVGs, and 600 dpi PNGs in
figures/generated/supplemental. Every point comes from the CSV files listed in
data/supplemental_final/README.md. The font is DejaVu Sans, bundled with Matplotlib.
Reference exports are supplied in figures/reference/supplemental.

To additionally export legend bounds and data/text intersection diagnostics:

    python scripts/reproduce_final_si.py --data data/supplemental_final --output figures/generated/supplemental --legend-audit

These JSON diagnostics complement visual inspection of the final PDF.

To compile the single-column supplement with an existing LaTeX installation:

    cd supplemental
    pdflatex -interaction=nonstopmode -halt-on-error Supplemental_Material.tex
    pdflatex -interaction=nonstopmode -halt-on-error Supplemental_Material.tex

The source uses standard article, geometry, newtx, graphicx, booktabs, float and hyperref
packages. The scientific supplement is also supplied as a compiled reference PDF.

## Archived analyses

The ell/R grid is 1,2,4,8 on 48 historical bases at W=0.79,0.80,0.81.
Acceptance counts are 48,48,48,42; six 8R failures retain diagnostics and missing SCC
outcomes. All scales use KS/W1 <=0.03, relative rate change <=0.03, changed edges >=0.25,
and preserved outdegrees and orientation multisets. Streams recur across W and ell;
pooled summaries are descriptive. No new base realization is required for plotting.

The deletion archive contains all 192 eligible graphs at the four selected conditions:
n(W)=96,48,24,24 for W=0.805,0.815,0.86,0.90, respectively.
Each has widths 2R,4R,8R. Complete-graph feedback labels are restricted by induction.
Connections are not refilled. Fractions use the retained core node count.

The full internal-recurrence cohort contains 1,528 graphs at 24 W values.
Largest-component geometry uses 784 archived graphs; direction controls use 24
matched bases at three orientations and co-rotating classifier coordinates.
The full finite-size cohort contains 4,634 graphs and only bracketed interpolated crossings.

Optional computational replay of existing frozen identities remains available through
scripts/replay_final_hardening.py and archived manifests in data/supplemental.
It is separate from source-table figure reproduction; see the script's --help.
The production feedback classifier is sparse and matrix free through order two.
The tiny validation graphs are archived toy definitions, not new physical realizations.
