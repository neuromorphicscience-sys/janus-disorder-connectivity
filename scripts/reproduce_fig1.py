"""Reproduce manuscript Figure 1 from supplied source data."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
from _figure_common import DATA, OUTPUT, COLORS as C, apply_style, panel, save
from _fig1_schematic import draw_schematic

apply_style()
plt.rcParams.update({"legend.fontsize": 7, "lines.linewidth": 1.15, "lines.markersize": 3.5,
                    "axes.linewidth": .7, "xtick.major.size": 3, "ytick.major.size": 3,
                    "xtick.minor.size": 1.7, "ytick.minor.size": 1.7,
                    "xtick.major.width": .65, "ytick.major.width": .65,
                    "xtick.minor.width": .5, "ytick.minor.width": .5,
                    "xtick.major.pad": 3, "ytick.major.pad": 3,
                    "axes.labelpad": 2, "legend.frameon": False, "legend.handlelength": 1.6})
styles = [(64, "o", "-"), (128, "s", "--"), (256, "^", "-."), (512, "D", ":")]

definition = json.loads((DATA / "figure1" / "schematic_vector_definition.json").read_text())
ref_w, ref_h = definition["canvas"]
fig = plt.figure(figsize=(86 / 25.4, 90 / 25.4))
a = fig.add_axes([1 / 86, 52 / 90, 84 / 86, (84 * ref_h / ref_w) / 90])
b = fig.add_axes([0.155, 9.96 / 90, 0.355, 33.615 / 90])
c = fig.add_axes([0.665, 9.96 / 90, 0.315, 33.615 / 90])
geometry = draw_schematic(a, definition, C)

panel(b, "b")
panel(c, "c")
size = pd.read_csv(DATA / "figure1" / "restoration_size_summary.csv")
scan = pd.read_csv(DATA / "figure1" / "restoration_W_scan.csv").sort_values("W")
for L, marker, linestyle in styles:
    z = size[size.L.eq(L)].sort_values("W")
    x, y = z.W.to_numpy(), z.median__S.to_numpy()
    if L == 128:
        low = scan[(scan.W >= 0.75) & (scan.W < z.W.min())]
        x = np.r_[low.W.to_numpy(), x]
        y = np.r_[low.S_median.to_numpy(), y]
    b.plot(x, y, color=C["full"], marker=marker, mfc="white", ms=3.3, ls=linestyle, label=str(L))
b.set(xlabel=r"$W$", ylabel=r"$S$", xlim=(0.748, 0.904), ylim=(-0.025, 1.045), xticks=[0.75, 0.80, 0.85, 0.90], yticks=[0, 0.5, 1])
b.legend(loc="lower right", title=r"$L/R$", title_fontsize=7, fontsize=6.4, handlelength=1.35, labelspacing=0.14)

validation = pd.read_csv(DATA / "figure1" / "bulk_marginal_validation.csv")
prediction = pd.read_csv(DATA / "figure1" / "bulk_fixedN_predictions.csv")
alignment = prediction.local_alignment_fixedN.iloc[0]
w = np.linspace(0.748, 0.904, 150)
c.plot(w, alignment + 0 * w, color=C["local_alignment"], ls="--")
c.plot(w, alignment * np.exp(-w * w / 2), color=C["drift"], ls="--")
c.plot(validation.W, validation.alignment_observed, "s", mfc="white", color=C["local_alignment"], ms=3.3)
c.plot(validation.W, validation.drift_observed, "o", mfc="white", color=C["drift"], ms=3.3)
c.set(xlabel=r"$W$", ylabel="directional order", xlim=(0.748, 0.904), ylim=(0.59, 1.035), xticks=[0.75, 0.80, 0.85, 0.90], yticks=[0.6, 0.8, 1])
c.legend(handles=[
    Line2D([], [], color=C["local_alignment"], ls="--", marker="s", mfc="white", ms=3.3, label=r"$A_{\rm loc}$"),
    Line2D([], [], color=C["drift"], ls="--", marker="o", mfc="white", ms=3.3, label=r"$D$"),
], loc="center left", bbox_to_anchor=(0.04, 0.63), fontsize=7.0, handlelength=1.6, handletextpad=0.4, labelspacing=0.25, borderpad=0.08, borderaxespad=0, frameon=False)

for ax in [b, c]:
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
fig.savefig(OUTPUT / "fig1.svg")
save(fig, 1)

width_mm, height_mm = definition["size_mm"]
figa = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
axa = figa.add_axes([0, 0, 1, 1])
draw_schematic(axa, definition, C)
for extension in ["pdf", "svg", "png"]:
    figa.savefig(OUTPUT / f"fig1a.{extension}", dpi=600)
plt.close(figa)
(OUTPUT / "fig1a_geometry_checks.json").write_text(json.dumps(geometry, indent=2) + "\n")
