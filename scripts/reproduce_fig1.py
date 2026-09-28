"""Reproduce manuscript Figure 1 from supplied source data."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
from PIL import Image
from _figure_common import DATA, COLORS as C, apply_style, panel, save

apply_style()
plt.rcParams.update({"legend.fontsize": 7, "lines.linewidth": 1.15, "lines.markersize": 3.5})
styles = [(64, "o", "-"), (128, "s", "--"), (256, "^", "-."), (512, "D", ":")]

reference_path = DATA / "figure1" / "schematic_source.png"
reference = np.asarray(Image.open(reference_path).convert("RGB"))
reference_display = reference.copy()
rgb = reference.astype(np.int16)
blue = (rgb[..., 2] - rgb[..., 0] > 8) & (rgb[..., 1] - rgb[..., 0] > 4) & (rgb[..., 2] > 100)
yy, xx = np.indices(blue.shape)
removed = blue & (
    ((xx >= 195) & (xx <= 330) & (yy >= 165) & (yy <= 365))
    | ((xx >= 1540) & (xx <= 1870) & (yy >= 325) & (yy <= 525))
)
reference_display[removed] = np.array([253, 253, 253], dtype=np.uint8)
ref_h, ref_w = reference.shape[:2]
canvas_w = 1800.0
canvas_h = canvas_w * ref_h / ref_w

fig = plt.figure(figsize=(86 / 25.4, 90 / 25.4))
a = fig.add_axes([1 / 86, 52 / 90, 84 / 86, (84 * ref_h / ref_w) / 90])
b = fig.add_axes([0.155, 9.96 / 90, 0.355, 33.615 / 90])
c = fig.add_axes([0.665, 9.96 / 90, 0.315, 33.615 / 90])
a.imshow(reference_display, extent=(0, canvas_w, canvas_h, 0), interpolation="none", zorder=0)
a.set(xlim=(0, canvas_w), ylim=(canvas_h, 0))
a.set_aspect("equal")
a.axis("off")
masks = [
    (0, 5, 58, 84), (0, 125, 395, 525), (406, 251, 603, 347),
    (775, 65, 925, 129), (1145, 292, 1265, 351), (151, 533, 269, 587),
    (110, 590, 311, 652), (744, 590, 985, 652), (1416, 590, 1725, 652),
]
for x0, y0, x1, y1 in masks:
    a.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="#FDFDFD", edgecolor="none", zorder=1))
a.text(17, 28, "(a)", ha="left", va="top", fontsize=8.5, fontweight="bold", zorder=2)
source = (205, 330)
selected = [(309, 246), (351, 330), (302, 430)]
unselected = [(66, 216), (38, 302), (69, 416), (157, 476)]
a.add_patch(Circle(source, 184, fill=False, edgecolor="#AEB4BA", linewidth=0.85, linestyle=(0, (4, 4)), zorder=2))
for target in unselected:
    a.plot([source[0], target[0]], [source[1], target[1]], color="#BCC1C6", linewidth=0.65, linestyle=(0, (3, 3)), zorder=2)
for target in selected:
    a.add_patch(FancyArrowPatch(source, target, arrowstyle="-|>", mutation_scale=6.5, linewidth=1.15, color=C["neutral"], zorder=3, shrinkA=2.2, shrinkB=1.8))
for target in unselected + selected:
    a.add_patch(Circle(target, 12, facecolor="#F2F2F2", edgecolor="#555A5E", linewidth=0.85, zorder=4))
a.add_patch(Circle(source, 16, facecolor=C["neutral"], edgecolor=C["neutral"], zorder=5))
a.add_patch(FancyArrowPatch((145, 198), (278, 198), arrowstyle="-|>", mutation_scale=9, linewidth=2.0, color=C["direction"], zorder=4))
a.text(211, 168, r"$\mathbf{e}_i$", ha="center", va="center", fontsize=6.5, color=C["direction"], zorder=4)
a.text(502, 298, "same\ntop-$q$ rule", ha="center", va="center", fontsize=7, linespacing=1.12, zorder=2)
a.text(848, 99, r"$W=0$", ha="center", va="center", fontsize=7.5, zorder=2)
a.text(1205, 324, r"$W\uparrow$", ha="center", va="center", fontsize=7.5, zorder=2)
a.text(210, 560, r"$q=3$", ha="center", va="center", fontsize=7, zorder=2)
a.text(210, 620, "local rule", ha="center", va="center", fontsize=7.3, zorder=2)
a.text(860, 620, "fragmented", ha="center", va="center", fontsize=7.3, zorder=2)
a.text(1570, 620, "restored GSCC", ha="center", va="center", fontsize=7.3, zorder=2)

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
], loc="center left", bbox_to_anchor=(0.04, 0.63), fontsize=7.0, handlelength=1.6, handletextpad=0.4, labelspacing=0.25, borderpad=0.08, frameon=False)

save(fig, 1)
