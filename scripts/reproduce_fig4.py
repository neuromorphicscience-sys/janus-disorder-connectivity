"""Reproduce manuscript Figure 4 from supplied source data."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _figure_common import DATA, COLORS as C, apply_style, panel, save

apply_style()
fig = plt.figure(figsize=(86 / 25.4, 94 / 25.4))
a = fig.add_axes([0.17, 0.60, 0.34, 0.32])
b = fig.add_axes([0.66, 0.60, 0.31, 0.32])
c = fig.add_axes([0.17, 0.13, 0.80, 0.31])
for ax, label in zip([a, b, c], "abc"):
    panel(ax, label)

profile = pd.read_csv(DATA / "figure4" / "boundary_depth_prediction.csv")
a.plot(profile.d_over_R, profile.mean_Z_up_fixedN_x_averaged, color=C["benchmark"], lw=1.8)
a.text(0.93, 0.88, "theory", transform=a.transAxes, ha="right", fontsize=7)
a.set(xlabel=r"$d/R$", ylabel=r"$\langle Z_\uparrow(d)\rangle$", xlim=(0, 0.22), ylim=(-0.1, 7.3), xticks=[0, 0.1, 0.2], yticks=[0, 3, 6])

clean = pd.read_csv(DATA / "figure4" / "clean_feedback_counts.csv")
prediction = json.loads((DATA / "figure4" / "clean_feedback_prediction.json").read_text())
values = np.sort(clean.F_count.to_numpy())
counts, positions = {}, []
for value in values:
    j = counts.get(value, 0)
    positions.append(((-1) ** j) * (j // 2 + 1) * 0.035)
    counts[value] = j + 1
b.scatter(np.array(positions) - 0.13, values, s=10, facecolors="white", edgecolors=C["full"], lw=0.7)
b.errorbar(0.21, values.mean(), yerr=values.std(ddof=1) / np.sqrt(len(values)), fmt="o", color=C["full"], ms=3.8, capsize=2, lw=1.1)
b.axhline(prediction["fixed_geometric_integral"], color=C["benchmark"], lw=1.8, ls="--")
b.text(0.95, 0.1, "theory", transform=b.transAxes, ha="right", color=C["neutral"], fontsize=7)
b.set(xlim=(-0.36, 0.36), ylim=(870, 935), ylabel=r"$|F|$", xticks=[0])
b.set_xticklabels([r"$W=0\ (n=24)$"], fontsize=7)

deletion = pd.read_csv(DATA / "figure4" / "boundary_deletion_checks.csv")
z = deletion[(deletion.cut_kind.eq("all_sides")) & (deletion.removed_width.eq(4))].sort_values("W")
assert len(z) == 4
for j, (_, row) in enumerate(z.iterrows()):
    c.plot([j, j], [row.G_fraction_remaining, row.G2_fraction_remaining], color=C["connector"], lw=1, zorder=1)
c.plot(range(4), z.G_fraction_remaining, "o", color=C["full"], ms=4, mfc="white", label=r"$G[I]$")
c.plot(range(4), z.G2_fraction_remaining, "^", color=C["g2"], ms=4, mfc="white", label=r"$G_2[I]$")
c.set(yscale="log", ylabel=r"largest SCC / $|I|$", xlabel=r"realization ($W$)", xlim=(-0.4, 3.4), ylim=(0.00085, 1.7), xticks=range(4))
c.set_xticklabels([".805", ".815", ".860", ".900"])
c.legend(loc="lower right", ncol=2, fontsize=7, handlelength=1.0, columnspacing=1.0)

save(fig, 4)
