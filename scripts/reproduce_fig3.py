"""Reproduce manuscript Figure 3 from supplied source data."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from _figure_common import DATA, COLORS as C, apply_style, panel, save

apply_style()
styles = [(64, "o", "-"), (128, "s", "--"), (256, "^", "-."), (512, "D", ":")]
fig, axes = plt.subplots(2, 4, figsize=(178 / 25.4, 92 / 25.4))
fig.subplots_adjust(left=0.064, right=0.990, bottom=0.125, top=0.945, wspace=0.57, hspace=0.59)
a, b, c, d, e, f, g, h = axes.ravel()
for ax, label in zip(axes.ravel(), "abcdefgh"):
    panel(ax, label)

scan = pd.read_csv(DATA / "figure3" / "filtered_W_scan.csv").sort_values("W")
for column, key, marker, linestyle, label in [
    ("S_median", "full", "o", "-", r"$S$"),
    ("S1_median", "g1", "s", "--", r"$S_1$"),
    ("S2_median", "g2", "^", "-.", r"$S_2$"),
]:
    z = scan[scan.W >= 0.75]
    a.plot(z.W, z[column], color=C[key], ls=linestyle, marker=marker, markevery=3, mfc="white", ms=3.3, label=label)
a.set(xlabel=r"$W$", ylabel="largest SCC fraction", xlim=(0.748, 0.904), ylim=(0.0008, 1.7), yscale="log", xticks=[0.75, 0.80, 0.85, 0.90])
a.legend(loc="lower right", ncol=3, fontsize=6.5, handlelength=0.85, columnspacing=0.5, handletextpad=0.25, borderpad=0.1)
a.axvspan(0.7975, 0.815, color=C["reference"], alpha=0.065, zorder=-1)

scale = pd.read_csv(DATA / "figure3" / "low_order_scaling.csv")
z = scale[scale.W.eq(0.805)].sort_values("L")
for k, key, marker, linestyle in [(1, "g1", "s", "--"), (2, "g2", "^", "-.")]:
    med = z[f"median__L_times_S{k}"]
    b.plot(z.L, med, color=C[key], ls=linestyle, marker=marker, mfc="white", label=rf"$k={k}$")
    b.fill_between(z.L, z[f"Q1__L_times_S{k}"], z[f"Q3__L_times_S{k}"], color=C[key], alpha=0.10)
b.set(xlabel=r"$L/R$", ylabel=r"$LS_k/R$", xscale="log", ylim=(0.5, 2.1), yticks=[0.5, 1, 1.5, 2])
b.set_xticks([64, 128, 256, 512], ["64", "128", "256", "512"])
b.text(0.035, 0.91, r"$W=0.805$", transform=b.transAxes, fontsize=6.5)
b.legend(loc="upper right", fontsize=6.5, handlelength=1.1, borderpad=0.1)

depth = pd.read_csv(DATA / "figure3" / "depth_summary.csv").sort_values("L")
for metric, key, marker, linestyle, label in [
    ("xiG_over_L", "full", "o", "-", r"$G$"),
    ("xi1_over_L", "g1", "s", "--", r"$G_1$"),
    ("xi2_over_L", "g2", "^", "-.", r"$G_2$"),
]:
    med = depth["median__" + metric]
    c.errorbar(depth.L, med, yerr=np.array([med - depth["Q1__" + metric], depth["Q3__" + metric] - med]), color=C[key], marker=marker, ls=linestyle, mfc="white", ms=3.3, capsize=1.5, elinewidth=0.8, label=label)
c.set(xlabel=r"$L/R$", ylabel=r"$d_{90}/L$", xscale="log", yscale="log", ylim=(0.0013, 1.6))
c.set_xticks([64, 128, 256, 512], ["64", "128", "256", "512"])
c.legend(loc="center left", bbox_to_anchor=(0.01, 0.55), ncol=3, fontsize=6.5, handlelength=0.85, columnspacing=0.45, handletextpad=0.25, borderpad=0.1)
c.text(0.035, 0.035, r"$W=0.805$", transform=c.transAxes, fontsize=6.5)

direction = pd.read_csv(DATA / "figure3" / "direction_pairs.csv")
summary = pd.read_csv(DATA / "figure3" / "direction_summary.csv")
transforms = ["ORIGINAL", "REVERSED", "ROTATED_90"]
for k, key, marker, linestyle, offset in [(1, "g1", "s", "--", -0.045), (2, "g2", "^", "-.", 0.045)]:
    wide = direction.pivot(index="graph_pair_id", columns="orientation_transform", values=f"DLOC_G{k}")[transforms]
    jitter = np.linspace(-0.085, 0.085, len(wide))
    for j, values in enumerate(wide.to_numpy()):
        d.plot(np.arange(3) + offset + jitter[j], values, color=C[key], lw=0.40, alpha=0.075, zorder=1)
        d.scatter(np.arange(3) + offset + jitter[j], values, s=3.2, color=C[key], alpha=0.35, zorder=2)
    z = summary.set_index("orientation_transform").loc[transforms]
    metric = f"DLOC_G{k}"
    med = z["median__" + metric]
    d.errorbar(np.arange(3) + offset, med, yerr=np.array([med - z["Q1__" + metric], z["Q3__" + metric] - med]), color=C[key], marker=marker, ls=linestyle, mfc="white", ms=3.7, capsize=2, elinewidth=1.1, lw=1.2, zorder=4, label=rf"$k={k}$")
d.set(ylabel=r"$\mathrm{DLOC}_k$", ylim=(-0.065, 1.08), xlim=(-0.38, 2.38), yticks=[0, 0.5, 1], xticks=[0, 1, 2])
d.set_xticklabels(["Original\nbottom", "Reversed\ntop", r"$+90^\circ$" + "\nright"], fontsize=6.2, linespacing=1.15)
d.xaxis.set_minor_locator(NullLocator())
d.legend(loc="center right", bbox_to_anchor=(1.01, 0.48), fontsize=6.5, handlelength=1.1, borderpad=0.1)

profile = pd.read_csv(DATA / "figure3" / "profile_occupancy_summary.csv")
for ax, k, key in [(e, 1, "g1"), (f, 2, "g2")]:
    for L, marker, linestyle in styles:
        z = profile[(profile.L.eq(L)) & (profile.component.eq(f"G{k}"))].sort_values("d_over_R_bin_left")
        ax.plot((z.d_over_R_bin_left + z.d_over_R_bin_right) / 2, z.median_occupancy, color=C[key], ls=linestyle, marker=marker, mfc="white", markevery=10, ms=3, label=str(L))
    ax.set(xlabel=r"$d/R$", ylabel=rf"$p_{k}(d/R)$", xlim=(0, 5), ylim=(0, 1.03), xticks=[0, 2, 4], yticks=[0, 0.5, 1])
e.legend(loc="upper right", ncol=1, title=r"$L/R$", fontsize=6.2, title_fontsize=6.5, handlelength=1.4, labelspacing=0.15, borderpad=0.1)
f.text(0.965, 0.91, r"$W=0.805$", transform=f.transAxes, fontsize=6.5, ha="right")

rec = pd.read_csv(DATA / "figure3" / "internal_recurrence_summary.csv")
for ax, quantity, ylabel in [(g, "f_internal_rec", r"$f^{\rm int}_{k,\rm rec}$"), (h, "eta_internal", r"$\eta_k^{\rm int}$")]:
    for k, key, marker, linestyle in [(1, "g1", "s", "--"), (2, "g2", "^", "-.")]:
        z = rec[(rec.k.eq(k)) & (rec.quantity.eq(quantity))].sort_values("W")
        ax.plot(z.W, z["median"], color=C[key], ls=linestyle, marker=marker, mfc="white", markevery=3, ms=3.3, label=rf"$k={k}$")
        ax.fill_between(z.W, z.Q25, z.Q75, color=C[key], alpha=0.10)
    ax.set(xlabel=r"$W$", ylabel=ylabel, xlim=(0.748, 0.904), xticks=[0.75, 0.80, 0.85, 0.90])
    ax.axvspan(0.7975, 0.815, color=C["reference"], alpha=0.065, zorder=-1)
g.set_ylim(0, 1.035)
g.legend(loc="upper left", fontsize=6.5, handlelength=1.1)
h.set(yscale="log", ylim=(0.0007, 1.5))

save(fig, 3)
