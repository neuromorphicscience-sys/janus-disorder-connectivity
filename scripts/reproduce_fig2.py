"""Reproduce manuscript Figure 2 from supplied source data."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from _figure_common import DATA, COLORS as C, apply_style, panel, save

apply_style()
fig, axes = plt.subplots(2, 2, figsize=(86 / 25.4, 91 / 25.4))
fig.subplots_adjust(left=0.155, right=0.98, bottom=0.12, top=0.94, wspace=0.55, hspace=0.56)
a, b, c, d = axes.ravel()
for ax, label in zip(axes.ravel(), "abcd"):
    panel(ax, label)

pairs = pd.read_csv(DATA / "figure2" / "paired_reassignment.csv").sort_values("pair_id")
assert len(pairs) == 48
jitter = np.linspace(-0.06, 0.06, 48)
for j, (_, row) in enumerate(pairs.iterrows()):
    a.plot([jitter[j], 1 + jitter[j]], [row.original_S, row.reassigned_S], color=C["connector"], lw=0.55, alpha=0.75, zorder=1)
a.scatter(jitter, pairs.original_S, s=8, c=C["full"], zorder=3)
a.scatter(1 + jitter, pairs.reassigned_S, s=8, facecolors="white", edgecolors=C["reassigned"], lw=0.7, zorder=3)
a.set(yscale="log", ylim=(5e-5, 1.3), xlim=(-0.22, 1.22), ylabel=r"$S$", xticks=[0, 1])
a.set_xticklabels(["Original", "Reassigned"], fontsize=6.4)
a.text(0.06, 0.09, "48/48 pairs", transform=a.transAxes, fontsize=7)

metrics = [(r"$S$", "S"), (r"$A_{\rm loc}$", "A_local"), (r"$D$", "D_global"), (r"$f_\uparrow$", "upstream_fraction")]
for j, (label, column) in enumerate(metrics):
    values = np.log10(pairs["reassigned_" + column] / pairs["original_" + column])
    q = np.quantile(values, [0.25, 0.5, 0.75])
    y = 3 - j
    color = C[["full", "local_alignment", "drift", "control"][j]]
    b.scatter(values, np.full(48, y) + np.linspace(-0.13, 0.13, 48), s=2.5, c=color, alpha=0.22)
    b.plot(q[[0, 2]], [y, y], color=color, lw=1.8)
    b.plot(q[1], y, "o", color=color, ms=3.5)
b.axvline(0, color=C["reassigned"], ls=":", lw=0.7)
b.set(xlabel=r"$\log_{10}(X_r/X_o)$", xlim=(-4.7, 0.45), ylim=(-0.5, 3.5), xticks=[-4, -2, 0], yticks=[3, 2, 1, 0])
b.set_yticklabels([x[0] for x in metrics])
b.yaxis.set_minor_locator(NullLocator())

dy = pd.read_csv(DATA / "figure2" / "displacement_histogram.csv")
for state, key, linestyle in [("original", "full", "-"), ("rearranged", "reassigned", "--")]:
    z = dy[dy.state.eq(state)].groupby(["bin_left", "bin_right"])["count"].sum().reset_index().sort_values("bin_right")
    x = np.r_[z.bin_left.iloc[0], z.bin_right]
    y = np.r_[0, z["count"].cumsum() / z["count"].sum()]
    c.step(x, y, where="post", color=C[key], ls=linestyle, label="Original" if state == "original" else "Reassigned", lw=1.25)
c.set(xlabel=r"$\Delta y/R$", ylabel="binned CDF", xlim=(-2, 2), ylim=(0, 1.025), xticks=[-2, 0, 2], yticks=[0, 0.5, 1])
c.legend(loc="lower right", fontsize=6.2, handlelength=1.2, labelspacing=0.25)

length = pd.read_csv(DATA / "figure2" / "edge_length_quantiles.csv")
probs = np.array([0.10, 0.25, 0.50, 0.75, 0.90, 0.95])
qcols = [f"q{int(v * 100)}" for v in probs]
keys = ["pair_id", "W"]
original = length[length.state.eq("original")][keys + qcols]
reassigned = length[length.state.eq("rearranged")][keys + qcols]
paired = original.merge(reassigned, on=keys, validate="one_to_one", suffixes=("_o", "_r")).sort_values(keys)
deltas = np.column_stack([paired[name + "_r"] - paired[name + "_o"] for name in qcols])
quartiles = np.quantile(deltas, [0.25, 0.5, 0.75], axis=0)
d.axhline(0, color=C["reference"], ls=":", lw=0.8, zorder=1)
for j, u in enumerate(probs):
    d.scatter(np.full(len(paired), u), deltas[:, j], s=2.3, color=C["control"], alpha=0.24, zorder=2)
d.errorbar(probs, quartiles[1], yerr=np.array([quartiles[1] - quartiles[0], quartiles[2] - quartiles[1]]), fmt="o-", color=C["control"], ms=3.4, mfc="white", capsize=2, elinewidth=1.25, lw=1.1, zorder=4)
d.set(xlabel="quantile level $u$", ylabel=r"$\Delta\ell_u/R$", xlim=(0.06, 0.99), ylim=(-0.0016, 0.0016), xticks=[0.1, 0.5, 0.9], yticks=[-0.001, 0, 0.001])
d.ticklabel_format(axis="y", style="sci", scilimits=(-3, -3), useMathText=True)

save(fig, 2)
