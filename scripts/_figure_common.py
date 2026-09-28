"""Shared manuscript-figure style and paths."""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUT = ROOT / "figures" / "generated"
OUTPUT.mkdir(parents=True, exist_ok=True)

COLORS = {
    "full": "#D81B60",
    "g1": "#004D40",
    "g2": "#1E88E5",
    "benchmark": "#FFC107",
    "reassigned": "#777777",
    "control": "#777777",
    "local_alignment": "#004D40",
    "drift": "#1E88E5",
    "direction": "#1E88E5",
    "neutral": "#111111",
    "reference": "#777777",
    "connector": "#D8D8D8",
}


def apply_style():
    mpl.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "lines.linewidth": 1.15,
            "lines.markersize": 3.5,
            "text.color": COLORS["neutral"],
            "axes.edgecolor": COLORS["neutral"],
            "axes.labelcolor": COLORS["neutral"],
            "xtick.color": COLORS["neutral"],
            "ytick.color": COLORS["neutral"],
            "savefig.transparent": False,
            "savefig.facecolor": "white",
        }
    )


def box_axes(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)


def panel(ax, label):
    ax.text(
        -0.17,
        1.04,
        f"({label})",
        transform=ax.transAxes,
        fontsize=8.5,
        fontweight="bold",
        va="bottom",
    )


def save(fig, number):
    for ax in fig.axes:
        if ax.axison:
            box_axes(ax)
    pdf = OUTPUT / f"fig{number}.pdf"
    png = OUTPUT / f"fig{number}.png"
    fig.savefig(pdf, dpi=600)
    fig.savefig(png, dpi=350)
    plt.close(fig)
    print(pdf)
    print(png)
