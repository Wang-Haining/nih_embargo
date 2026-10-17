"""Draw Fig. 1: cumulative PMC availability with 95% bands (A, B) and policy estimates by horizon (C).

Availability is dated by the day each article went live in PMC. Panel C also shows, as open markers,
the estimates obtained with the release date PMC reports.

Run from the repository root:  python src/make_figure.py
Inputs are the tables in results/ written by src/estimate.py.
Outputs are figures/fig1.pdf and figures/fig1.png.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

DATA = Path("results")
FIGURES = Path("figures")
CM = 1 / 2.54
NIH_COLOR = "#0072B2"  # Okabe-Ito blue
COMPARISON_COLOR = "#7F7F7F"
PRIMARY_HORIZON = 7
MAX_DAY = 120
HORIZONS = (0, 7, 30, 60, 90, 120)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 600,
    }
)


def load_curves() -> pd.DataFrame:
    """Cumulative share in PMC by day, with journal-cluster bootstrap bands."""
    curves = pd.read_csv(DATA / "curves.csv")
    return curves[curves["day"] <= MAX_DAY]


def load_estimates() -> list[tuple[str, tuple | None]]:
    """Rows for panel C: (label, (live est, lo, hi, release est or None)) or (group header, None)."""
    did = pd.read_csv(DATA / "did_by_horizon.csv")
    did = did[did["stratum"] == "all"].set_index(["horizon", "outcome"])
    numbers = json.loads((DATA / "numbers.json").read_text())

    def row(h: int) -> tuple:
        live = did.loc[(h, "live")]
        return (live["est"], live["lo"], live["hi"], did.loc[(h, "rel"), "est"])

    def placebo(h: int) -> tuple:
        r = numbers[f"placebo_same_journals_{h}d_live"]
        return (r["est"], r["lo"], r["hi"], None)

    labels = {
        0: "By publication",
        7: "Within 7 days",
        30: "Within 30 days",
        60: "Within 60 days",
        90: "Within 90 days",
        120: "Within 120 days",
    }
    return (
        [("Policy, 2025 vs 2026", None)]
        + [(labels[h], row(h)) for h in HORIZONS]
        + [
            ("Pre-policy, 2024 vs 2025", None),
            ("Within 7 days", placebo(7)),
            ("Within 120 days", placebo(120)),
        ]
    )


def draw_curves(ax: Axes, curves: pd.DataFrame, in_doaj: bool, title: str) -> None:
    subset = curves[curves["in_doaj"] == in_doaj]
    for year, linestyle in ((2025, "--"), (2026, "-")):
        for nih, color, width in ((False, COMPARISON_COLOR, 0.9), (True, NIH_COLOR, 1.3)):
            arm = subset[(subset["year"] == year) & (subset["nih"] == nih)].sort_values("day")
            group = "NIH-supported" if nih else "Comparison"
            ax.fill_between(
                arm["day"], 100 * arm["ci_lower"], 100 * arm["ci_upper"], color=color, alpha=0.12, lw=0
            )
            ax.plot(
                arm["day"], 100 * arm["share"], linestyle, color=color, lw=width, label=f"{group}, {year}"
            )
    ax.set_title(title, loc="left", pad=3)
    ax.set_xlim(0, MAX_DAY)
    ax.set_ylim(0, 100)
    ax.set_xticks([0, 30, 60, 90, 120])
    ax.set_xlabel("Days since publication")
    ax.axvline(PRIMARY_HORIZON, color="#BBBBBB", lw=0.5, ls=":", zorder=0)
    ax.text(PRIMARY_HORIZON + 1.5, 97, "7 d", ha="left", va="top", fontsize=7.5, color="#888888")


def annotate_nih_shares(ax: Axes, curves: pd.DataFrame) -> None:
    """Label the NIH-supported shares at 120 days, 2025 and 2026, in the non-DOAJ panel."""
    day = curves[(~curves["in_doaj"]) & (curves["day"] == MAX_DAY) & curves["nih"]]
    for year in (2026, 2025):
        share = 100 * day.loc[day["year"] == year, "share"].item()
        ax.annotate(
            f"{share:.0f}%",
            xy=(MAX_DAY, share),
            xytext=(MAX_DAY + 2, share),
            va="center",
            fontsize=8,
            color=NIH_COLOR,
            annotation_clip=False,
        )


def draw_estimates(ax: Axes, rows: list[tuple[str, tuple[float, float, float] | None]]) -> None:
    placebo_header = next(i for i, (label, est) in enumerate(rows) if label.startswith("Pre-policy"))
    for position, (index, (label, est)) in enumerate(reversed(list(enumerate(rows)))):
        y_axis = ax.get_yaxis_transform()
        if est is None:
            ax.text(
                -0.03,
                position,
                label,
                transform=y_axis,
                ha="right",
                va="center",
                fontsize=8,
                fontweight="bold",
            )
            continue
        estimate, lower, upper, release = est
        is_placebo = index > placebo_header
        color = "#999999" if is_placebo else "#333333"
        ax.plot([lower, upper], [position, position], color=color, lw=1.2, solid_capstyle="butt", zorder=4)
        ax.plot(estimate, position, "|", color=color, ms=5, mew=0.8, zorder=5)  # tick: narrow CIs visible
        if release is not None:
            ax.plot(release, position, "o", mfc="none", mec="#999999", mew=0.8, ms=3.4, zorder=3)
        ax.text(-0.03, position, label, transform=y_axis, ha="right", va="center", fontsize=8)
    ax.axvline(0, color="#999999", lw=0.6)
    ax.set_yticks([])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlim(-5, 27)
    ax.set_xticks([0, 10, 20])
    ax.set_xlabel("Percentage points")
    ax.spines["left"].set_visible(False)
    ax.set_title("Estimates", loc="left", pad=3)
    key = len(rows)  # marker key on its own row above the "Policy" header (rows are drawn bottom-up)
    ax.set_ylim(-0.6, len(rows) + 0.6)
    handles = [
        Line2D([], [], ls="none", marker="|", color="#333333", ms=5, mew=0.8, label="go-live"),
        Line2D(
            [], [], ls="none", marker="o", mfc="none", mec="#999999", mew=0.8, ms=3.4, label="release date"
        ),
    ]
    ax.legend(  # key in the label column, right-aligned like the row labels
        handles=handles,
        loc="center right",
        bbox_to_anchor=(-0.01, key),
        bbox_transform=ax.get_yaxis_transform(),
        ncol=2,
        frameon=False,
        fontsize=7,
        handlelength=0.8,
        handletextpad=0.3,
        columnspacing=0.9,
        borderaxespad=0,
        borderpad=0,
    )


def main() -> None:
    curves = load_curves()
    rows = load_estimates()

    fig = plt.figure(figsize=(17.8 * CM, 7.0 * CM))
    grid = gridspec.GridSpec(
        1, 3, width_ratios=[1, 1, 0.75], wspace=0.25, left=0.07, right=0.955, top=0.90, bottom=0.18
    )
    ax_a, ax_b, ax_c = (fig.add_subplot(grid[i]) for i in range(3))

    draw_curves(ax_a, curves, in_doaj=False, title="Journals not in DOAJ")
    draw_curves(ax_b, curves, in_doaj=True, title="Journals in DOAJ")
    ax_a.set_ylabel("Articles in PubMed Central (%)")
    ax_b.set_yticklabels([])
    ax_a.legend(frameon=False, loc="upper right", handlelength=2.2)
    annotate_nih_shares(ax_a, curves)

    draw_estimates(ax_c, rows)
    box = ax_c.get_position()
    ax_c.set_position([box.x0 + 0.12, box.y0, box.width - 0.12, box.height])

    for ax, label, dx in ((ax_a, "A", -0.20), (ax_b, "B", -0.10), (ax_c, "C", -1.05)):
        ax.text(dx, 1.10, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "fig1.pdf", metadata={"CreationDate": None})
    fig.savefig(FIGURES / "fig1.png")


if __name__ == "__main__":
    main()
