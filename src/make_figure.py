"""Draw Fig. 1: difference-in-differences estimates by horizon, by journal type and one year earlier (A),
cumulative PMC availability with 95% bands outside fully open-access journals (B), and a timeline of the two
deposit routes before and after the policy (C).

Availability is dated by the day each article went live in PMC. Every label in panel C is read from results/.

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

DATA = Path("results")
FIGURES = Path("figures")
CM = 1 / 2.54
NIH_COLOR = "#0072B2"  # Okabe-Ito blue: NIH-funded, in every panel
COMPARISON_COLOR = "#7F7F7F"  # comparison articles
NIH_2025_COLOR = "#6FA8D2"  # the same hue, lighter: before the policy
COMPARISON_2025_COLOR = "#A8A8A8"
PRIMARY_HORIZON = 7
MAX_DAY = 120
HORIZONS = (0, 7, 30, 60, 90, 120)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
        "font.size": 8.5,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.bf": "Arial:bold",
        "axes.labelsize": 8.5,
        "axes.titlesize": 8.5,
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
    """Rows for panel A: (label, (est, lo, hi)) or (group header, None)."""
    did = pd.read_csv(DATA / "did_by_horizon.csv").set_index(["stratum", "horizon", "outcome"])
    numbers = json.loads((DATA / "numbers.json").read_text())

    def row(h: int, stratum: str = "all") -> tuple:
        live = did.loc[(stratum, h, "live")]
        return (live["est"], live["lo"], live["hi"])

    def placebo(h: int) -> tuple:
        r = numbers[f"placebo_same_journals_{h}d_live"]
        return (r["est"], r["lo"], r["hi"])

    labels = {
        0: "By publication",
        7: "Within 7 days",
        30: "Within 30 days",
        60: "Within 60 days",
        90: "Within 90 days",
        120: "Within 120 days",
    }
    return (
        [("2026 vs 2025: policy in effect", None)]
        + [(labels[h], row(h)) for h in (120, 30, 7)]
        + [
            ("120-day effect by journal type", None),
            ("Not fully OA", row(120, "non_doaj")),
            ("Fully OA, 88 to 94% already free", row(120, "doaj")),
            ("2025 vs 2024: no policy yet", None),
            ("Within 120 days", placebo(120)),
            ("Within 7 days", placebo(7)),
        ]
    )


def draw_curves(ax: Axes, curves: pd.DataFrame, in_doaj: bool, title: str) -> None:
    subset = curves[curves["in_doaj"] == in_doaj]
    arms = (  # year, nih, color, linestyle, width: the 2026 NIH-funded curve is the only heavy line
        (2025, False, COMPARISON_2025_COLOR, "--", 0.9),
        (2025, True, NIH_2025_COLOR, "--", 1.1),
        (2026, False, COMPARISON_COLOR, "-", 0.9),
        (2026, True, NIH_COLOR, "-", 1.6),
    )
    for year, nih, color, linestyle, width in arms:
        arm = subset[(subset["year"] == year) & (subset["nih"] == nih)].sort_values("day")
        group = "NIH-funded" if nih else "Not NIH-funded"
        ax.fill_between(
            arm["day"], 100 * arm["ci_lower"], 100 * arm["ci_upper"], color=color, alpha=0.12, lw=0
        )
        ax.plot(arm["day"], 100 * arm["share"], linestyle, color=color, lw=width, label=f"{group}, {year}")
    ax.set_title(title, loc="left", pad=4)
    ax.set_xlim(0, MAX_DAY)
    ax.set_ylim(0, 100)
    ax.set_xticks([0, 30, 60, 90, 120])
    ax.set_xlabel("Days since publication")
    ax.axvline(PRIMARY_HORIZON, color="#BBBBBB", lw=0.5, ls=":", zorder=0)
    ax.text(PRIMARY_HORIZON + 1.5, 97, "7 days", ha="left", va="top", fontsize=7.5)


def annotate_shares(ax: Axes, curves: pd.DataFrame) -> None:
    """Label all four 120-day shares in the non-DOAJ panel; the DID at 120 days is their double difference."""
    day = curves[(~curves["in_doaj"]) & (curves["day"] == MAX_DAY)]
    colors = {
        (2026, True): NIH_COLOR,
        (2025, True): NIH_2025_COLOR,
        (2026, False): COMPARISON_COLOR,
        (2025, False): COMPARISON_2025_COLOR,
    }
    shares = {
        (y, n): 100 * day.loc[(day["year"] == y) & (day["nih"] == n), "share"].item() for y, n in colors
    }
    low, high = sorted((shares[(2025, False)], shares[(2026, False)]))
    labels = (  # the two comparison curves share one label
        (shares[(2026, True)], f"{shares[(2026, True)]:.0f}%"),
        (shares[(2025, True)], f"{shares[(2025, True)]:.0f}%"),
        ((low + high) / 2, f"{low:.0f} to {high:.0f}%"),
    )
    for y, text in labels:
        ax.text(MAX_DAY + 2, y, text, va="center", fontsize=7.5, color="black", clip_on=False)


def load_schematic() -> dict[str, float | str]:
    """Panel C values: matched NIH-funded articles before (2025) and after (2026) the policy."""
    n = json.loads((DATA / "numbers.json").read_text())
    did = pd.read_csv(DATA / "did_by_horizon.csv").set_index(["stratum", "horizon", "outcome"])
    route = pd.read_csv(DATA / "route_decomposition.csv").set_index(["stratum", "h", "route"])
    share = route.loc[("all", 120, "manuscript"), "est"] / did.loc[("all", 120, "live"), "est"]
    assert 0.6 <= share < 0.7, share  # "two thirds of the gain"
    return {
        "final_2025": n["lag_live_median_2025_final"],
        "final_2026": n["lag_live_median_2026_final"],
        "ms_2025": n["ms25_lag_live_q50"],
        "ms_2026": n["ms26_lag_live_q50"],
        "accepted_day": -n["acc_acc_to_pub_q50"],
    }


def draw_schematic(ax: Axes, v: dict[str, float | str]) -> None:
    """Timeline: the policy asks for day 0; each route before (2025) and after (2026) the policy."""
    a_max = 390  # panel C runs to a year so the 374-day median is on the axis
    ax.set_xlim(-150, a_max)
    ax.set_ylim(0, 3.2)
    ax.set_yticks([])
    for side in ("left", "top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_bounds(v["accepted_day"], a_max)
    ax.set_xticks([0, 100, 200, 300, 365])
    ax.set_xticklabels(["0", "100", "200", "300", "1 year"])
    ax.set_xlabel("Days after publication", labelpad=1)
    ax.plot(v["accepted_day"], 0, "|", color="black", ms=4, mew=0.5, clip_on=False, zorder=3)
    ax.text(v["accepted_day"] - 3, 0.1, "Accepted", ha="right", va="bottom", fontsize=7.5)
    ax.text(
        v["accepted_day"] - 3,
        0.02,
        "(author's deposit due)",
        ha="right",
        va="top",
        fontsize=7.5,
        clip_on=False,
    )
    ax.text(4, 0.1, "Published", ha="left", va="bottom", fontsize=7.5)
    ax.axvline(0, color="#999999", lw=0.6, zorder=2)
    ax.text(
        1.5,
        3.15,
        "NIH policy since July 2025: free to read on publication day",
        ha="left",
        va="top",
        fontsize=7.5,
    )
    routes = (  # label, y of the 2026 line, (2025 median, 2026 median)
        (r"Deposited by $\bf{publisher}$" + "\n(final version)", 2.0, v["final_2025"], v["final_2026"]),
        (r"Deposited by $\bf{author}$" + "\n(accepted version)", 0.75, v["ms_2025"], v["ms_2026"]),
    )
    for label, y, before, after in routes:
        ax.text(-20, y + 0.2, label, ha="right", va="center", fontsize=7.5, linespacing=1.15)
        for day, dy, color, style, tag in (
            (before, 0.4, NIH_2025_COLOR, "--", "2025"),
            (after, 0.0, NIH_COLOR, "-", "2026"),
        ):
            ax.plot([v["accepted_day"], day], [y + dy, y + dy], style, color=color, lw=1.6, zorder=3)
            ax.plot(day, y + dy, "o", color=color, ms=5, zorder=4)
            note = f"{tag}: median {day:.0f} days"
            if y < 1 and dy == 0:
                note += ", down from about a year"
            if day > a_max - 60:  # label to the left of a dot near the right edge
                ax.text(day - 6, y + dy + 0.32, note, ha="right", va="bottom", fontsize=7.5)
            else:
                ax.text(day + 6, y + dy, note, ha="left", va="center", fontsize=7.5)


def draw_estimates(ax: Axes, rows: list[tuple[str, tuple[float, float, float] | None]]) -> None:
    placebo_header = next(i for i, (label, est) in enumerate(rows) if label.startswith("2025 vs 2024"))
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
        estimate, lower, upper = est
        is_placebo = index > placebo_header
        color = COMPARISON_2025_COLOR if is_placebo else NIH_COLOR
        ax.plot([lower, upper], [position, position], color=color, lw=1.2, solid_capstyle="butt", zorder=4)
        ax.plot(estimate, position, "|", color=color, ms=5, mew=0.8, zorder=5)  # tick: narrow CIs visible
        ax.text(-0.03, position, label, transform=y_axis, ha="right", va="center", fontsize=8)
    ax.axvline(0, color="#999999", lw=0.6)
    ax.set_yticks([])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlim(-5, 36)
    ax.set_xticks([0, 10, 20, 30])
    ax.set_xlabel("Policy effect (percentage points)")
    ax.spines["left"].set_visible(False)
    ax.text(  # title starts to the right of the panel letter, above the label column
        -0.58,
        1.04,
        "Policy effect: one more paper in five free\nby 120 days, none detectable within 7 days",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=plt.rcParams["axes.titlesize"],
        clip_on=False,
    )


def main() -> None:
    curves = load_curves()
    rows = load_estimates()

    fig = plt.figure(figsize=(17.8 * CM, 9.2 * CM))
    top = gridspec.GridSpec(
        1, 2, width_ratios=[0.72, 1], wspace=0.55, left=0.25, right=0.92, top=0.89, bottom=0.60
    )
    bottom = gridspec.GridSpec(1, 1, left=0.07, right=0.955, top=0.40, bottom=0.11)
    ax_c, ax_a = (fig.add_subplot(top[i]) for i in range(2))
    ax_s = fig.add_subplot(bottom[0])

    draw_schematic(ax_s, load_schematic())
    ax_s.set_title(
        "Two routes to free access for NIH-funded papers, and both still arrived weeks after publication",
        loc="left",
        pad=4,
    )
    draw_curves(
        ax_a,
        curves,
        in_doaj=False,
        title=("NIH papers free to read by 120 days rose,\nothers barely moved (fully OA journals excluded)"),
    )
    ax_a.set_ylabel("Free to read in PMC (%)")
    lines = {line.get_label(): line for line in ax_a.get_lines()}
    order = [
        "NIH-funded, 2026",
        "NIH-funded, 2025",
        "Not NIH-funded, 2026",
        "Not NIH-funded, 2025",
    ]
    ax_a.legend(
        [lines[k] for k in order],
        order,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.22, 1.06),
        handlelength=1.8,
        labelspacing=0.2,
    )
    annotate_shares(ax_a, curves)

    draw_estimates(ax_c, rows)

    # Panel letters in figure coordinates: A and C share a left edge, A and B share a baseline.
    top = max(ax_c.get_position().y1, ax_a.get_position().y1) + 0.085
    for label, x, y in (
        ("A", 0.015, top),
        ("B", ax_a.get_position().x0 - 0.075, top),
        ("C", 0.015, ax_s.get_position().y1 + 0.075),
    ):
        fig.text(x, y, label, fontsize=11, fontweight="bold", va="top")

    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "fig1.pdf", metadata={"CreationDate": None})
    fig.savefig(FIGURES / "fig1.png")


if __name__ == "__main__":
    main()
