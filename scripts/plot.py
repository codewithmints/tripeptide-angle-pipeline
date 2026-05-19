import argparse
import logging
import os
import sys

import matplotlib
matplotlib.use("Agg")                      
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

try:
    import seaborn as sns
except ImportError:
    sys.exit("[ERROR] seaborn is not installed. Run: pip install seaborn")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── Visual design ─────────────────────────────────────────────────────────────
CATEGORY_ORDER = ["Tiny", "Small", "Intermediate", "Large", "Bulky"]

PALETTE = {
    "Tiny":         "#4C9BE8",   # steel blue
    "Small":        "#54B472",   # medium green
    "Intermediate": "#E8A23C",   # amber
    "Large":        "#9B59B6",   # purple
    "Bulky":        "#E05C5C",   # coral red
}

LINE_STYLES = {
    "Tiny":         "solid",
    "Small":        "dashed",
    "Intermediate": "dashdot",
    "Large":        "dotted",
    "Bulky":        (0, (5, 1)),  # densely dashed
}

# Minimum samples required to draw a KDE curve
MIN_SAMPLES = 10


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot angle distribution KDE by residue size category."
    )
    p.add_argument("--input",  required=True, help="Path to aggregated angles.csv.")
    p.add_argument("--output", required=True, help="Output PNG path (e.g. results/plot.png).")
    p.add_argument("--xmin",   type=float, default=-180, help="KDE x-axis minimum (degrees).")
    p.add_argument("--xmax",   type=float, default=180,  help="KDE x-axis maximum (degrees).")
    return p.parse_args()


def _empty_plot(output: str, message: str) -> None:
    """Save a blank placeholder plot with a message when no data exist."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5, message, ha="center", va="center",
            fontsize=14, color="gray", transform=ax.transAxes)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)
    logger.warning("Saved empty placeholder plot → %s", output)


def main() -> None:
    args = parse_args()

    # ── Load data ─────────────────────────────────────────────────────────
    if not os.path.isfile(args.input):
        sys.exit(f"[ERROR] Input file not found: {args.input}")

    try:
        df = pd.read_csv(args.input)
    except Exception as exc:
        sys.exit(f"[ERROR] Cannot read CSV: {exc}")

    if df.empty or "angle" not in df.columns or "category" not in df.columns:
        _empty_plot(args.output, "No data available.")
        return

    df = df.dropna(subset=["angle", "category"])
    logger.info("Loaded %d records across %d categories.", len(df), df["category"].nunique())

    # ── Plot ──────────────────────────────────────────────────────────────
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.25)
    fig, axes = plt.subplots(
        2, 1,
        figsize=(13, 10),
        gridspec_kw={"height_ratios": [3, 1]},
    )

    ax_kde  = axes[0]
    ax_box  = axes[1]

    present_cats = [c for c in CATEGORY_ORDER if c in df["category"].values]

    # ── KDE panel ─────────────────────────────────────────────────────────
    plotted_any = False
    for cat in present_cats:
        subset = df.loc[df["category"] == cat, "angle"].dropna()
        if len(subset) < MIN_SAMPLES:
            logger.warning("'%s' has only %d samples – skipping KDE.", cat, len(subset))
            continue

        sns.kdeplot(
            data      = subset,
            ax        = ax_kde,
            label     = f"{cat}  (n = {len(subset):,})",
            color     = PALETTE.get(cat, "gray"),
            linestyle = LINE_STYLES.get(cat, "solid"),
            linewidth = 2.2,
            fill      = True,
            alpha     = 0.12,
            bw_adjust = 0.8,
        )
        plotted_any = True

    if not plotted_any:
        _empty_plot(args.output, "Insufficient data for KDE (< 10 samples per category).")
        return

    # Reference line at 0°
    ax_kde.axvline(0, color="black", linestyle="--", linewidth=1.0, alpha=0.45, label="0°")

    ax_kde.set_xlim(args.xmin, args.xmax)
    ax_kde.set_xlabel("")
    ax_kde.set_ylabel("Density", fontsize=13)
    ax_kde.set_title(
        "Signed CA→Centroid Angle Distribution\n"
        "X – Arg(helix) – X Tripeptides, Grouped by Left-Residue Size Category",
        fontsize=14, fontweight="bold", pad=12,
    )
    ax_kde.xaxis.set_major_locator(ticker.MultipleLocator(45))
    ax_kde.xaxis.set_minor_locator(ticker.MultipleLocator(15))
    ax_kde.legend(
        title="Left Residue Category",
        title_fontsize=11,
        fontsize=10,
        framealpha=0.85,
        loc="upper right",
    )
    sns.despine(ax=ax_kde)

    # ── Box / strip panel ─────────────────────────────────────────────────
    plot_df = df[df["category"].isin(present_cats)].copy()
    plot_df["category"] = pd.Categorical(plot_df["category"], categories=present_cats, ordered=True)

    sns.boxplot(
        data      = plot_df,
        x         = "category",
        y         = "angle",
        order     = present_cats,
        palette   = PALETTE,
        ax        = ax_box,
        linewidth = 1.2,
        fliersize = 1.5,
        flierprops= dict(alpha=0.3),
        width     = 0.55,
    )
    ax_box.axhline(0, color="black", linestyle="--", linewidth=0.9, alpha=0.45)
    ax_box.set_ylim(args.xmin, args.xmax)
    ax_box.set_xlabel("Left Residue Size Category", fontsize=12)
    ax_box.set_ylabel("Angle (°)", fontsize=12)
    ax_box.yaxis.set_major_locator(ticker.MultipleLocator(90))
    sns.despine(ax=ax_box)

    # ── Annotate sample counts below each box ─────────────────────────────
    for idx, cat in enumerate(present_cats):
        n = int((plot_df["category"] == cat).sum())
        ax_box.text(
            idx, args.xmin + 5,
            f"n={n:,}",
            ha="center", va="bottom",
            fontsize=8.5, color="dimgray",
        )

    # ── Finalise ──────────────────────────────────────────────────────────
    fig.tight_layout(h_pad=0.4)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    fig.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logger.info("Plot saved → %s", args.output)

    # ── Print descriptive stats to stdout ─────────────────────────────────
    stats = (
        plot_df.groupby("category", observed=True)["angle"]
        .agg(n="count", mean="mean", median="median", std="std",
             q25=lambda x: x.quantile(0.25),
             q75=lambda x: x.quantile(0.75))
        .round(2)
    )
    print("\n── Angle statistics by category ──")
    print(stats.to_string())
    print()


if __name__ == "__main__":
    main()
