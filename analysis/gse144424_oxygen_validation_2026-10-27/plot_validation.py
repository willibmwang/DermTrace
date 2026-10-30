"""Plot saved donor-level GSE144424 marker results without refitting the analysis."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator


OUTPUT = Path(__file__).resolve().parent
PANELS = (
    ("B-A", "B  ·  Six hours of hypoxia", "1% oxygen for 6 h − matched 10% reference"),
    ("C-A", "C  ·  Six hours after reoxygenation", "6 h at 1%, then 6 h at 10% − earlier reference"),
    ("D-A", "D  ·  Twenty-four hours after reoxygenation", "6 h at 1%, then 24 h at 10% − earlier reference"),
)


def main() -> None:
    plan = json.loads((OUTPUT / "analysis_plan.json").read_text())
    genes = plan["markers"]
    deltas = pd.read_csv(OUTPUT / "donor_deltas.csv", dtype={"donor_id": str})
    summary = pd.read_csv(OUTPUT / "marker_summary.csv")
    primary = deltas.loc[
        deltas["analysis"].eq("primary") & deltas["normalization"].eq("median_ratio")
    ].copy()
    forest = summary.loc[
        summary["analysis"].eq("primary")
        & summary["normalization"].eq("median_ratio")
        & summary["contrast"].eq("B-A")
    ].set_index("gene")
    if forest.index.has_duplicates or set(forest.index) != set(genes):
        raise ValueError("Expected one primary B-A result for every specified marker.")
    forest = forest.loc[genes]
    donors = sorted(primary["donor_id"].unique())
    if len(donors) != 15 or not forest["n_donors"].eq(15).all():
        raise ValueError("The primary plot requires the planned 15 complete donors.")
    if not primary["n_experiments"].eq(1).all():
        raise ValueError("Primary results must contain one matched experiment per donor.")

    matrices = {}
    for contrast, _, _ in PANELS:
        rows = primary.loc[primary["contrast"].eq(contrast)]
        if rows.duplicated(["donor_id", "gene"]).any():
            raise ValueError(f"Duplicate donor-marker result in {contrast}.")
        matrix = rows.pivot(index="donor_id", columns="gene", values="delta")
        if set(matrix.columns) != set(genes) or set(matrix.index) != set(donors):
            raise ValueError(f"Incomplete marker or donor axis in {contrast}.")
        matrix = matrix.loc[donors, genes]
        if not np.isfinite(matrix.to_numpy()).all():
            raise ValueError(f"Nonfinite donor differences in {contrast}.")
        matrices[contrast] = matrix
    values = forest[["median_delta", "ci_low", "ci_high", "q_bh"]].to_numpy()
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite primary estimates, intervals or adjusted tests.")
    if not (forest["ci_low"] <= forest["ci_high"]).all():
        raise ValueError("Inverted bootstrap interval.")

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.labelcolor": "#273440", "text.color": "#273440",
        "xtick.color": "#46515b", "ytick.color": "#46515b",
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(16.5, 13.5), facecolor="white")
    grid = fig.add_gridspec(
        2, 3, width_ratios=[1, 1, 0.035], height_ratios=[1, 1],
        left=0.075, right=0.94, bottom=0.16, top=0.83, hspace=0.45, wspace=0.31,
    )
    ax = fig.add_subplot(grid[0, 0])
    supported = int((forest["median_delta"].gt(0) & forest["q_bh"].le(0.05)).sum())
    fig.suptitle(
        f"{supported} of 8 markers support the expected hypoxia direction across donors",
        x=0.075, y=0.97, ha="left", fontsize=20, fontweight="bold",
    )
    fig.text(
        0.075, 0.935,
        "GSE144424 · 15 human donors · iPSC-derived cardiomyocytes · exploratory external validation",
        ha="left", fontsize=12, color="#586875",
    )
    fig.text(
        0.075, 0.91,
        "Primary selection: one complete A–D experiment per donor; repeated differentiations do not add donors.",
        ha="left", fontsize=11, color="#586875",
    )

    ax.set_title("A  ·  Hypoxia versus matched reference", loc="left", pad=32,
                 fontsize=14, fontweight="bold")
    ax.text(0, 1.035, "Median donor difference and unadjusted 95% bootstrap interval",
            transform=ax.transAxes, fontsize=10.5, color="#586875")
    lower = min(0.0, float(forest["ci_low"].min()))
    upper = max(0.0, float(forest["ci_high"].max()))
    span = max(upper - lower, 1.0)
    ax.set_xlim(lower - 0.10 * span, upper + 0.08 * span)
    ax.set_ylim(len(genes) - 0.4, -0.85)
    ax.axvline(0, color="#6d7b85", linestyle=(0, (3, 3)), linewidth=1.2, zorder=1)
    for y, gene in enumerate(genes):
        row = forest.loc[gene]
        if y % 2 == 0:
            ax.axhspan(y - 0.45, y + 0.45, color="#f3f6f8", zorder=0)
        color = "#a92b38" if row["median_delta"] > 0 else "#27638d"
        ax.plot([row["ci_low"], row["ci_high"]], [y, y], color=color, linewidth=2.2)
        ax.scatter(row["median_delta"], y, s=65, color=color, edgecolors="white", zorder=3)
        ax.text(1.025, y, f"{int(row['n_positive'])}/15   q={row['q_bh']:.3g}",
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=10)
    ax.text(1.025, -0.74, "Positive donors · BH q", transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=9.5, color="#586875")
    ax.set_yticks(np.arange(len(genes)), genes)
    ax.set_xlabel("Paired difference in log₂(normalized count + 1)", labelpad=12)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#b9c3ca")
    ax.tick_params(axis="y", length=0, pad=10)
    ax.grid(axis="x", color="#edf0f2", zorder=0)
    # Leave a right-hand column for the donor support annotations.
    position = ax.get_position()
    ax.set_position([position.x0, position.y0, position.width * 0.74, position.height])

    vmax = max(0.5, np.ceil(max(np.abs(m.to_numpy()).max() for m in matrices.values()) * 2) / 2)
    norm = Normalize(vmin=-vmax, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")
    heat_axes = [fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 0]),
                 fig.add_subplot(grid[1, 1])]
    for heat, (contrast, title, subtitle) in zip(heat_axes, PANELS):
        matrix = matrices[contrast]
        im = heat.imshow(matrix.to_numpy(), cmap=cmap, norm=norm, aspect="auto")
        heat.set_title(title, loc="left", pad=32, fontsize=14, fontweight="bold")
        heat.text(0, 1.035, subtitle, transform=heat.transAxes, fontsize=10.5, color="#586875")
        heat.set_xticks(np.arange(len(genes)), genes, rotation=35, ha="right")
        heat.set_yticks(np.arange(len(donors)), donors)
        heat.set_ylabel("Donor", labelpad=8)
        heat.tick_params(axis="both", length=0, pad=6, labelsize=10)
        heat.set_xticks(np.arange(-0.5, len(genes), 1), minor=True)
        heat.set_yticks(np.arange(-0.5, len(donors), 1), minor=True)
        heat.grid(which="minor", color="white", linewidth=0.7)
        heat.tick_params(which="minor", bottom=False, left=False)
        for y, x in np.ndindex(matrix.shape):
            delta = matrix.iat[y, x]
            color = "white" if abs(delta) > 0.55 * vmax else "#25313b"
            heat.text(x, y, f"{delta:+.1f}", ha="center", va="center", color=color, fontsize=9.5)
        for spine in heat.spines.values():
            spine.set_visible(False)

    color_ax = fig.add_subplot(grid[:, 2])
    colorbar = fig.colorbar(im, cax=color_ax)
    colorbar.set_label("Donor-paired difference in log₂(normalized count + 1)", labelpad=14)
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(length=0, pad=6)
    fig.text(
        0.075, 0.105,
        "Heatmaps use the same donor order and symmetric scale: red = higher, blue = lower. "
        "Labels are rounded; signed zeros can represent small nonzero differences.",
        ha="left", fontsize=10.5,
    )
    fig.text(
        0.075, 0.078,
        "Panel A: one-sided exact sign tests; Benjamini–Hochberg q across eight planned markers. "
        "Intervals: 20,000 donor resamples, conditional on normalization.",
        ha="left", fontsize=10.5,
    )
    fig.text(
        0.075, 0.051,
        "Reoxygenation is compared with the earlier baseline; later time-matched reference controls are absent. "
        "Expression proximity does not establish functional recovery.",
        ha="left", fontsize=10.5,
    )
    fig.text(
        0.075, 0.024,
        "Median-of-ratios normalization within this dataset only. These differences are not on the Tsao TPM scale. "
        "Frozen DermTrace datasets, models and benchmark are unchanged.",
        ha="left", fontsize=10.5, color="#586875",
    )
    for extension in ("png", "pdf"):
        path = OUTPUT / f"donor_oxygen_response.{extension}"
        metadata = {"CreationDate": None, "ModDate": None} if extension == "pdf" else None
        fig.savefig(path, dpi=200, facecolor="white", bbox_inches="tight", metadata=metadata)
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
