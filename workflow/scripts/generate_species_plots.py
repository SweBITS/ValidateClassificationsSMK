#!/usr/bin/env python3
"""
generate_species_plots.py — Stage 5: Publication-Grade 3-Panel Diagnostic Figure Generator
==========================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
Summary correlation numbers alone cannot reveal *why* a putative eDNA detection succeeded
or failed validation. A visual diagnostic figure for each target species allows researchers
and manuscript reviewers to immediately inspect three complementary lines of evidence:

  - **Panel A (Genomic Uniformity & Outlier Regression Scatter)**:
    Plots Mapped eDNA Read Pairs ($Y$) against Bin Classifiable Region Length ($X$) across
    all 100 variable-size bins ($0.3\\bar{L}$ to $1.7\\bar{L}$). Retained bins (`●`) are
    distinguished from pruned outlier bins (`▲`, annotated with removal step `#1`, `#2`, ...),
    alongside both the Initial Unpruned OLS Regression Line (dashed grey) and the Final
    Validated OLS Regression Line (solid color with an exact $t$-distribution 95% confidence
    interval band).

  - **Panel B (Bidirectional Outlier Pruning Trajectory)**:
    Dual-axis trajectory showing how Pearson correlation $R$ (left axis) increases and
    Retained CR Sequence % (right axis) decreases at each iterative outlier removal step,
    relative to the $R = 0.90$ validation target and the $50\\%$ CR retention floor.

  - **Panel C (10-Year / 516-Week Longitudinal Time Series Comparison, 2013–2022)**:
    Compares weekly Raw Kraken 2 Clade Reads against All-CR Mapped Reads and Validated
    Clean-Bin Reads across the 516-week Ljungbyhed monitoring series, revealing whether
    seasonal ecological phenology is preserved or if sporadic false-positive spikes are
    eliminated after removing artifact bins.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import ValidationStatus, setup_logger


STATUS_PALETTE: Dict[str, Dict[str, str]] = {
    ValidationStatus.VALIDATED_CLEAN: {
        "primary": "#0f766e",   # Deep teal
        "fill": "#99f6e4",
        "badge_bg": "#d1fae5",
        "badge_fg": "#065f46",
    },
    ValidationStatus.VALIDATED_PRUNED: {
        "primary": "#0284c7",   # Ocean blue
        "fill": "#bae6fd",
        "badge_bg": "#e0f2fe",
        "badge_fg": "#075985",
    },
    ValidationStatus.FAILED_LOW_CORRELATION: {
        "primary": "#d97706",   # Amber
        "fill": "#fde68a",
        "badge_bg": "#fef3c7",
        "badge_fg": "#92400e",
    },
    ValidationStatus.FAILED_LOW_OCCUPANCY: {
        "primary": "#dc2626",   # Crimson red
        "fill": "#fecaca",
        "badge_bg": "#fee2e2",
        "badge_fg": "#991b1b",
    },
}


def compute_ols_ci_band(
    x_obs: np.ndarray,
    y_obs: np.ndarray,
    x_grid: np.ndarray,
    confidence: float = 0.95,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute OLS fitted line `y_hat` and exact two-sided $t$-distribution confidence interval
    `(ci_low, ci_high)` across `x_grid` for retained bins `(x_obs, y_obs)`.
    """
    n = len(x_obs)
    if n < 3 or np.std(x_obs) == 0.0:
        mean_y = float(np.mean(y_obs)) if n > 0 else 0.0
        flat = np.full_like(x_grid, mean_y, dtype=np.float64)
        return flat, flat, flat

    lr = stats.linregress(x_obs, y_obs)
    slope = float(lr.slope) if not math.isnan(lr.slope) else 0.0
    intercept = float(lr.intercept) if not math.isnan(lr.intercept) else 0.0
    y_hat_grid = intercept + slope * x_grid

    y_hat_obs = intercept + slope * x_obs
    dof = max(1, n - 2)
    s_err = math.sqrt(float(np.sum((y_obs - y_hat_obs) ** 2)) / dof)
    x_mean = float(np.mean(x_obs))
    ssx = float(np.sum((x_obs - x_mean) ** 2))
    if ssx <= 0.0 or s_err <= 0.0:
        return y_hat_grid, y_hat_grid, y_hat_grid

    t_val = float(stats.t.ppf((1.0 + confidence) / 2.0, dof))
    se_line = s_err * np.sqrt((1.0 / n) + ((x_grid - x_mean) ** 2) / ssx)
    margin = t_val * se_line
    return y_hat_grid, y_hat_grid - margin, y_hat_grid + margin


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 5: Generate 3-panel publication diagnostic plot (PNG + PDF) for a target species."
    )
    parser.add_argument("--taxid", type=int, required=True, help="Target species TaxID")
    parser.add_argument("--bin-regression-tsv", required=True, help="Path to {t_id}.bin_regression_data.tsv")
    parser.add_argument("--trajectory-tsv", required=True, help="Path to {t_id}.pruning_trajectory.tsv")
    parser.add_argument("--weekly-tsv", required=True, help="Path to {t_id}.weekly_timeseries_516w.tsv")
    parser.add_argument("--metrics-json", required=True, help="Path to {t_id}.validation_metrics.json")
    parser.add_argument("--target-r", type=float, default=0.90, help="Target Pearson R threshold")
    parser.add_argument("--min-retained-cr-frac", type=float, default=0.50, help="Minimum retained CR fraction")
    parser.add_argument("--out-png", required=True, help="Output PNG path (300 DPI)")
    parser.add_argument("--out-pdf", required=True, help="Output vector PDF path")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger(f"plot_species_{args.taxid}", args.log)

    # 1. Load metrics JSON and tabular inputs
    with open(args.metrics_json, "r", encoding="utf-8") as jfh:
        metrics: Dict[str, Any] = json.load(jfh)

    with open(args.bin_regression_tsv, "r", encoding="utf-8") as bfh:
        bin_rows = list(csv.DictReader(bfh, delimiter="\t"))

    with open(args.trajectory_tsv, "r", encoding="utf-8") as tfh:
        traj_rows = list(csv.DictReader(tfh, delimiter="\t"))

    with open(args.weekly_tsv, "r", encoding="utf-8") as wfh:
        weekly_rows = list(csv.DictReader(wfh, delimiter="\t"))

    status = metrics.get("validation_status", ValidationStatus.FAILED_LOW_CORRELATION)
    palette = STATUS_PALETTE.get(
        status,
        {"primary": "#475569", "fill": "#cbd5e1", "badge_bg": "#f1f5f9", "badge_fg": "#334155"},
    )

    # Set clean publication rcParams
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "axes.edgecolor": "#cbd5e1",
        "axes.linewidth": 0.9,
        "grid.color": "#e2e8f0",
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
    })

    fig = plt.figure(figsize=(15, 9.5), dpi=300, facecolor="white")
    gs = gridspec.GridSpec(
        2, 2,
        figure=fig,
        height_ratios=[1.15, 1.0],
        width_ratios=[1.15, 1.0],
        hspace=0.34,
        wspace=0.24,
        left=0.06,
        right=0.95,
        top=0.88,
        bottom=0.07,
    )

    ax_scatter = fig.add_subplot(gs[0, 0])
    ax_traj = fig.add_subplot(gs[0, 1])
    ax_ts = fig.add_subplot(gs[1, :])

    # --- Header Banner --- #
    sci_name = metrics.get("species_name", f"TaxID {args.taxid}")
    fig.suptitle(
        f"{sci_name} (TaxID: {args.taxid})  —  [{status}]",
        fontsize=15,
        fontweight="bold",
        color="#0f172a",
        y=0.96,
    )
    subtitle_str = (
        f"Total CR: {metrics.get('total_cr_bp', 0):,} bp ({metrics.get('retained_cr_pct', 0):.1f}% retained)   |   "
        f"Mapped Pairs: {metrics.get('total_mapped_read_pairs', 0):,} ({metrics.get('retained_mapped_reads_pct', 0):.1f}% in clean bins)   |   "
        f"Initial R = {metrics.get('initial_pearson_r', 0):.3f}  →  Final R = {metrics.get('final_pearson_r', 0):.3f} "
        f"(ρ = {metrics.get('final_spearman_rho', 0):.3f}, {metrics.get('bins_pruned_count', 0)} bins pruned)"
    )
    fig.text(0.5, 0.922, subtitle_str, ha="center", va="center", fontsize=10, color="#334155")

    # ======================================================================== #
    # PANEL A: Regression Scatter Plot (Retained vs. Pruned Bins + 95% CI)
    # ======================================================================== #
    x_all = np.array([float(r["actual_cr_bp"]) for r in bin_rows], dtype=np.float64)
    y_all = np.array([float(r["mapped_reads"]) for r in bin_rows], dtype=np.float64)
    is_ret = np.array([r["is_retained"].lower() == "true" for r in bin_rows], dtype=bool)
    pruned_steps = [int(r["pruned_step"]) for r in bin_rows]

    ax_scatter.grid(True, alpha=0.7, zorder=0)

    if len(x_all) > 0 and np.max(x_all) > 0:
        x_grid = np.linspace(float(np.min(x_all)) * 0.95, float(np.max(x_all)) * 1.03, 200)

        # Initial regression line across all 100 bins (dashed slate grey)
        init_slope = float(metrics.get("initial_ols_slope", 0.0))
        init_intercept = float(metrics.get("initial_ols_intercept", 0.0))
        y_init_grid = init_intercept + init_slope * x_grid
        ax_scatter.plot(
            x_grid,
            y_init_grid,
            linestyle="--",
            linewidth=1.5,
            color="#64748b",
            label=f"Initial OLS (all {len(x_all)} bins, R={metrics.get('initial_pearson_r', 0):.3f})",
            zorder=2,
        )

        # Final validated regression line + 95% CI band on retained bins
        x_ret = x_all[is_ret]
        y_ret = y_all[is_ret]
        if len(x_ret) >= 3:
            y_fin_grid, ci_low, ci_high = compute_ols_ci_band(x_ret, y_ret, x_grid)
            ax_scatter.fill_between(
                x_grid,
                ci_low,
                ci_high,
                color=palette["fill"],
                alpha=0.45,
                label="95% CI (retained bins)",
                zorder=1,
            )
            ax_scatter.plot(
                x_grid,
                y_fin_grid,
                linestyle="-",
                linewidth=2.2,
                color=palette["primary"],
                label=f"Final OLS ({len(x_ret)} retained bins, R={metrics.get('final_pearson_r', 0):.3f})",
                zorder=3,
            )

        # Plot retained bins (●)
        ax_scatter.scatter(
            x_ret,
            y_ret,
            s=42,
            c=palette["primary"],
            edgecolors="white",
            linewidths=0.6,
            alpha=0.88,
            marker="o",
            label=f"Retained bins (n={len(x_ret)})",
            zorder=4,
        )

        # Plot pruned outlier bins (▲) with step order labels
        x_pruned = x_all[~is_ret]
        y_pruned = y_all[~is_ret]
        if len(x_pruned) > 0:
            ax_scatter.scatter(
                x_pruned,
                y_pruned,
                s=68,
                c="#e11d48",
                edgecolors="#881337",
                linewidths=0.8,
                alpha=0.92,
                marker="^",
                label=f"Pruned outlier bins (n={len(x_pruned)})",
                zorder=5,
            )
            # Label up to top 15 pruned bins by step order so plot stays legible
            pruned_indices = [i for i, flag in enumerate(is_ret) if not flag]
            pruned_indices.sort(key=lambda idx: pruned_steps[idx])
            for idx in pruned_indices[:15]:
                step_num = pruned_steps[idx]
                ax_scatter.annotate(
                    f"#{step_num}",
                    xy=(x_all[idx], y_all[idx]),
                    xytext=(5, 4),
                    textcoords="offset points",
                    fontsize=7.5,
                    fontweight="bold",
                    color="#9f1239",
                    zorder=6,
                )

    ax_scatter.set_title("A. Bin Classifiable Region Length vs. Mapped eDNA Reads", fontsize=11, fontweight="bold", color="#0f172a", loc="left", pad=8)
    ax_scatter.set_xlabel("Bin Classifiable Region Length (bp)", fontsize=9.5, color="#1e293b")
    ax_scatter.set_ylabel("Mapped eDNA Read Pairs", fontsize=9.5, color="#1e293b")
    ax_scatter.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax_scatter.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax_scatter.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.9, fontsize=8)

    # ======================================================================== #
    # PANEL B: Dual-Trace Outlier Pruning Trajectory Curve
    # ======================================================================== #
    iters = [int(r["iteration"]) for r in traj_rows]
    r_trace = [float(r["pearson_r"]) for r in traj_rows]
    cr_pct_trace = [float(r["retained_cr_pct"]) for r in traj_rows]

    ax_traj.grid(True, alpha=0.7, zorder=0)
    ax_traj_right = ax_traj.twinx()

    line_r = ax_traj.plot(
        iters,
        r_trace,
        color="#2563eb",
        marker="o",
        markersize=5,
        linewidth=2.0,
        label="Pearson R",
        zorder=4,
    )
    thresh_r = ax_traj.axhline(
        args.target_r,
        color="#16a34a",
        linestyle="--",
        linewidth=1.3,
        label=f"Target R = {args.target_r:.2f}",
        zorder=2,
    )

    line_cr = ax_traj_right.plot(
        iters,
        cr_pct_trace,
        color="#d97706",
        marker="s",
        markersize=4.5,
        linewidth=1.8,
        linestyle="-",
        label="Retained CR (%)",
        zorder=3,
    )
    thresh_cr = ax_traj_right.axhline(
        100.0 * args.min_retained_cr_frac,
        color="#dc2626",
        linestyle=":",
        linewidth=1.3,
        label=f"Min CR Retention ({100.0 * args.min_retained_cr_frac:.0f}%)",
        zorder=2,
    )

    ax_traj.set_title("B. Recursive Outlier Pruning Trajectory", fontsize=11, fontweight="bold", color="#0f172a", loc="left", pad=8)
    ax_traj.set_xlabel("Outlier Bins Pruned (Iteration)", fontsize=9.5, color="#1e293b")
    ax_traj.set_ylabel("Pearson Correlation (R)", fontsize=9.5, color="#2563eb")
    ax_traj_right.set_ylabel("Retained CR Sequence (%)", fontsize=9.5, color="#d97706")

    min_r_plot = min(0.0, min(r_trace) - 0.05) if r_trace else 0.0
    ax_traj.set_ylim(min_r_plot, 1.03)
    ax_traj_right.set_ylim(max(0.0, min(cr_pct_trace) - 10.0) if cr_pct_trace else 0.0, 105.0)
    ax_traj.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    # Combine legends from both axes
    lines_all = line_r + [thresh_r] + line_cr + [thresh_cr]
    labels_all = [l.get_label() for l in lines_all]
    ax_traj.legend(lines_all, labels_all, loc="lower left", frameon=True, facecolor="white", framealpha=0.9, fontsize=8)

    # ======================================================================== #
    # PANEL C: 516-Week Time Series Comparison (2013–2022)
    # ======================================================================== #
    n_weeks = len(weekly_rows)
    w_idx = np.arange(n_weeks)
    raw_reads = np.array([int(r["raw_kraken_reads"]) for r in weekly_rows], dtype=np.float64)
    all_cr_reads = np.array([int(r["all_cr_mapped_reads"]) for r in weekly_rows], dtype=np.float64)
    clean_reads = np.array([int(r["validated_clean_bin_reads"]) for r in weekly_rows], dtype=np.float64)
    years = [int(r["year"]) for r in weekly_rows]

    ax_ts.grid(True, alpha=0.6, zorder=0)

    # Alternating annual background bands and year tick locations
    year_ticks: List[int] = []
    year_labels: List[str] = []
    if years:
        unique_years = sorted(set(years))
        for y_i, yr in enumerate(unique_years):
            yr_positions = [i for i, val in enumerate(years) if val == yr]
            if yr_positions:
                start_pos, end_pos = yr_positions[0], yr_positions[-1]
                year_ticks.append(start_pos)
                year_labels.append(str(yr))
                if y_i % 2 == 1:
                    ax_ts.axvspan(start_pos - 0.5, end_pos + 0.5, color="#f8fafc", zorder=0)

    ax_ts.fill_between(
        w_idx,
        0,
        raw_reads,
        color="#94a3b8",
        alpha=0.30,
        label=f"Raw Kraken 2 Clade Reads (sum = {int(np.sum(raw_reads)):,})",
        zorder=1,
    )
    ax_ts.plot(
        w_idx,
        raw_reads,
        color="#64748b",
        linewidth=0.8,
        alpha=0.75,
        zorder=2,
    )
    ax_ts.plot(
        w_idx,
        all_cr_reads,
        color="#38bdf8",
        linewidth=1.0,
        alpha=0.85,
        label=f"All-CR Mapped Reads (sum = {int(np.sum(all_cr_reads)):,})",
        zorder=3,
    )
    ax_ts.plot(
        w_idx,
        clean_reads,
        color=palette["primary"],
        linewidth=1.5,
        label=f"Validated Clean-Bin Reads (sum = {int(np.sum(clean_reads)):,})",
        zorder=4,
    )

    ax_ts.set_title(
        "C. 10-Year Weekly eDNA Time Series (516 Weeks, 2013–2022): Raw Kraken 2 vs. Validated Clean-Bin Signal",
        fontsize=11,
        fontweight="bold",
        color="#0f172a",
        loc="left",
        pad=8,
    )
    ax_ts.set_xlabel("Monitoring Year (Weekly Air-Filter Samples, Ljungbyhed)", fontsize=9.5, color="#1e293b")
    ax_ts.set_ylabel("Weekly Read Pairs", fontsize=9.5, color="#1e293b")
    if year_ticks:
        ax_ts.set_xticks(year_ticks)
        ax_ts.set_xticklabels(year_labels, fontsize=9)
    ax_ts.set_xlim(0, max(1, n_weeks - 1))
    ax_ts.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax_ts.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, fontsize=8.5)

    # Save PNG and PDF
    out_png = Path(args.out_png)
    out_pdf = Path(args.out_pdf)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    logger.info(f"Saved 3-panel diagnostic figure for TaxID {args.taxid} to {out_png.name} and {out_pdf.name}")


if __name__ == "__main__":
    main()
