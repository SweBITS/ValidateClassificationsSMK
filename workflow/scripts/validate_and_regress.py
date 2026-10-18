#!/usr/bin/env python3
"""
validate_and_regress.py — Stage 5: Occupancy Gate, Bidirectional OLS Pruning & 516-Week Reconstruction
======================================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
Even when a species is genuinely present in an environmental sample, a small subset of its
100 Classifiable Region (CR) bins can deviate sharply from the linear genome-wide relationship
between Mapped eDNA Reads ($Y$) and Bin CR Length ($X$):
  - **Positive Outliers (Above the Regression Line, $Y_i \\gg \\hat{Y}_i$)**:
    Bins containing high-copy ribosomal RNA arrays, transposable elements, mitochondrial/plastid
    insertions (NUMTs/NUPTs), or ultra-conserved house-keeping exons that attract cross-mapping
    reads from abundant co-occurring organisms.
  - **Negative Outliers (Below the Regression Line, $Y_i \\ll \\hat{Y}_i$)**:
    Bins containing misassembled contaminant scaffolds in the NCBI reference assembly or
    hyper-variable strain-specific islands absent from the local population in Ljungbyhed.

Conversely, a **false-positive taxon** (absent from the environment, detected solely via
cross-mapping to a few conserved loci) exhibits empty CRs across most of its genome with all
reads piled into 1–10 localized spikes.

HOW OUR TWO-STEP VALIDATION ENGINE SEPARATES TRUE POSITIVES FROM ARTIFACTS:
---------------------------------------------------------------------------
1. **Pre-Regression Sanity Gate (Minimum Bin Occupancy & Positive Initial Slope)**:
   - Before any outlier pruning is allowed, at least `min_occupied_bins_frac` (default 50%,
     i.e. >= 50 of 100 bins) must contain >= 1 mapped read pair, AND the initial OLS slope
     $\\beta_1$ must be strictly positive ($\\beta_1 > 0$).
   - **Why this gate is essential**: Without an upfront occupancy requirement, a false-positive
     species with 95 empty bins ($Y=0$) and 5 artifact spikes could iteratively prune away
     bins or fit a flat/negative line. Failing this gate immediately classifies the species as
     `FAILED_LOW_OCCUPANCY`.

2. **Bidirectional Recursive OLS Outlier Pruning ($R \\ge 0.90$, Retained CR $\\ge 50\\%$)**:
   - We fit an Ordinary Least Squares (OLS) line $\\hat{Y} = \\beta_0 + \\beta_1 X$ across all
     100 bins and compute Pearson $R$ (and Spearman $\\rho$).
   - If $R \\ge 0.90$ and $\\beta_1 > 0$ initially, no pruning is needed (`VALIDATED_CLEAN`).
   - Otherwise, at each iteration we identify the single bin $i^*$ with the **largest absolute
     vertical residual** $|Y_i - \\hat{Y}_i|$ from the current regression line—removing both
     spikes above the line and contaminant/dropout bins below the line.
   - If removing bin $i^*$ would drop retained CR sequence below `min_retained_cr_frac * L_CR`
     (default 50% of total CR length), or if $\\beta_1 \\le 0$, pruning terminates and the
     species is flagged as `FAILED_LOW_CORRELATION`.
   - If Pearson $R$ reaches $\\ge 0.90$ with $\\beta_1 > 0$ while retaining $\\ge 50\\%$ of total
     CR length, pruning stops immediately and the species is classified as `VALIDATED_PRUNED`.

3. **Cleaned 516-Week Time Series Reconstruction**:
   - Each mapped read pair (`read_id`) is deduplicated at the fragment level (taking the highest
     MAPQ / longest aligned mate so 1 Kraken-classified read pair contributes at most 1 mapped
     fragment) and decomposed by its embedded `|sample=<sample_id>` tag across all 516 weeks:
       * `raw_kraken_reads`: Baseline weekly clade read pairs from Kraken 2.
       * `all_cr_mapped_reads`: Weekly read pairs mapping to all 100 CR bins.
       * `validated_clean_bin_reads`: Weekly read pairs mapping **strictly to the retained
         (non-pruned) CR bins** (set to 0 in a separate `strictly_validated_signal` column if
         the species failed validation, while preserving `retained_bin_reads` for inspection).
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import polars as pl
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import ValidationStatus, setup_logger


def fit_ols_and_correlations(
    x: np.ndarray,
    y: np.ndarray,
) -> Dict[str, float]:
    """
    Compute Ordinary Least Squares (OLS) slope, intercept, fitted values, Pearson R (with p-value),
    and Spearman rho (with p-value) for arrays `x` (Bin CR length) and `y` (Mapped reads).
    """
    n = len(x)
    if n < 3 or np.std(x) == 0.0:
        return {
            "slope": 0.0,
            "intercept": float(np.mean(y)) if n > 0 else 0.0,
            "pearson_r": 0.0,
            "pearson_p": 1.0,
            "spearman_rho": 0.0,
            "spearman_p": 1.0,
        }

    # Check if y is constant (e.g. all zeros)
    if np.std(y) == 0.0:
        return {
            "slope": 0.0,
            "intercept": float(y[0]),
            "pearson_r": 0.0,
            "pearson_p": 1.0,
            "spearman_rho": 0.0,
            "spearman_p": 1.0,
        }

    lr = stats.linregress(x, y)
    sp_rho, sp_p = stats.spearmanr(x, y)

    slope = float(lr.slope) if not math.isnan(lr.slope) else 0.0
    intercept = float(lr.intercept) if not math.isnan(lr.intercept) else 0.0
    r_val = float(lr.rvalue) if not math.isnan(lr.rvalue) else 0.0
    p_val = float(lr.pvalue) if not math.isnan(lr.pvalue) else 1.0
    rho_val = float(sp_rho) if not math.isnan(sp_rho) else 0.0
    rho_p = float(sp_p) if not math.isnan(sp_p) else 1.0

    return {
        "slope": slope,
        "intercept": intercept,
        "pearson_r": r_val,
        "pearson_p": p_val,
        "spearman_rho": rho_val,
        "spearman_p": rho_p,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 5: Run occupancy gate, bidirectional recursive OLS outlier pruning, and 516-week timeseries reconstruction."
    )
    parser.add_argument("--taxid", type=int, required=True, help="Target species TaxID")
    parser.add_argument("--species-name", required=True, help="Scientific name of species")
    parser.add_argument("--species-slug", required=True, help="Sanitized species slug")
    parser.add_argument("--bins-tsv", required=True, help="Path to Stage 3 {t_id}.bins_metadata.tsv")
    parser.add_argument("--alignments-parquet", required=True, help="Path to Stage 4 {t_id}.passing_alignments.parquet")
    parser.add_argument("--mapping-stats-tsv", required=True, help="Path to Stage 4 {t_id}.mapping_stats.tsv")
    parser.add_argument("--raw-weekly-tsv", required=True, help="Path to Stage 0 {t_id}.raw_weekly.tsv")
    parser.add_argument("--out-bin-regression-tsv", required=True, help="Output TSV: {t_id}.bin_regression_data.tsv")
    parser.add_argument("--out-trajectory-tsv", required=True, help="Output TSV: {t_id}.pruning_trajectory.tsv")
    parser.add_argument("--out-weekly-tsv", required=True, help="Output TSV: {t_id}.weekly_timeseries_516w.tsv")
    parser.add_argument("--out-metrics-json", required=True, help="Output JSON: {t_id}.validation_metrics.json")
    parser.add_argument("--min-occupied-frac", type=float, default=0.50, help="Minimum fraction of bins with >= 1 read (default: 0.50)")
    parser.add_argument("--target-r", type=float, default=0.90, help="Target Pearson R threshold (default: 0.90)")
    parser.add_argument("--min-retained-cr-frac", type=float, default=0.50, help="Minimum fraction of total CR bp retained (default: 0.50)")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger(f"validate_regress_{args.taxid}", args.log)

    # 1. Load Stage 3 bins metadata
    bins_df = pl.read_csv(args.bins_tsv, separator="\t")
    n_bins = len(bins_df)
    total_cr_bp = int(bins_df["actual_cr_bp"].sum()) if n_bins > 0 else 0

    # 2. Load Stage 4 alignments & collapse each strict 2-mate proper pair on the same CR to 1 fragment count
    #    Why:
    #    Kraken 2 classifies paired-end fragments (1 count per R1+R2 pair). We require BOTH mates
    #    (R1 and R2) of a read pair to map as a proper pair within the exact same CR (`cr_id`) and
    #    both to pass all Stage 4 alignment filters. Because `filter_sam_to_parquet.py` (and `minimap2`)
    #    emits the two mates of each passing pair as adjacent rows `(row i, row i+1)` with the same
    #    `read_id`, same `cr_id`, and complementary `mate` (`1` vs `2`), we can extract one row per
    #    complete same-CR proper pair in O(N) time and low memory (~7.4 GB RSS even on 52M rows).
    aln_df = pl.read_parquet(
        args.alignments_parquet,
        columns=["read_id", "year_week", "mate", "bin_id", "cr_id"],
    )
    if len(aln_df) > 0:
        dedup_pairs_df = aln_df.filter(
            (pl.col("read_id") == pl.col("read_id").shift(1))
            & (pl.col("cr_id") == pl.col("cr_id").shift(1))
            & (pl.col("mate") != pl.col("mate").shift(1))
        )
    else:
        dedup_pairs_df = aln_df

    total_mapped_pairs = len(dedup_pairs_df)

    # Count mapped read pairs per bin_id
    if total_mapped_pairs > 0:
        bin_counts_df = dedup_pairs_df.group_by("bin_id").agg(pl.len().alias("mapped_reads"))
        bin_count_map: Dict[str, int] = dict(
            zip(bin_counts_df["bin_id"].to_list(), bin_counts_df["mapped_reads"].to_list())
        )
    else:
        bin_count_map = {}

    bin_ids: List[str] = bins_df["bin_id"].to_list()
    bin_indices: List[int] = [int(x) for x in bins_df["bin_index"].to_list()]
    x_all = np.array([float(x) for x in bins_df["actual_cr_bp"].to_list()], dtype=np.float64)
    y_all = np.array([float(bin_count_map.get(bid, 0)) for bid in bin_ids], dtype=np.float64)

    occupied_bins = int(np.sum(y_all >= 1.0))
    occupied_frac = (occupied_bins / float(n_bins)) if n_bins > 0 else 0.0

    # Compute initial OLS fit across all 100 bins
    init_fit = fit_ols_and_correlations(x_all, y_all)
    init_r = init_fit["pearson_r"]
    init_p = init_fit["pearson_p"]
    init_rho = init_fit["spearman_rho"]
    init_rho_p = init_fit["spearman_p"]
    init_slope = init_fit["slope"]
    init_intercept = init_fit["intercept"]

    # Track pruning state
    retained_mask = np.ones(n_bins, dtype=bool)
    pruned_order: Dict[int, int] = {}         # array_idx -> iteration_removed (1-based)
    pruned_residual: Dict[int, float] = {}    # array_idx -> vertical residual (Y - Y_hat) at removal
    trajectory_rows: List[Dict[str, Any]] = []

    # Record Iteration 0 (Initial unpruned state)
    trajectory_rows.append({
        "iteration": 0,
        "bin_id_removed": "NONE",
        "bin_index_removed": 0,
        "cr_bp_removed": 0,
        "reads_removed": 0,
        "vertical_residual": 0.0,
        "retained_bins_count": n_bins,
        "retained_cr_bp": total_cr_bp,
        "retained_cr_pct": 100.0 if total_cr_bp > 0 else 0.0,
        "retained_reads": total_mapped_pairs,
        "retained_reads_pct": 100.0 if total_mapped_pairs > 0 else 0.0,
        "pearson_r": round(init_r, 6),
        "pearson_p": init_p,
        "spearman_rho": round(init_rho, 6),
        "ols_slope": init_slope,
        "ols_intercept": init_intercept,
    })

    # 3. Evaluate Stage 5 Gates & Bidirectional Recursive Linear Regression
    curr_fit = dict(init_fit)
    curr_retained_bp = total_cr_bp
    curr_retained_reads = total_mapped_pairs

    if total_cr_bp < n_bins or n_bins == 0:
        status = ValidationStatus.EXCLUDED_INSUFFICIENT_CR_LENGTH
        reason = f"Total Classifiable Region sequence ({total_cr_bp:,} bp) is insufficient to form {n_bins} bins"
    elif occupied_frac + 1e-9 < args.min_occupied_frac:
        status = ValidationStatus.FAILED_LOW_OCCUPANCY
        reason = (
            f"Bin occupancy ({occupied_bins}/{n_bins} bins = {100.0 * occupied_frac:.1f}%) "
            f"is below min_occupied_bins_frac ({100.0 * args.min_occupied_frac:.1f}%)"
        )
    elif init_slope <= 0.0:
        status = ValidationStatus.FAILED_LOW_OCCUPANCY
        reason = (
            f"Initial regression slope across {n_bins} bins is non-positive "
            f"(slope = {init_slope:.4e}, R = {init_r:.4f})"
        )
    elif init_r >= args.target_r and init_slope > 0.0:
        status = ValidationStatus.VALIDATED_CLEAN
        reason = (
            f"Validated with 0 bins pruned: initial Pearson R = {init_r:.4f} >= {args.target_r:.2f} "
            f"(Spearman rho = {init_rho:.4f}, {occupied_bins}/{n_bins} occupied bins)"
        )
    else:
        # Iterative bidirectional residual pruning loop
        iteration = 0
        min_required_bp = args.min_retained_cr_frac * float(total_cr_bp)

        while (
            curr_fit["pearson_r"] < args.target_r
            and curr_fit["slope"] > 0.0
            and int(np.sum(retained_mask)) > 3
        ):
            active_indices = np.where(retained_mask)[0]
            x_sub = x_all[active_indices]
            y_sub = y_all[active_indices]

            y_hat = curr_fit["intercept"] + curr_fit["slope"] * x_sub
            residuals = y_sub - y_hat
            abs_residuals = np.abs(residuals)

            # Select bin with maximum absolute vertical residual |Y_i - Y_hat_i|
            # Break ties deterministically by larger Y_i, then smaller bin_index
            best_local_idx = int(
                 max(
                    range(len(active_indices)),
                    key=lambda idx: (
                        round(float(abs_residuals[idx]), 9),
                        float(y_sub[idx]),
                        -int(bin_indices[active_indices[idx]]),
                    ),
                )
            )
            remove_global_idx = int(active_indices[best_local_idx])
            bp_to_remove = int(x_all[remove_global_idx])
            reads_to_remove = int(y_all[remove_global_idx])
            resid_val = float(residuals[best_local_idx])

            # Check if removing this bin would breach the 50% CR length retention floor
            if (curr_retained_bp - bp_to_remove) < min_required_bp - 1e-6:
                logger.info(
                    f"Stopping pruning before iteration {iteration + 1}: removing {bin_ids[remove_global_idx]} "
                    f"({bp_to_remove:,} bp) would drop retained CR sequence to "
                    f"{100.0 * (curr_retained_bp - bp_to_remove) / total_cr_bp:.2f}% (< {100.0 * args.min_retained_cr_frac:.1f}%)."
                )
                break

            iteration += 1
            r_before = curr_fit["pearson_r"]
            retained_mask[remove_global_idx] = False
            pruned_order[remove_global_idx] = iteration
            pruned_residual[remove_global_idx] = resid_val
            curr_retained_bp -= bp_to_remove
            curr_retained_reads -= reads_to_remove

            new_active = np.where(retained_mask)[0]
            curr_fit = fit_ols_and_correlations(x_all[new_active], y_all[new_active])
            r_after = curr_fit["pearson_r"]

            retained_cr_pct = (100.0 * curr_retained_bp / total_cr_bp) if total_cr_bp > 0 else 0.0
            retained_reads_pct = (100.0 * curr_retained_reads / total_mapped_pairs) if total_mapped_pairs > 0 else 0.0

            trajectory_rows.append({
                "iteration": iteration,
                "bin_id_removed": bin_ids[remove_global_idx],
                "bin_index_removed": bin_indices[remove_global_idx],
                "cr_bp_removed": bp_to_remove,
                "reads_removed": reads_to_remove,
                "vertical_residual": round(resid_val, 4),
                "retained_bins_count": int(len(new_active)),
                "retained_cr_bp": curr_retained_bp,
                "retained_cr_pct": round(retained_cr_pct, 4),
                "retained_reads": curr_retained_reads,
                "retained_reads_pct": round(retained_reads_pct, 4),
                "pearson_r": round(r_after, 6),
                "pearson_p": curr_fit["pearson_p"],
                "spearman_rho": round(curr_fit["spearman_rho"], 6),
                "ols_slope": curr_fit["slope"],
                "ols_intercept": curr_fit["intercept"],
            })

            logger.info(
                f"Pruning step {iteration:02d}: removed {bin_ids[remove_global_idx]} "
                f"(CR={bp_to_remove:,} bp, reads={reads_to_remove:,}, residual={resid_val:+.1f}) "
                f"-> R: {r_before:.4f} -> {r_after:.4f} (retained CR: {retained_cr_pct:.2f}%)"
            )

        n_pruned = len(pruned_order)
        final_r = curr_fit["pearson_r"]
        final_rho = curr_fit["spearman_rho"]
        final_slope = curr_fit["slope"]
        final_cr_pct = (100.0 * curr_retained_bp / total_cr_bp) if total_cr_bp > 0 else 0.0
        final_reads_pct = (100.0 * curr_retained_reads / total_mapped_pairs) if total_mapped_pairs > 0 else 0.0

        if final_r >= args.target_r and final_slope > 0.0 and final_cr_pct + 1e-6 >= (100.0 * args.min_retained_cr_frac):
            status = ValidationStatus.VALIDATED_PRUNED if n_pruned > 0 else ValidationStatus.VALIDATED_CLEAN
            reason = (
                f"Validated after pruning {n_pruned} outlier bin(s): Pearson R improved from "
                f"{init_r:.4f} to {final_r:.4f} >= {args.target_r:.2f} "
                f"(retained {final_cr_pct:.1f}% CR bp and {final_reads_pct:.1f}% mapped reads)"
            )
        else:
            status = ValidationStatus.FAILED_LOW_CORRELATION
            reason = (
                f"Failed correlation threshold: reached R = {final_r:.4f} (< {args.target_r:.2f}, "
                f"slope = {final_slope:.4e}) after pruning {n_pruned} bin(s) down to "
                f"{final_cr_pct:.1f}% retained CR length (min retention floor: {100.0 * args.min_retained_cr_frac:.1f}%)"
            )

    logger.info(f"[{status}] TaxID {args.taxid} ({args.species_name}): {reason}")

    # 4. Write {t_id}.bin_regression_data.tsv (per-bin table for scatter plot & inspection)
    out_bin_path = Path(args.out_bin_regression_tsv)
    out_bin_path.parent.mkdir(parents=True, exist_ok=True)
    num_segs_list = [int(x) for x in bins_df["num_cr_segments"].to_list()]
    num_contigs_list = [int(x) for x in bins_df["num_distinct_contigs"].to_list()]
    summary_list = bins_df["source_contigs_summary"].to_list()

    with open(out_bin_path, "w", encoding="utf-8", newline="") as bfh:
        writer = csv.DictWriter(
            bfh,
            fieldnames=[
                "t_id", "bin_id", "bin_index", "actual_cr_bp", "mapped_reads",
                "is_retained", "pruned_step", "vertical_residual_at_removal",
                "initial_y_hat", "final_y_hat", "num_cr_segments",
                "num_distinct_contigs", "source_contigs_summary",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for idx in range(n_bins):
            x_val = float(x_all[idx])
            y_val = int(y_all[idx])
            init_yhat = init_intercept + init_slope * x_val
            fin_yhat = curr_fit["intercept"] + curr_fit["slope"] * x_val
            writer.writerow({
                "t_id": args.taxid,
                "bin_id": bin_ids[idx],
                "bin_index": bin_indices[idx],
                "actual_cr_bp": int(x_all[idx]),
                "mapped_reads": y_val,
                "is_retained": bool(retained_mask[idx]),
                "pruned_step": pruned_order.get(idx, 0),
                "vertical_residual_at_removal": round(pruned_residual.get(idx, 0.0), 4),
                "initial_y_hat": round(init_yhat, 4),
                "final_y_hat": round(fin_yhat, 4),
                "num_cr_segments": num_segs_list[idx],
                "num_distinct_contigs": num_contigs_list[idx],
                "source_contigs_summary": summary_list[idx],
            })

    # 5. Write {t_id}.pruning_trajectory.tsv
    out_traj_path = Path(args.out_trajectory_tsv)
    out_traj_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_traj_path, "w", encoding="utf-8", newline="") as tfh:
        writer = csv.DictWriter(
            tfh,
            fieldnames=list(trajectory_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(trajectory_rows)

    # 6. Reconstruct 516-Week Time Series ({t_id}.weekly_timeseries_516w.tsv)
    retained_bin_set: Set[str] = {
        bin_ids[idx] for idx in range(n_bins) if retained_mask[idx]
    }

    all_cr_weekly_map: Dict[str, int] = {}
    clean_bin_weekly_map: Dict[str, int] = {}

    if total_mapped_pairs > 0:
        w_all_df = dedup_pairs_df.group_by("year_week").agg(pl.len().alias("cnt"))
        all_cr_weekly_map = dict(zip(w_all_df["year_week"].to_list(), w_all_df["cnt"].to_list()))

        w_clean_df = (
            dedup_pairs_df.filter(pl.col("bin_id").is_in(list(retained_bin_set)))
            .group_by("year_week")
            .agg(pl.len().alias("cnt"))
        )
        clean_bin_weekly_map = dict(zip(w_clean_df["year_week"].to_list(), w_clean_df["cnt"].to_list()))

    out_weekly_path = Path(args.out_weekly_tsv)
    out_weekly_path.parent.mkdir(parents=True, exist_ok=True)
    with open(args.raw_weekly_tsv, "r", encoding="utf-8") as rfh, open(
        out_weekly_path, "w", encoding="utf-8", newline=""
    ) as wfh:
        rreader = csv.DictReader(rfh, delimiter="\t")
        wwriter = csv.DictWriter(
            wfh,
            fieldnames=[
                "t_id", "t_species", "year_week", "year", "week",
                "raw_kraken_reads", "all_cr_mapped_reads", "validated_clean_bin_reads",
                "validation_status",
            ],
            delimiter="\t",
        )
        wwriter.writeheader()
        for rrow in rreader:
            yw = rrow["year_week"]
            raw_cnt = int(rrow["raw_kraken_reads"])
            all_cr_cnt = int(all_cr_weekly_map.get(yw, 0))
            clean_cnt = int(clean_bin_weekly_map.get(yw, 0))
            wwriter.writerow({
                "t_id": args.taxid,
                "t_species": args.species_name,
                "year_week": yw,
                "year": int(rrow["year"]),
                "week": int(rrow["week"]),
                "raw_kraken_reads": raw_cnt,
                "all_cr_mapped_reads": all_cr_cnt,
                "validated_clean_bin_reads": clean_cnt,
                "validation_status": status,
            })

    # 7. Write {t_id}.validation_metrics.json for fast Master Report ingestion
    map_stats: Dict[str, Any] = {}
    if Path(args.mapping_stats_tsv).is_file():
        with open(args.mapping_stats_tsv, "r", encoding="utf-8") as msfh:
            rows = list(csv.DictReader(msfh, delimiter="\t"))
            if rows:
                map_stats = rows[0]

    retained_cr_pct_final = (100.0 * curr_retained_bp / total_cr_bp) if total_cr_bp > 0 else 0.0
    retained_reads_pct_final = (100.0 * curr_retained_reads / total_mapped_pairs) if total_mapped_pairs > 0 else 0.0

    metrics = {
        "t_id": args.taxid,
        "species_name": args.species_name,
        "species_slug": args.species_slug,
        "n_bins": n_bins,
        "occupied_bins": occupied_bins,
        "occupied_bins_frac": round(occupied_frac, 4),
        "total_cr_bp": total_cr_bp,
        "retained_cr_bp": curr_retained_bp,
        "retained_cr_pct": round(retained_cr_pct_final, 4),
        "total_mapped_read_pairs": total_mapped_pairs,
        "retained_mapped_read_pairs": curr_retained_reads,
        "retained_mapped_reads_pct": round(retained_reads_pct_final, 4),
        "bins_pruned_count": len(pruned_order),
        "initial_pearson_r": round(init_r, 6),
        "initial_pearson_p": init_p,
        "initial_spearman_rho": round(init_rho, 6),
        "initial_spearman_p": init_rho_p,
        "initial_ols_slope": init_slope,
        "initial_ols_intercept": init_intercept,
        "final_pearson_r": round(curr_fit["pearson_r"], 6),
        "final_pearson_p": curr_fit["pearson_p"],
        "final_spearman_rho": round(curr_fit["spearman_rho"], 6),
        "final_spearman_p": curr_fit["spearman_p"],
        "final_ols_slope": curr_fit["slope"],
        "final_ols_intercept": curr_fit["intercept"],
        "validation_status": status,
        "status_reason": reason,
        "mapping_stats": map_stats,
    }

    out_json_path = Path(args.out_metrics_json)
    out_json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json_path, "w", encoding="utf-8") as jfh:
        json.dump(metrics, jfh, indent=2)

    logger.info(f"Wrote Stage 5 regression tables, 516-week timeseries, and metrics JSON for TaxID {args.taxid}.")


if __name__ == "__main__":
    main()
