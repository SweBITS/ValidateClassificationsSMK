#!/usr/bin/env python3
"""
generate_enhanced_species_dashboard.py — Widescreen 9-Panel Publication Diagnostic Suite
========================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
Produces a widescreen (29 x 17 in, 300 DPI PNG + vector PDF) 9-panel static diagnostic
dashboard for each target species, designed for widescreen monitors and manuscript figures:

  - **Row 1 (Genomic Bin Regression, Pruning Mechanics, Outlier Quadrant & Attrition)**:
    * **Panel A (Bin CR Length vs. Mapped Reads + Slim Contig/CR-Length Band)**:
      - Slim top marginal strip (14% height) locked to the exact same X-axis limits as the
        main scatter plot (with the colorbar in a dedicated side axis so bars and dots align
        pixel-for-pixel).
      - Vertical whiskers on each bin dot show the **length-weighted intra-bin CR coverage
        interquartile range (`p25`–`p75`)** extrapolated to bin width $W_b$.
      - Dots (`●` retained, `▲` pruned) are colored by **516-week Temporal Cosine Similarity**
        (`0.0` to `1.0`, Viridis) against the leave-one-out clean-genome consensus vector.
    * **Panel B (Recursive Outlier Pruning Trajectory)**:
      - Dual-axis trajectory of Pearson $R$ and Retained CR Sequence % across pruning steps.
    * **Panel D (Outlier Diagnostic Quadrant: Excess Depth vs. Temporal Cosine Similarity)**:
      - Plots $\\log_2(\\text{Observed} / \\text{Expected Clean OLS Reads})$ vs. 516-week
        Temporal Cosine Similarity, with marker area scaled by the **Length-Adjusted Top-1 CR
        Concentration Score** ($\\Delta f_{\\max, b} = \\max_{c \\in b}(\\text{Reads}_c/\\text{Reads}_b - \\text{Length}_c/W_b)$)
        and an explicit 3-size reference legend (`5%`, `25%`, `≥50%`).
    * **Panel E1 (True Symmetrical Read Attrition Funnel)** & **Panel E2 (`[90%, 100%]` Alignment Identity Histogram)**:
      - E1 renders a centered trapezoidal funnel across the 5 validation stages.
      - E2 renders the 1%-resolution sequence identity histogram (`<90%`, `90%..98%` rejected
        vs. `99%` and `100%` passing, split by Retained vs. Pruned bins) for reads passing
        `MAPQ ≥ 30` and `≥ 80%` aligned read length.

  - **Row 2 & Row 3 Left/Center (Wide Stacked Genomic & Temporal Panels + Log 1:1 + Simplified Phenology)**:
    * **Panel C (Wide Genome-Wide CR Manhattan Plot, `bin_001` → `bin_100`)**:
      - Every CR segment across the assembly in FASTA length order (`shortest → longest`),
        plotting `Reads / kb` on a `symlog` scale (`linear 0–1, log10 > 1`, equivalent to
        $\\log_{10}(\\text{Reads/kb} + 1)$) and colored by **516-week Temporal Cosine Similarity**
        (computed for all CRs with $\\ge 1$ read pair).
    * **Panel G (Wide 10-Year / 516-Week Longitudinal Time Series, stacked below Panel C)**:
      - Compares weekly Raw Kraken Clade Pairs, All-CR Mapped Pairs, and Validated Clean-Bin Pairs.
    * **Panel F (Standalone 1:1 Weekly Sample Correlation on `log10(Reads + 1)` Scale)**:
      - Reveals all 516 weekly samples across 6 orders of magnitude alongside the `1:1` line,
        proportionality curve ($Y = \\beta X$), linear $R$, log-scale $R_{\\log}$, and Spearman $\\rho$.
    * **Panel H (Simplified 10-Year Seasonal Phenology, Weeks 1–52, 3w Circular Smooth)**:
      - Compares only the 3-week circular-smoothed (`[0.25, 0.50, 0.25]` wrapping Week 52 $\\leftrightarrow$ Week 1)
        Pre-Validation vs. Post-Validation Clean seasonal curves, with a dual Y-axis showing
        both **Share of Annual Signal (`% / week`)** (left) and **Mean Weekly Clean Reads (`pairs / week`)** (right).

  - **Row 2 & Row 3 Right (`Panel I`: Hierarchical Clustering Heatmap of CR Bins Over Time)**:
    * Clusters all 100 CR bins (`bin_001..bin_100`) by their temporal profiles across active
      monitoring weeks using hierarchical linkage (`scipy.cluster.hierarchy`), displaying:
      1. Left dendrogram of the 100 bins,
      2. Two annotation columns (`1. Cosine Similarity` and `2. Pruned Status`), and
      3. The `100 Bins × Active Weeks` heatmap of $\\log_{10}(\\text{Reads / kb} + 1)$.
"""

import argparse
import csv
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.transforms as mtransforms
from matplotlib.patches import Polygon, Patch
import numpy as np
import polars as pl
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import ValidationStatus, setup_logger


STATUS_PALETTE: Dict[str, Dict[str, str]] = {
    ValidationStatus.VALIDATED_CLEAN: {
        "primary": "#047857",
        "fill": "#a7f3d0",
        "badge_bg": "#d1fae5",
        "badge_fg": "#064e3b",
        "badge_border": "#047857",
    },
    ValidationStatus.VALIDATED_PRUNED: {
        "primary": "#10b981",
        "fill": "#d1fae5",
        "badge_bg": "#ecfdf5",
        "badge_fg": "#047857",
        "badge_border": "#10b981",
    },
    ValidationStatus.FAILED_LOW_CORRELATION: {
        "primary": "#dc2626",
        "fill": "#fecaca",
        "badge_bg": "#fee2e2",
        "badge_fg": "#991b1b",
        "badge_border": "#dc2626",
    },
    ValidationStatus.FAILED_LOW_OCCUPANCY: {
        "primary": "#dc2626",
        "fill": "#fecaca",
        "badge_bg": "#fee2e2",
        "badge_fg": "#991b1b",
        "badge_border": "#dc2626",
    },
}


def compute_ols_ci_band(
    x_obs: np.ndarray,
    y_obs: np.ndarray,
    x_grid: np.ndarray,
    confidence: float = 0.95,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute OLS fitted line and exact t-distribution 95% CI band."""
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


def weighted_quantiles(values: np.ndarray, weights: np.ndarray, qs: List[float]) -> List[float]:
    """
    Compute length-weighted quantiles (e.g. p25, median, p75) of CR coverage within a bin
    so that hundreds of tiny 100 bp zero-read fragments do not artificially collapse the
    coverage IQR of a 2 Mb bin to zero.
    """
    if len(values) == 0:
        return [0.0 for _ in qs]
    order = np.argsort(values)
    v_sorted = values[order]
    w_sorted = np.maximum(1e-9, weights[order])
    cum_w = np.cumsum(w_sorted)
    tot_w = cum_w[-1]
    norm_w = (cum_w - 0.5 * w_sorted) / tot_w
    return [float(np.interp(q, norm_w, v_sorted)) for q in qs]


def parse_cr_fasta_headers(cr_fasta_path: Path) -> pl.DataFrame:
    """Fast scan of `{tid}.cr_contigs.fa` headers for exact ordered CR coordinates and lengths."""
    cr_ids: List[str] = []
    bin_ids: List[str] = []
    contigs: List[str] = []
    starts: List[int] = []
    ends: List[int] = []
    lengths: List[int] = []
    cum_mids: List[float] = []

    running_bp = 0
    with open(cr_fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip().split()[0]
            parts = header.split("|")
            cr_id = parts[0]
            bin_id = "UNKNOWN"
            contig = "UNKNOWN"
            s0, e0 = 0, 0
            for p in parts[1:]:
                if p.startswith("bin="):
                    bin_id = p[4:]
                elif p.startswith("contig="):
                    contig = p[7:]
                elif p.startswith("coords="):
                    c_str = p[7:]
                    if "-" in c_str:
                        s_str, e_str = c_str.split("-", 1)
                        s0, e0 = int(s_str), int(e_str)
            seg_len = max(1, e0 - s0)
            mid_bp = running_bp + 0.5 * seg_len
            running_bp += seg_len

            cr_ids.append(cr_id)
            bin_ids.append(bin_id)
            contigs.append(contig)
            starts.append(s0)
            ends.append(e0)
            lengths.append(seg_len)
            cum_mids.append(mid_bp)

    return pl.DataFrame({
        "cr_id": cr_ids,
        "bin_id": bin_ids,
        "contig": contigs,
        "start0": starts,
        "end0": ends,
        "cr_len_bp": lengths,
        "cum_mid_bp": cum_mids,
    })


def circular_smooth_52w(vec_52: np.ndarray) -> np.ndarray:
    """
    Apply a 3-week circular kernel `[0.25, 0.50, 0.25]` across weeks 1..52,
    wrapping Week 52 <-> Week 1 across the New Year boundary.
    """
    padded = np.concatenate([[vec_52[-1]], vec_52, [vec_52[0]]])
    return 0.25 * padded[:-2] + 0.50 * padded[1:-1] + 0.25 * padded[2:]


def compute_species_enhanced_data(
    taxid: int,
    slug: str,
    project_root: Path,
    logger: Any,
) -> Dict[str, Any]:
    """Compute all bin-level, CR-level, funnel/identity, phenology, and clustermap metrics."""
    val_dir = project_root / "results" / "04_validation" / f"{taxid}_{slug}"
    cr_dir = project_root / "results" / "02_classifiable_regions" / str(taxid)
    map_dir = project_root / "results" / "03_mapping" / str(taxid)

    metrics_json_path = val_dir / f"{taxid}.validation_metrics.json"
    bin_reg_path = val_dir / f"{taxid}.bin_regression_data.tsv"
    traj_path = val_dir / f"{taxid}.pruning_trajectory.tsv"
    weekly_path = val_dir / f"{taxid}.weekly_timeseries_516w.tsv"
    cr_fasta_path = cr_dir / f"{taxid}.cr_contigs.fa"
    map_stats_path = map_dir / f"{taxid}.mapping_stats.tsv"
    parquet_path = map_dir / f"{taxid}.passing_alignments.parquet"
    lt99_hist_path = map_dir / f"{taxid}.identity_histogram_lt99.tsv"

    with open(metrics_json_path, "r", encoding="utf-8") as fh:
        metrics = json.load(fh)
    with open(bin_reg_path, "r", encoding="utf-8") as fh:
        bin_rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(traj_path, "r", encoding="utf-8") as fh:
        traj_rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(weekly_path, "r", encoding="utf-8") as fh:
        weekly_rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(map_stats_path, "r", encoding="utf-8") as fh:
        map_stats_rows = list(csv.DictReader(fh, delimiter="\t"))
        map_stats = map_stats_rows[0] if map_stats_rows else {}

    active_sp_path = project_root / "results" / "00_ledger" / "active_species.tsv"
    if active_sp_path.exists():
        with open(active_sp_path, "r", encoding="utf-8") as fh:
            for arow in csv.DictReader(fh, delimiter="\t"):
                if str(arow.get("t_id")) == str(taxid):
                    metrics["assembly_accession"] = arow.get("chosen_accession", arow.get("clean_ncbi_accession", "NA"))
                    metrics["kept_assembly_bp"] = int(arow.get("kept_assembly_bp") or 0)
                    break
    stage1_audit_path = project_root / "results" / "00_ledger" / "assembly_filtering_summary.tsv"
    if stage1_audit_path.exists():
        with open(stage1_audit_path, "r", encoding="utf-8") as fh:
            for s1row in csv.DictReader(fh, delimiter="\t"):
                if str(s1row.get("t_id")) == str(taxid):
                    metrics["raw_assembly_bp"] = int(s1row.get("raw_assembly_bp") or 0)
                    if not metrics.get("kept_assembly_bp"):
                        metrics["kept_assembly_bp"] = int(s1row.get("kept_assembly_bp") or 0)
                    break
    total_raw_from_weekly = sum(int(r["raw_kraken_reads"]) for r in weekly_rows)
    metrics["raw_eDNA_reads"] = max(int(metrics.get("raw_eDNA_reads", 0)), total_raw_from_weekly)

    year_weeks = [r["year_week"] for r in weekly_rows]
    yw_to_idx = {yw: i for i, yw in enumerate(year_weeks)}
    n_weeks = len(year_weeks)

    # 1. Parse CR headers in FASTA order
    logger.info(f"[{taxid}] Scanning CR headers from {cr_fasta_path.name}...")
    cr_df = parse_cr_fasta_headers(cr_fasta_path)

    # 2. Aggregate alignments per CR
    logger.info(f"[{taxid}] Aggregating alignments from {parquet_path.name}...")
    aln_lazy = pl.scan_parquet(parquet_path)

    cr_aln_df = (
        aln_lazy.group_by("cr_id")
        .agg([
            (pl.len() // 2).cast(pl.Int64).alias("mapped_pairs"),
            (pl.col("mismatch_frac").mean() * 100.0).alias("mean_mismatch_pct"),
            ((pl.col("nm") == 0).mean() * 100.0).alias("exact_match_pct"),
        ])
        .collect()
    )

    cr_full = (
        cr_df.join(cr_aln_df, on="cr_id", how="left")
        .with_columns([
            pl.col("mapped_pairs").fill_null(0),
            pl.col("mean_mismatch_pct").fill_null(0.0),
            pl.col("exact_match_pct").fill_null(100.0),
        ])
        .with_columns([
            (pl.col("mapped_pairs").cast(pl.Float64) / pl.col("cr_len_bp").cast(pl.Float64)).alias("reads_per_bp"),
            (pl.col("mapped_pairs").cast(pl.Float64) * 1000.0 / pl.col("cr_len_bp").cast(pl.Float64)).alias("reads_per_kb"),
        ])
    )

    # 3. Build 516-week matrix for all 100 bins (100 x 516)
    bin_yw_df = (
        aln_lazy.group_by(["bin_id", "year_week"])
        .agg((pl.len() // 2).cast(pl.Int64).alias("pairs"))
        .collect()
    )
    bin_ids_ordered = [r["bin_id"] for r in bin_rows]
    bin_to_idx = {bid: i for i, bid in enumerate(bin_ids_ordered)}
    bin_ts_matrix = np.zeros((len(bin_ids_ordered), n_weeks), dtype=np.float64)

    for row in bin_yw_df.iter_rows(named=True):
        b_i = bin_to_idx.get(row["bin_id"])
        w_i = yw_to_idx.get(row["year_week"])
        if b_i is not None and w_i is not None:
            bin_ts_matrix[b_i, w_i] = float(row["pairs"])

    is_retained_mask = np.array([r["is_retained"].lower() == "true" for r in bin_rows], dtype=bool)
    clean_consensus_516w = np.sum(bin_ts_matrix[is_retained_mask], axis=0)
    clean_norm = float(np.linalg.norm(clean_consensus_516w))

    bin_cos_sim = np.zeros(len(bin_ids_ordered), dtype=np.float64)
    bin_active_weeks = np.sum(bin_ts_matrix > 0, axis=1).astype(int)
    for b_i in range(len(bin_ids_ordered)):
        vec_b = bin_ts_matrix[b_i]
        norm_b = float(np.linalg.norm(vec_b))
        if norm_b <= 0.0:
            bin_cos_sim[b_i] = 0.0
            continue
        if is_retained_mask[b_i] and np.sum(is_retained_mask) > 1:
            ref_vec = clean_consensus_516w - vec_b
        else:
            ref_vec = clean_consensus_516w
        ref_norm = float(np.linalg.norm(ref_vec))
        if ref_norm > 0.0:
            bin_cos_sim[b_i] = float(np.clip(np.dot(vec_b, ref_vec) / (norm_b * ref_norm), 0.0, 1.0))

    # 4. Compute 516-week Temporal Cosine Similarity for ALL CRs with >= 1 mapped read pair
    cr_cos_map: Dict[str, float] = {}
    if clean_norm > 0.0:
        yw_weights_df = pl.DataFrame({
            "year_week": year_weeks,
            "ref_w": (clean_consensus_516w / clean_norm).tolist(),
        })
        cr_dots = (
            aln_lazy.group_by(["cr_id", "year_week"])
            .agg((pl.len() // 2).cast(pl.Float64).alias("pairs"))
            .join(yw_weights_df.lazy(), on="year_week", how="inner")
            .group_by("cr_id")
            .agg([
                (pl.col("pairs") * pl.col("ref_w")).sum().alias("dot_ref"),
                (pl.col("pairs") ** 2).sum().sqrt().alias("cr_norm"),
            ])
            .with_columns(
                pl.when(pl.col("cr_norm") > 0)
                .then((pl.col("dot_ref") / pl.col("cr_norm")).clip(0.0, 1.0))
                .otherwise(0.0)
                .alias("cos_sim")
            )
            .collect()
        )
        for r in cr_dots.iter_rows(named=True):
            cr_cos_map[r["cr_id"]] = float(r["cos_sim"])

    cr_full = cr_full.with_columns(
        pl.col("cr_id").replace_strict(cr_cos_map, default=0.0).alias("cos_sim_516w")
    )

    # 5. Per-bin identity breakdown
    bin_id_stats_df = (
        aln_lazy.group_by("bin_id")
        .agg([
            (pl.col("mismatch_frac").mean() * 100.0).alias("mean_mismatch_pct"),
            ((pl.col("nm") == 0).mean() * 100.0).alias("exact_match_pct"),
            (pl.col("align_frac").mean() * 100.0).alias("mean_align_pct"),
        ])
        .collect()
    )
    bin_id_stats_map = {r["bin_id"]: r for r in bin_id_stats_df.iter_rows(named=True)}

    # 6. Per-bin CR structural & length-weighted coverage quartiles
    bin_cr_stats_df = (
        cr_full.group_by("bin_id", maintain_order=True)
        .agg([
            pl.col("cr_len_bp").quantile(0.25, interpolation="linear").alias("cr_len_p25"),
            pl.col("cr_len_bp").median().alias("cr_len_median"),
            pl.col("cr_len_bp").quantile(0.75, interpolation="linear").alias("cr_len_p75"),
            pl.col("cr_len_bp").mean().alias("cr_len_mean"),
            pl.col("reads_per_bp").max().alias("cr_cov_max"),
            pl.col("mapped_pairs").max().alias("top1_cr_reads"),
            pl.col("cr_id").sort_by("mapped_pairs", descending=True).first().alias("top1_cr_id"),
            pl.col("contig").sort_by("mapped_pairs", descending=True).first().alias("top1_cr_contig"),
            pl.col("cr_len_bp").sort_by("mapped_pairs", descending=True).first().alias("top1_cr_len_bp"),
            pl.col("cum_mid_bp").min().alias("bin_start_cum_bp"),
            pl.col("cum_mid_bp").max().alias("bin_end_cum_bp"),
        ])
    )
    bin_cr_map = {r["bin_id"]: r for r in bin_cr_stats_df.iter_rows(named=True)}

    # Length-weighted intra-bin coverage quartiles (p25, median, p75) and Length-Weighted Gini Coefficient
    bin_w_quartiles: Dict[str, Tuple[float, float, float]] = {}
    bin_gini_map: Dict[str, float] = {}
    for bid_tuple, sub_df in cr_full.group_by("bin_id", maintain_order=True):
        bid_str = bid_tuple[0] if isinstance(bid_tuple, tuple) else str(bid_tuple)
        lens = sub_df["cr_len_bp"].to_numpy().astype(np.float64)
        pairs = sub_df["mapped_pairs"].to_numpy().astype(np.float64)
        # Partition the bin's CRs into ~12 equal-rank sub-windows to measure intra-bin coverage IQR cleanly
        n_segs = len(lens)
        n_chunks = min(12, max(4, n_segs // 8))
        if n_segs >= 4:
            chunks_l = np.array_split(lens, n_chunks)
            chunks_p = np.array_split(pairs, n_chunks)
            c_lens = np.array([np.sum(cl) for cl in chunks_l], dtype=np.float64)
            c_covs = np.array([np.sum(cp) / max(1.0, np.sum(cl)) for cp, cl in zip(chunks_p, chunks_l)], dtype=np.float64)
            q25, q50, q75 = weighted_quantiles(c_covs, c_lens, [0.25, 0.50, 0.75])
        else:
            covs = pairs / np.maximum(1.0, lens)
            q25, q50, q75 = weighted_quantiles(covs, lens, [0.25, 0.50, 0.75])
        bin_w_quartiles[bid_str] = (q25, q50, q75)

        # Length-Weighted Gini Coefficient G_b in [0, 1] across CRs in bin b:
        # G_b = 0 if read density (reads/bp) is uniform across all CRs in the bin;
        # G_b -> 1 if a tiny fraction of CR base pairs monopolizes the bin's reads (whether 1 CR or 5 CRs).
        tot_l = float(np.sum(lens))
        tot_p = float(np.sum(pairs))
        if tot_l > 0.0 and tot_p > 0.0 and n_segs >= 2:
            dens = pairs / np.maximum(1.0, lens)
            ord_g = np.argsort(dens)
            w_l = lens[ord_g] / tot_l
            w_p = pairs[ord_g] / tot_p
            cum_p = np.cumsum(w_p)
            cum_p_prev = np.concatenate([[0.0], cum_p[:-1]])
            lorenz_area = float(np.sum(0.5 * (cum_p_prev + cum_p) * w_l))
            gini_b = float(np.clip(1.0 - 2.0 * lorenz_area, 0.0, 1.0))
        else:
            gini_b = 0.0
        bin_gini_map[bid_str] = gini_b

    # Excess top-1 CR share Delta f_{max, b}
    cr_excess_df = (
        cr_full.join(
            pl.DataFrame({
                "bin_id": [r["bin_id"] for r in bin_rows],
                "bin_cr_bp": [float(r["actual_cr_bp"]) for r in bin_rows],
                "bin_reads": [float(r["mapped_reads"]) for r in bin_rows],
            }),
            on="bin_id",
            how="left",
        )
        .with_columns([
            pl.when(pl.col("bin_reads") > 0)
            .then(pl.col("mapped_pairs").cast(pl.Float64) / pl.col("bin_reads"))
            .otherwise(0.0)
            .alias("read_share"),
            (pl.col("cr_len_bp").cast(pl.Float64) / pl.col("bin_cr_bp")).alias("len_share"),
        ])
        .with_columns((pl.col("read_share") - pl.col("len_share")).alias("excess_share"))
        .group_by("bin_id")
        .agg([
            pl.col("read_share").max().alias("f_max_raw"),
            pl.col("excess_share").max().alias("delta_f_max"),
        ])
    )
    cr_excess_map = {r["bin_id"]: r for r in cr_excess_df.iter_rows(named=True)}

    # Length-proportional clean CR read density (reads per bp across retained bins) for Panel D:
    # Enforces E(L=0) = 0 and strictly non-negative expectation even when unconstrained OLS has
    # a negative slope (e.g. single-bin spike failures) or a negative Y-intercept.
    ret_cr_bp_sum = sum(float(r["actual_cr_bp"]) for r in bin_rows if r["is_retained"].lower() == "true")
    ret_reads_sum = sum(float(r["mapped_reads"]) for r in bin_rows if r["is_retained"].lower() == "true")
    if ret_cr_bp_sum <= 0:
        ret_cr_bp_sum = max(1.0, sum(float(r["actual_cr_bp"]) for r in bin_rows))
        ret_reads_sum = sum(float(r["mapped_reads"]) for r in bin_rows)
    clean_density_per_bp = ret_reads_sum / max(1.0, ret_cr_bp_sum)

    enhanced_bins: List[Dict[str, Any]] = []
    for b_i, r in enumerate(bin_rows):
        bid = r["bin_id"]
        actual_bp = float(r["actual_cr_bp"])
        mapped_r = float(r["mapped_reads"])
        exp_reads_prop = max(1.0, actual_bp * clean_density_per_bp)
        cr_st = bin_cr_map.get(bid, {})
        id_st = bin_id_stats_map.get(bid, {})
        ex_st = cr_excess_map.get(bid, {})
        cov_p25, cov_med, cov_p75 = bin_w_quartiles.get(bid, (0.0, 0.0, 0.0))
        cov_max = float(cr_st.get("cr_cov_max", 0.0))
        gini_val = float(bin_gini_map.get(bid, 0.0))

        log2_oe = math.log2(max(0.125, (mapped_r + 1.0) / (exp_reads_prop + 1.0)))
        occ_pct = 100.0 * float(bin_active_weeks[b_i]) / max(1, n_weeks)

        rec = dict(r)
        rec.update({
            "cr_len_p25_bp": round(float(cr_st.get("cr_len_p25", 0.0)), 1),
            "cr_len_median_bp": round(float(cr_st.get("cr_len_median", 0.0)), 1),
            "cr_len_p75_bp": round(float(cr_st.get("cr_len_p75", 0.0)), 1),
            "cr_len_mean_bp": round(float(cr_st.get("cr_len_mean", 0.0)), 1),
            "cr_cov_p25_per_kb": round(cov_p25 * 1000.0, 4),
            "cr_cov_median_per_kb": round(cov_med * 1000.0, 4),
            "cr_cov_p75_per_kb": round(cov_p75 * 1000.0, 4),
            "cr_cov_max_per_kb": round(cov_max * 1000.0, 4),
            "extrap_reads_p25": round(actual_bp * cov_p25, 1),
            "extrap_reads_median": round(actual_bp * cov_med, 1),
            "extrap_reads_p75": round(actual_bp * cov_p75, 1),
            "gini_coeff": round(gini_val, 4),
            "f_max_raw": round(float(ex_st.get("f_max_raw", 0.0)), 4),
            "delta_f_max": round(max(0.0, float(ex_st.get("delta_f_max", 0.0))), 4),
            "top1_cr_id": str(cr_st.get("top1_cr_id", "NA")),
            "top1_cr_contig": str(cr_st.get("top1_cr_contig", "NA")),
            "top1_cr_reads": int(cr_st.get("top1_cr_reads", 0)),
            "top1_cr_len_bp": int(cr_st.get("top1_cr_len_bp", 0)),
            "cos_sim_516w": round(float(bin_cos_sim[b_i]), 4),
            "active_weeks": int(bin_active_weeks[b_i]),
            "sample_occupancy_pct": round(occ_pct, 3),
            "mean_mismatch_pct": round(float(id_st.get("mean_mismatch_pct", 0.0)), 4),
            "exact_match_pct": round(float(id_st.get("exact_match_pct", 100.0)), 2),
            "log2_obs_exp": round(log2_oe, 4),
            "bin_start_cum_mb": round(float(cr_st.get("bin_start_cum_bp", 0.0)) / 1e6, 4),
            "bin_end_cum_mb": round(float(cr_st.get("bin_end_cum_bp", 0.0)) / 1e6, 4),
        })
        enhanced_bins.append(rec)

    enhanced_tsv_path = val_dir / f"{taxid}.enhanced_bin_metrics.tsv"
    if enhanced_bins:
        with open(enhanced_tsv_path, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(enhanced_bins[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(enhanced_bins)

    # 7. Read Attrition Funnel & [90%, 100%] Identity Histogram (Panel E1 & E2)
    raw_pairs = int(metrics.get("raw_eDNA_reads", 0))
    total_sam = int(map_stats.get("total_sam_records", raw_pairs * 2))
    prim_mapped_mates = int(map_stats.get("primary_mapped_records", 0))
    failed_mapq_mates = int(map_stats.get("failed_mapq", 0))
    failed_align_mates = int(map_stats.get("failed_align_frac", 0))
    failed_mm_mates = int(map_stats.get("failed_mismatch_frac", 0))
    passing_pairs = int(map_stats.get("unique_passing_read_pairs", 0))
    retained_pairs = int(metrics.get("retained_mapped_read_pairs", 0))

    stage1_raw_pairs = float(max(raw_pairs, total_sam // 2))
    stage2_mapq30_pairs = (prim_mapped_mates - failed_mapq_mates) / 2.0
    stage3_align80_pairs = (prim_mapped_mates - failed_mapq_mates - failed_align_mates) / 2.0
    stage4_both_mates_pairs = float(passing_pairs)
    stage5_clean_pairs = float(retained_pairs)
    rejected_lt99_pairs = failed_mm_mates / 2.0

    retained_bin_set = {r["bin_id"] for r in bin_rows if r["is_retained"].lower() == "true"}
    id_split_df = (
        aln_lazy.with_columns([
            pl.col("bin_id").is_in(list(retained_bin_set)).alias("is_ret"),
            (pl.col("nm") == 0).alias("is_exact"),
        ])
        .group_by(["is_ret", "is_exact"])
        .agg((pl.len() / 2.0).alias("pair_equiv"))
        .collect()
    )
    id_split = {(bool(r["is_ret"]), bool(r["is_exact"])): float(r["pair_equiv"]) for r in id_split_df.iter_rows(named=True)}

    # Load <99% identity distribution [89 (<90%), 90, 91, ..., 98]:
    # Prefer exact 516-week columns from mapping_stats.tsv if present; otherwise load sidecar TSV
    lt99_bins = {b: 0.0 for b in range(89, 99)}
    if "id_lt90_pairs" in map_stats:
        lt99_bins[89] = float(map_stats.get("id_lt90_pairs", 0))
        for b in range(90, 99):
            lt99_bins[b] = float(map_stats.get(f"id_{b}_pairs", 0))
    elif lt99_hist_path.exists():
        with open(lt99_hist_path, "r", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                b = int(r.get("id_bin", r.get("id_bin_floor", 0)))
                if b in lt99_bins:
                    lt99_bins[b] = float(r.get("pair_count", r.get("sampled_mates", 0.0)))
    tot_sampled_lt99 = sum(lt99_bins.values())
    if tot_sampled_lt99 > 0 and rejected_lt99_pairs > 0:
        for b in lt99_bins:
            lt99_bins[b] = rejected_lt99_pairs * (lt99_bins[b] / tot_sampled_lt99)

    funnel_data = {
        "stage_labels": [
            "1. Raw Kraken Clade Pairs",
            "2. CR Mapped (MAPQ ≥ 30)",
            "3. Passed ≥ 80% Aligned",
            "4. Passed ≥ 99% Identity",
            "5. Clean Validated Bins",
        ],
        "stage_counts": [
            stage1_raw_pairs,
            stage2_mapq30_pairs,
            stage3_align80_pairs,
            stage4_both_mates_pairs,
            stage5_clean_pairs,
        ],
        "rejected_lt99_pairs": rejected_lt99_pairs,
        "lt99_hist_pairs": lt99_bins,
        "ret_exact_100_pairs": id_split.get((True, True), 0.0),
        "ret_1mm_99_pairs": id_split.get((True, False), 0.0),
        "pruned_exact_100_pairs": id_split.get((False, True), 0.0),
        "pruned_1mm_99_pairs": id_split.get((False, False), 0.0),
    }

    # 8. Simplified 10-Year Seasonal Phenology (Weeks 1..52, Circular 3-Week Smoothing)
    years_all = sorted(set(int(r["year"]) for r in weekly_rows))
    raw_yr_wk = np.zeros((len(years_all), 52), dtype=np.float64)
    clean_yr_wk = np.zeros((len(years_all), 52), dtype=np.float64)
    yr_to_i = {y: i for i, y in enumerate(years_all)}

    for r in weekly_rows:
        yi = yr_to_i[int(r["year"])]
        wi = min(51, max(0, int(r["week"]) - 1))
        raw_yr_wk[yi, wi] += float(r["raw_kraken_reads"])
        clean_yr_wk[yi, wi] += float(r["validated_clean_bin_reads"])

    raw_52_sum = circular_smooth_52w(np.sum(raw_yr_wk, axis=0))
    clean_52_sum = circular_smooth_52w(np.sum(clean_yr_wk, axis=0))

    raw_share_52 = 100.0 * raw_52_sum / max(1.0, float(np.sum(raw_52_sum)))
    clean_share_52 = 100.0 * clean_52_sum / max(1.0, float(np.sum(clean_52_sum)))
    clean_mean_reads_per_wk = clean_52_sum / max(1, len(years_all))

    if np.std(raw_share_52) > 0 and np.std(clean_share_52) > 0:
        r_pheno = float(stats.pearsonr(raw_share_52, clean_share_52)[0])
    else:
        r_pheno = 0.0

    phenology_data = {
        "weeks": list(range(1, 53)),
        "raw_share": raw_share_52.tolist(),
        "clean_share": clean_share_52.tolist(),
        "clean_mean_reads": clean_mean_reads_per_wk.tolist(),
        "r_pheno": r_pheno,
    }

    return {
        "taxid": taxid,
        "slug": slug,
        "metrics": metrics,
        "enhanced_bins": enhanced_bins,
        "traj_rows": traj_rows,
        "weekly_rows": weekly_rows,
        "cr_full": cr_full,
        "bin_ts_matrix": bin_ts_matrix,
        "funnel_data": funnel_data,
        "phenology_data": phenology_data,
        "val_dir": val_dir,
    }


def fmt_compact(v: float) -> str:
    """Format large counts compactly (e.g. 7.15M, 428k, 520)."""
    av = abs(v)
    if av >= 1e6:
        return f"{v / 1e6:.2f}M"
    if av >= 1e4:
        return f"{v / 1e3:.0f}k"
    if av >= 1e3:
        return f"{v / 1e3:.1f}k"
    return f"{int(round(v)):,}"


def fmt_bp(bp: float) -> str:
    """Format genomic base-pair sizes cleanly (e.g. 1,035.9 Mb, 319.2 Mb, 420.5 kb)."""
    abp = abs(bp)
    if abp >= 1e6:
        return f"{bp / 1e6:,.1f} Mb"
    if abp >= 1e3:
        return f"{bp / 1e3:,.1f} kb"
    return f"{int(round(bp)):,} bp"


def render_10panel_figure(
    data: Dict[str, Any],
    out_png: Path,
    out_pdf: Path,
    target_r: float = 0.90,
    min_retained_cr_frac: float = 0.50,
) -> None:
    """
    Render the 10-panel diagnostic dashboard (`{t_id}.diagnostic_10panel.{png,pdf}`, 32.0 x 18.2 in, 300 DPI).
    Uses a unified 3x4 GridSpec so that:
      - Panel A is the exact same width as Panel C and Panel G (Col 0)
      - Panel B is the exact same width as Panel F (Col 1)
      - Panel D is the exact same width as Panel I (Col 2)
      - Panel H spans Col 1 + Col 2 (Row 2, Cols 1..2)
      - Panel E and Panel J share Col 3 (Row 0 vs. Rows 1..2)
    """
    taxid = data["taxid"]
    metrics = data["metrics"]
    bin_rows = data["enhanced_bins"]
    traj_rows = data["traj_rows"]
    weekly_rows = data["weekly_rows"]
    cr_full: pl.DataFrame = data["cr_full"]
    bin_ts_matrix: np.ndarray = data["bin_ts_matrix"]
    funnel = data["funnel_data"]
    pheno = data["phenology_data"]

    status = metrics.get("validation_status", ValidationStatus.FAILED_LOW_CORRELATION)
    palette = STATUS_PALETTE.get(
        status,
        {
            "primary": "#475569",
            "fill": "#cbd5e1",
            "badge_bg": "#f1f5f9",
            "badge_fg": "#475569",
            "badge_border": "#cbd5e1",
        },
    )

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "axes.edgecolor": "#cbd5e1",
        "axes.linewidth": 0.9,
        "grid.color": "#e2e8f0",
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
    })

    fig = plt.figure(figsize=(32.0, 18.2), dpi=300, facecolor="white")

    # Unified 3x4 Master GridSpec guarantees exact column width alignment across all 3 rows:
    #   Col 0 (1.54): A (Row 0), C (Row 1), G (Row 2)
    #   Col 1 (0.78): B (Row 0), F (Row 1)
    #   Col 2 (0.84): D (Row 0), I (Row 1); H spans Cols 1..2 on Row 2
    #   Col 3 (1.34): E (Row 0), J (Rows 1..2)
    master_gs = gridspec.GridSpec(
        3, 4,
        figure=fig,
        height_ratios=[1.08, 0.96, 0.96],
        width_ratios=[1.54, 0.78, 0.84, 1.34],
        hspace=0.43,
        wspace=0.20,
        left=0.038,
        right=0.978,
        top=0.856,
        bottom=0.048,
    )

    def set_panel_header(
        ax: plt.Axes,
        title_str: str,
        subtitle_str: str = "",
        title_pts: float = 17.5,
        sub_pts: float = 4.2,
        title_size: float = 12.0,
        sub_size: float = 8.7,
    ) -> None:
        """Render a large bold panel title and a concise explanatory subtitle at fixed physical point offsets."""
        if subtitle_str:
            trans_title = ax.transAxes + mtransforms.ScaledTranslation(0.0, title_pts / 72.0, fig.dpi_scale_trans)
            trans_sub = ax.transAxes + mtransforms.ScaledTranslation(0.0, sub_pts / 72.0, fig.dpi_scale_trans)
            ax.text(
                0.0, 1.0, title_str, transform=trans_title,
                fontsize=title_size, fontweight="bold", color="#0f172a", ha="left", va="bottom", clip_on=False
            )
            ax.text(
                0.0, 1.0, subtitle_str, transform=trans_sub,
                fontsize=sub_size, fontweight="normal", color="#475569", ha="left", va="bottom", clip_on=False
            )
        else:
            trans_single = ax.transAxes + mtransforms.ScaledTranslation(0.0, 7.0 / 72.0, fig.dpi_scale_trans)
            ax.text(
                0.0, 1.0, title_str, transform=trans_single,
                fontsize=title_size, fontweight="bold", color="#0f172a", ha="left", va="bottom", clip_on=False
            )

    # --- Header Banner (Italic Species Name + TaxID + Report-Matched Colored Status Pill) --- #
    sci_name = metrics.get("species_name", f"TaxID {taxid}")
    raw_w = np.array([float(r["raw_kraken_reads"]) for r in weekly_rows], dtype=np.float64)
    all_cr_w = np.array([float(r["all_cr_mapped_reads"]) for r in weekly_rows], dtype=np.float64)
    clean_w = np.array([float(r["validated_clean_bin_reads"]) for r in weekly_rows], dtype=np.float64)
    r_wk_lin = float(stats.pearsonr(raw_w, clean_w)[0]) if (np.std(raw_w) > 0 and np.std(clean_w) > 0) else 0.0

    total_cr_bp = int(metrics.get("total_cr_bp", 0))
    retained_cr_bp = int(metrics.get("retained_cr_bp", 0))
    total_asm_bp = int(metrics.get("raw_assembly_bp") or metrics.get("kept_assembly_bp") or total_cr_bp)

    y_hdr = 0.962
    t_sp = fig.text(
        0.5, y_hdr, sci_name,
        fontsize=20.0, fontstyle="italic", fontweight="bold", color="#0f172a", va="center", ha="left"
    )
    t_tid = fig.text(
        0.5, y_hdr, f"  (TaxID: {taxid})   —   ",
        fontsize=18.0, fontweight="bold", color="#334155", va="center", ha="left"
    )
    t_badge = fig.text(
        0.5, y_hdr, f" {status} ",
        fontsize=14.5, fontweight="bold", color=palette["badge_fg"], va="center", ha="left",
        bbox=dict(
            boxstyle="round,pad=0.28,rounding_size=0.40",
            facecolor=palette["badge_bg"],
            edgecolor=palette["badge_border"],
            linewidth=1.6,
        ),
    )
    # Center the 3-piece header [Italic Species + TaxID + Colored Status Pill] at x = 0.50
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv_fig = fig.transFigure.inverted()
    w_sp = inv_fig.transform_bbox(t_sp.get_window_extent(renderer=renderer)).width
    w_tid = inv_fig.transform_bbox(t_tid.get_window_extent(renderer=renderer)).width
    w_badge = inv_fig.transform_bbox(t_badge.get_window_extent(renderer=renderer)).width
    tot_hdr_w = w_sp + w_tid + w_badge
    x_start_hdr = 0.50 - tot_hdr_w / 2.0
    t_sp.set_position((x_start_hdr, y_hdr))
    t_tid.set_position((x_start_hdr + w_sp, y_hdr))
    t_badge.set_position((x_start_hdr + w_sp + w_tid, y_hdr))

    subtitle_1 = (
        f"Reference Genome: {metrics.get('assembly_accession', 'NA')} ({fmt_bp(total_asm_bp)})   |   "
        f"Classifiable Sequence (CR): {fmt_bp(total_cr_bp)} ({100.0 * total_cr_bp / max(1, total_asm_bp):.1f}% of genome; "
        f"{metrics.get('retained_cr_pct', 0):.1f}% retained across "
        f"{100 - int(metrics.get('bins_pruned_count', 0))}/100 clean bins; {cr_full.height:,} CR segments)   |   "
        f"eDNA Read Pairs: {int(metrics.get('raw_eDNA_reads', 0)):,} raw  →  "
        f"{int(metrics.get('total_mapped_read_pairs', 0)):,} CR-mapped  →  "
        f"{int(metrics.get('retained_mapped_read_pairs', 0)):,} clean ({metrics.get('retained_mapped_reads_pct', 0):.1f}%)"
    )
    subtitle_2 = (
        f"Genomic Proportionality: Initial R = {metrics.get('initial_pearson_r', 0):.3f}  →  "
        f"Final R = {metrics.get('final_pearson_r', 0):.3f} (Spearman ρ = {metrics.get('final_spearman_rho', 0):.3f}; "
        f"{metrics.get('bins_pruned_count', 0)} outlier bins removed)   |   "
        f"Weekly Signal Retention: R = {r_wk_lin:.3f}   |   "
        f"Seasonal Phenology Concordance: r = {pheno['r_pheno']:.3f}"
    )
    fig.text(0.5, 0.927, subtitle_1, ha="center", va="center", fontsize=11.0, color="#1e293b")
    fig.text(0.5, 0.903, subtitle_2, ha="center", va="center", fontsize=10.4, color="#475569")

    cmap_cos = plt.cm.viridis
    norm_cos = mcolors.Normalize(vmin=0.0, vmax=1.0)

    # ======================================================================== #
    # PANEL A: Genomic Proportionality Regression (Row 0, Col 0 — Same Width as C & G)
    # ======================================================================== #
    gs_a = gridspec.GridSpecFromSubplotSpec(
        3, 1,
        subplot_spec=master_gs[0, 0],
        height_ratios=[0.145, 0.725, 0.130],
        hspace=0.04,
    )
    ax_a_top = fig.add_subplot(gs_a[0, 0])
    ax_a = fig.add_subplot(gs_a[1, 0], sharex=ax_a_top)

    x_all = np.array([float(r["actual_cr_bp"]) for r in bin_rows], dtype=np.float64)
    y_all = np.array([float(r["mapped_reads"]) for r in bin_rows], dtype=np.float64)
    is_ret = np.array([r["is_retained"].lower() == "true" for r in bin_rows], dtype=bool)
    pruned_steps = [int(r["pruned_step"]) for r in bin_rows]
    cos_all = np.array([float(r["cos_sim_516w"]) for r in bin_rows], dtype=np.float64)
    n_contigs = np.array([int(r["num_distinct_contigs"]) for r in bin_rows], dtype=np.float64)
    med_cr_kb = np.array([float(r["cr_len_median_bp"]) / 1000.0 for r in bin_rows], dtype=np.float64)
    y_p25 = np.array([float(r["extrap_reads_p25"]) for r in bin_rows], dtype=np.float64)
    y_p75 = np.array([float(r["extrap_reads_p75"]) for r in bin_rows], dtype=np.float64)

    x_span = max(1.0, float(np.max(x_all) - np.min(x_all))) if len(x_all) > 1 else 1000.0
    x_min_a = float(np.min(x_all)) - 0.035 * x_span
    x_max_a = float(np.max(x_all)) + 0.035 * x_span

    # Slim Top Strip: Distinct Contigs (bars) & Median CR Length kb (blue line)
    ax_a_top.grid(True, alpha=0.45, zorder=0)
    bar_w = x_span / 135.0
    ax_a_top.bar(
        x_all, n_contigs, width=bar_w, align="center", color="#64748b", alpha=0.65, zorder=2
    )
    ax_a_top.set_ylabel("Contigs", fontsize=7.6, color="#475569")
    if np.max(n_contigs) > 12:
        ax_a_top.set_yscale("symlog", linthresh=2.0)
    ax_a_top.tick_params(axis="x", labelbottom=False, bottom=False)
    ax_a_top.tick_params(axis="y", labelsize=7.0)

    ax_a_top_r = ax_a_top.twinx()
    ax_a_top_r.plot(x_all, med_cr_kb, color="#0284c7", linewidth=1.3, linestyle="-", zorder=3)
    ax_a_top_r.set_ylabel("Med CR (kb)", fontsize=7.6, color="#0284c7")
    ax_a_top_r.tick_params(axis="y", labelsize=7.0, colors="#0284c7")
    set_panel_header(
        ax_a_top,
        "A. Genomic Proportionality Regression",
        "Mapped read pairs vs. classifiable region (CR) length per bin (top = contigs & median CR kb; whiskers = intra-bin IQR)",
    )

    # Main Scatter Plot in Panel A
    ax_a.grid(True, alpha=0.7, zorder=0)
    if len(x_all) > 0 and np.max(x_all) > 0:
        x_grid = np.linspace(x_min_a, x_max_a, 200)
        init_slope = float(metrics.get("initial_ols_slope", 0.0))
        init_intercept = float(metrics.get("initial_ols_intercept", 0.0))
        ax_a.plot(
            x_grid, init_intercept + init_slope * x_grid,
            linestyle="--", linewidth=1.4, color="#64748b",
            label=f"Initial OLS (R={metrics.get('initial_pearson_r', 0):.3f})", zorder=2
        )

        x_ret = x_all[is_ret]
        y_ret = y_all[is_ret]
        if len(x_ret) >= 3:
            y_fin_grid, ci_low, ci_high = compute_ols_ci_band(x_ret, y_ret, x_grid)
            ax_a.fill_between(
                x_grid, ci_low, ci_high, color=palette["fill"], alpha=0.45,
                label="95% CI (clean bins)", zorder=1
            )
            ax_a.plot(
                x_grid, y_fin_grid, linestyle="-", linewidth=2.1, color=palette["primary"],
                label=f"Final OLS (R={metrics.get('final_pearson_r', 0):.3f})", zorder=3
            )

        # Crisp, visible intra-bin p25-p75 coverage whiskers with horizontal caps
        err_low = np.maximum(0.0, y_all - np.minimum(y_all, y_p25))
        err_high = np.maximum(0.0, np.maximum(y_all, y_p75) - y_all)
        ax_a.errorbar(
            x_all, y_all, yerr=[err_low, err_high],
            fmt="none", ecolor="#334155", elinewidth=1.2, capsize=2.4, capthick=1.0,
            alpha=0.75, label="Intra-bin coverage IQR (p25–p75)", zorder=3.5
        )

        sc_a = ax_a.scatter(
            x_ret, y_ret, s=46, c=cos_all[is_ret], cmap=cmap_cos, norm=norm_cos,
            edgecolors="#0f172a", linewidths=0.5, alpha=0.92, marker="o",
            label=f"Retained bins (n={len(x_ret)})", zorder=4
        )

        if np.any(~is_ret):
            ax_a.scatter(
                x_all[~is_ret], y_all[~is_ret], s=74, c=cos_all[~is_ret], cmap=cmap_cos, norm=norm_cos,
                edgecolors="#e11d48", linewidths=1.4, alpha=0.96, marker="^",
                label=f"Pruned bins (n={int(np.sum(~is_ret))})", zorder=5
            )
            pruned_indices = sorted([i for i, flag in enumerate(is_ret) if not flag], key=lambda idx: pruned_steps[idx])
            for idx in pruned_indices[:15]:
                ax_a.annotate(
                    f"#{pruned_steps[idx]}",
                    xy=(x_all[idx], y_all[idx]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=7.3, fontweight="bold", color="#be123c", zorder=6
                )

        # Horizontal Colorbar positioned cleanly in the space between Panel A and Panel C
        ax_a_cbar = ax_a.inset_axes([0.12, -0.245, 0.76, 0.040])
        cbar_a = fig.colorbar(sc_a, cax=ax_a_cbar, orientation="horizontal")
        cbar_a.set_label("516-Week Temporal Concordance (Cosine Similarity — Shared by Panels A & C)", fontsize=8.0, fontweight="bold", color="#1e293b", labelpad=1.8)
        cbar_a.ax.tick_params(labelsize=7.0, pad=1)

    # Lock identical X-limits across ax_a, ax_a_top, and ax_a_top_r so bars & dots align pixel-for-pixel
    ax_a.set_xlim(x_min_a, x_max_a)
    ax_a_top.set_xlim(x_min_a, x_max_a)
    ax_a_top_r.set_xlim(x_min_a, x_max_a)

    ax_a.set_xlabel("Bin Classifiable Region Length (bp, sorted shortest → longest contig)", fontsize=9.0, color="#1e293b", labelpad=3)
    ax_a.set_ylabel("Mapped eDNA Read Pairs", fontsize=9.2, color="#1e293b")
    ax_a.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax_a.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    y_max_a = float(np.max(np.maximum(y_all, y_p75))) if len(y_all) > 0 else 1.0
    ax_a.set_ylim(bottom=min(-0.03 * y_max_a, 0.0), top=y_max_a * 1.22)
    ax_a.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.84, fontsize=7.5, ncol=2)

    # ======================================================================== #
    # PANEL B: Outlier Pruning Trajectory (Row 0, Col 1 — Same Width as Panel F)
    # ======================================================================== #
    ax_b = fig.add_subplot(master_gs[0, 1])
    iters = [int(r["iteration"]) for r in traj_rows]
    r_trace = [float(r["pearson_r"]) for r in traj_rows]
    cr_pct_trace = [float(r["retained_cr_pct"]) for r in traj_rows]

    ax_b.grid(True, alpha=0.7, zorder=0)
    ax_b_right = ax_b.twinx()

    l_r = ax_b.plot(
        iters, r_trace, color="#2563eb", marker="o", markersize=4.5, linewidth=2.0,
        label="Pearson R", zorder=4
    )
    t_r = ax_b.axhline(
        target_r, color="#16a34a", linestyle="--", linewidth=1.2,
        label=f"Target R ≥ {target_r:.2f}", zorder=2
    )
    l_cr = ax_b_right.plot(
        iters, cr_pct_trace, color="#d97706", marker="s", markersize=4.0, linewidth=1.8,
        label="Retained CR (%)", zorder=3
    )
    t_cr = ax_b_right.axhline(
        100.0 * min_retained_cr_frac, color="#dc2626", linestyle=":", linewidth=1.2,
        label=f"Min CR ({100.0 * min_retained_cr_frac:.0f}%)", zorder=2
    )

    set_panel_header(
        ax_b,
        "B. Outlier Pruning Trajectory",
        "Correlation gain vs. CR sequence retained per step",
    )
    ax_b.set_xlabel("Outlier Bins Removed (Step #)", fontsize=9.0, color="#1e293b")
    ax_b.set_ylabel("Linear Correlation (Pearson R)", fontsize=9.0, color="#2563eb")
    ax_b_right.set_ylabel("Retained CR Sequence (%)", fontsize=9.0, color="#d97706")
    ax_b.set_ylim(min(0.0, min(r_trace) - 0.05) if r_trace else 0.0, 1.04)
    y_min_cr_b = max(0.0, min(cr_pct_trace) - 10.0) if cr_pct_trace else 0.0
    ax_b_right.set_ylim(y_min_cr_b, 103.5)
    ax_b_right.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, integer=True))
    ticks_cr_b = [t for t in ax_b_right.get_yticks() if y_min_cr_b <= t <= 100.0]
    if not ticks_cr_b or ticks_cr_b[-1] < 100.0:
        ticks_cr_b.append(100.0)
    ax_b_right.set_yticks(sorted(set(ticks_cr_b)))
    ax_b_right.set_ylim(y_min_cr_b, 103.5)
    if iters and max(iters) == 0:
        ax_b.set_xlim(-0.4, 1.0)
        ax_b.set_xticks([0])
    else:
        ax_b.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    lines_b = l_r + [t_r] + l_cr + [t_cr]
    ax_b.legend(lines_b, [l.get_label() for l in lines_b], loc="lower left", frameon=True, facecolor="white", framealpha=0.88, fontsize=7.4)

    # ======================================================================== #
    # PANEL D: Outlier Classification Quadrant (Row 0, Col 2 — Same Width as Panel I)
    # ======================================================================== #
    ax_d = fig.add_subplot(master_gs[0, 2])
    ax_d.grid(True, alpha=0.7, zorder=0)

    log2_oe = np.array([float(r["log2_obs_exp"]) for r in bin_rows], dtype=np.float64)
    gini_all = np.array([float(r.get("gini_coeff", 0.0)) for r in bin_rows], dtype=np.float64)
    occ_pct_all = np.array([float(r.get("sample_occupancy_pct", 0.0)) for r in bin_rows], dtype=np.float64)
    sizes_d = 26.0 + 340.0 * (np.clip(gini_all, 0.0, 1.0) ** 1.45)

    ax_d.axvline(0.0, color="#64748b", linestyle="--", linewidth=1.1, zorder=1)
    ax_d.axvline(1.0, color="#94a3b8", linestyle=":", linewidth=1.0, zorder=1)
    ax_d.axhline(0.60, color="#94a3b8", linestyle=":", linewidth=1.0, zorder=1)

    x_max_d = max(2.4, float(np.max(log2_oe)) * 1.18) if len(log2_oe) > 0 else 2.5
    x_min_d = min(-1.5, float(np.min(log2_oe)) * 1.15) if len(log2_oe) > 0 else -1.5
    ax_d.text(
        x_max_d - 0.08, 0.74, "Endogenous Multi-Copy Repeat\n(High depth, matches species phenology)",
        ha="right", va="top", fontsize=7.0, fontweight="bold", color="#047857",
        bbox=dict(boxstyle="round,pad=0.22", fc="#ecfdf5", ec="#a7f3d0", alpha=0.88), zorder=2
    )
    ax_d.text(
        x_max_d - 0.08, 0.04, "Foreign Cross-Mapping Artifact\n(High depth, discordant time series)",
        ha="right", va="bottom", fontsize=7.0, fontweight="bold", color="#be123c",
        bbox=dict(boxstyle="round,pad=0.22", fc="#fff1f2", ec="#fecdd3", alpha=0.88), zorder=2
    )

    ax_d.scatter(
        log2_oe[is_ret], cos_all[is_ret], s=sizes_d[is_ret],
        c=palette["primary"], edgecolors="white", linewidths=0.6, alpha=0.82,
        marker="o", zorder=3
    )
    if np.any(~is_ret):
        ax_d.scatter(
            log2_oe[~is_ret], cos_all[~is_ret], s=sizes_d[~is_ret],
            c="#e11d48", edgecolors="#881337", linewidths=0.9, alpha=0.90,
            marker="^", zorder=4
        )
        pruned_d_list = sorted([i for i, f in enumerate(is_ret) if not f], key=lambda idx: pruned_steps[idx])[:6]
        x_mid_d = x_min_d + 0.58 * (x_max_d - x_min_d)
        dy_cycle = [6, -10, 13, -16, 8, -12]
        for rank_i, idx in enumerate(pruned_d_list):
            dy_d = dy_cycle[rank_i % len(dy_cycle)]
            is_right_side = log2_oe[idx] > x_mid_d
            dx_d = -5 if is_right_side else 5
            ha_d = "right" if is_right_side else "left"
            ax_d.annotate(
                f"#{pruned_steps[idx]} ({bin_rows[idx]['bin_id']})",
                xy=(log2_oe[idx], cos_all[idx]),
                xytext=(dx_d, dy_d), textcoords="offset points",
                ha=ha_d, fontsize=6.8, fontweight="bold", color="#9f1239", zorder=5
            )

    legend_handles_d = [
        plt.scatter([], [], s=50, c=palette["primary"], marker="o", edgecolors="white", label="Retained Bin"),
        plt.scatter([], [], s=70, c="#e11d48", marker="^", edgecolors="#881337", label="Pruned Bin"),
        plt.scatter([], [], s=26.0 + 340.0 * (0.20 ** 1.45), c="#94a3b8", marker="o", alpha=0.75, label="Gini = 0.20 (Uniform CRs)"),
        plt.scatter([], [], s=26.0 + 340.0 * (0.55 ** 1.45), c="#64748b", marker="o", alpha=0.75, label="Gini = 0.55 (Moderate)"),
        plt.scatter([], [], s=26.0 + 340.0 * (0.85 ** 1.45), c="#334155", marker="o", alpha=0.75, label="Gini ≥ 0.85 (CR Hotspot)"),
    ]
    ax_d.set_xlim(x_min_d, x_max_d)
    ax_d.set_ylim(-0.04, 1.15)
    set_panel_header(
        ax_d,
        "D. Outlier Classification Quadrant",
        "Excess read depth vs. temporal concordance (dot size = Gini inequality)",
    )
    ax_d.set_xlabel("Excess Read Depth [log₂(Observed / Expected)]", fontsize=9.0, color="#1e293b")
    ax_d.set_ylabel("Temporal Concordance (Cosine Sim)", fontsize=9.0, color="#1e293b")
    ax_d.legend(handles=legend_handles_d, loc="lower left", frameon=True, facecolor="white", framealpha=0.88, fontsize=6.9)

    # ======================================================================== #
    # PANEL E1, E2 & E3 (Row 0, Col 3 — Same Width as Panel J):
    #   - Left (E1): eDNA Read Attrition Funnel (5 stages)
    #   - Top-Right (E2, ~half height): Alignment Identity Spectrum ([90–100%])
    #   - Bottom-Right (E3, ~half height, under E2): Reference Sequence Funnel
    # ======================================================================== #
    gs_e = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=master_gs[0, 3], width_ratios=[0.94, 1.06], wspace=0.20
    )
    ax_e1 = fig.add_subplot(gs_e[0])

    gs_e_right = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_e[1], height_ratios=[1.05, 0.95], hspace=0.52
    )
    ax_e2 = fig.add_subplot(gs_e_right[0])
    ax_e3 = fig.add_subplot(gs_e_right[1])

    # E1: True Centered Symmetrical Trapezoidal Funnel Plot (Reads)
    stages = funnel["stage_labels"]
    counts = np.array(funnel["stage_counts"], dtype=np.float64)
    max_c = max(1.0, float(counts[0]))
    widths = 0.46 + 0.50 * np.clip(counts / max_c, 0.0, 1.0)
    n_st = len(stages)
    y_centers = np.arange(n_st - 1, -1, -1, dtype=np.float64)
    band_h = 0.52
    stage_colors = ["#64748b", "#0284c7", "#2563eb", "#0d9488", palette["primary"]]

    for i in range(n_st):
        yc = y_centers[i]
        w = widths[i]
        if i < n_st - 1:
            yc_next = y_centers[i + 1]
            w_next = widths[i + 1]
            trap_verts = [
                (-w / 2.0, yc - band_h / 2.0),
                (w / 2.0, yc - band_h / 2.0),
                (w_next / 2.0, yc_next + band_h / 2.0),
                (-w_next / 2.0, yc_next + band_h / 2.0),
            ]
            ax_e1.add_patch(Polygon(trap_verts, closed=True, facecolor=stage_colors[i], alpha=0.22, edgecolor="none", zorder=2))

        rect_verts = [
            (-w / 2.0, yc - band_h / 2.0),
            (w / 2.0, yc - band_h / 2.0),
            (w / 2.0, yc + band_h / 2.0),
            (-w / 2.0, yc + band_h / 2.0),
        ]
        ax_e1.add_patch(Polygon(rect_verts, closed=True, facecolor=stage_colors[i], alpha=0.92, edgecolor="white", linewidth=1.0, zorder=3))

        pct_tot = 100.0 * counts[i] / max_c
        ax_e1.text(
            0.0, yc + 0.08, stages[i],
            ha="center", va="center", fontsize=7.2, fontweight="bold",
            color="white" if w > 0.42 else "#0f172a",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75) if w <= 0.42 else None,
            zorder=5
        )
        ax_e1.text(
            0.0, yc - 0.13, f"{fmt_compact(counts[i])} ({pct_tot:.1f}%)",
            ha="center", va="center", fontsize=7.1, fontweight="bold",
            color="#f8fafc" if w > 0.42 else "#1e293b",
            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.75) if w <= 0.42 else None,
            zorder=5
        )

    ax_e1.set_xlim(-0.56, 0.56)
    ax_e1.set_ylim(-0.5, n_st - 0.4)
    ax_e1.axis("off")
    set_panel_header(
        ax_e1,
        "E1. eDNA Read Attrition",
        "Read pairs retained per filter stage",
    )

    # E2: [90%, 100%] Sequence Identity Histogram for CR Reads (cut to ~half height)
    ax_e2.grid(True, axis="y", alpha=0.6, zorder=0)
    lt99_map: Dict[int, float] = funnel["lt99_hist_pairs"]
    x_bins_e2 = list(range(89, 101))  # 89 (<90%), 90..98, 99, 100
    x_labels_e2 = ["<90", "90", "91", "92", "93", "94", "95", "96", "97", "98", "99", "100"]

    lt99_vals = [lt99_map.get(b, 0.0) for b in range(89, 99)] + [0.0, 0.0]
    ret_pass_vals = [0.0] * 10 + [funnel["ret_1mm_99_pairs"], funnel["ret_exact_100_pairs"]]
    prn_pass_vals = [0.0] * 10 + [funnel["pruned_1mm_99_pairs"], funnel["pruned_exact_100_pairs"]]

    idx_e2 = np.arange(len(x_bins_e2))
    ax_e2.bar(
        idx_e2[:10], lt99_vals[:10], color="#e11d48", width=0.74, alpha=0.88,
        label=f"Rejected <99% ({fmt_compact(funnel['rejected_lt99_pairs'])})", zorder=3
    )
    ax_e2.bar(
        idx_e2[10:], ret_pass_vals[10:], color=palette["primary"], width=0.74, alpha=0.92,
        label=f"Clean ≥99% ({fmt_compact(sum(ret_pass_vals))})", zorder=3
    )
    ax_e2.bar(
        idx_e2[10:], prn_pass_vals[10:], bottom=ret_pass_vals[10:], color="#f59e0b", width=0.74, alpha=0.90,
        label=f"Pruned ≥99% ({fmt_compact(sum(prn_pass_vals))})", zorder=3
    )
    ax_e2.axvline(9.5, color="#0f172a", linestyle="--", linewidth=1.1, zorder=4)
    ax_e2.set_xticks(idx_e2)
    ax_e2.set_xticklabels(x_labels_e2, fontsize=6.5, rotation=0)
    ax_e2.tick_params(axis="x", pad=1.5)
    ax_e2.tick_params(axis="y", labelsize=6.8)
    ax_e2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt_compact(v)))
    ax_e2.set_xlabel("Sequence Identity (%)", fontsize=8.0, color="#1e293b", labelpad=2)
    ax_e2.set_ylabel("Read Pairs", fontsize=8.0, color="#1e293b")
    set_panel_header(
        ax_e2,
        "E2. Alignment Identity Spectrum",
        "Read pairs by % identity (≥99% cutoff)",
        title_size=11.2,
        sub_size=8.3,
    )
    ax_e2.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.86, fontsize=6.2, borderpad=0.25, labelspacing=0.2)

    # E3: Reference Sequence Funnel (Total Assembly Size -> Total CR Size -> Total CR Size After Pruning)
    asm_counts = np.array([float(total_asm_bp), float(total_cr_bp), float(retained_cr_bp)], dtype=np.float64)
    max_asm = max(1.0, float(asm_counts[0]))
    widths_e3 = 0.46 + 0.50 * np.clip(asm_counts / max_asm, 0.0, 1.0)
    stages_e3 = [
        "1. Total Genome Assembly",
        "2. Classifiable Regions (CR)",
        "3. Clean CRs (Post-Pruning)",
    ]
    cr_pct_asm = 100.0 * asm_counts[1] / max_asm
    ret_pct_asm = 100.0 * asm_counts[2] / max_asm
    ret_pct_cr = 100.0 * asm_counts[2] / max(1.0, float(asm_counts[1]))
    sublabels_e3 = [
        f"{fmt_bp(asm_counts[0])} (100.0%)",
        f"{fmt_bp(asm_counts[1])} ({cr_pct_asm:.1f}% genome)",
        f"{fmt_bp(asm_counts[2])} ({ret_pct_asm:.1f}% genome | {ret_pct_cr:.1f}% CR)",
    ]
    colors_e3 = ["#475569", "#0284c7", palette["primary"]]
    y_centers_e3 = np.array([2.0, 1.0, 0.0], dtype=np.float64)
    band_h_e3 = 0.56

    for i in range(3):
        yc = y_centers_e3[i]
        w = widths_e3[i]
        if i < 2:
            yc_next = y_centers_e3[i + 1]
            w_next = widths_e3[i + 1]
            trap_verts = [
                (-w / 2.0, yc - band_h_e3 / 2.0),
                (w / 2.0, yc - band_h_e3 / 2.0),
                (w_next / 2.0, yc_next + band_h_e3 / 2.0),
                (-w_next / 2.0, yc_next + band_h_e3 / 2.0),
            ]
            ax_e3.add_patch(Polygon(trap_verts, closed=True, facecolor=colors_e3[i], alpha=0.22, edgecolor="none", zorder=2))

        rect_verts = [
            (-w / 2.0, yc - band_h_e3 / 2.0),
            (w / 2.0, yc - band_h_e3 / 2.0),
            (w / 2.0, yc + band_h_e3 / 2.0),
            (-w / 2.0, yc + band_h_e3 / 2.0),
        ]
        ax_e3.add_patch(Polygon(rect_verts, closed=True, facecolor=colors_e3[i], alpha=0.92, edgecolor="white", linewidth=1.0, zorder=3))

        ax_e3.text(
            0.0, yc + 0.09, stages_e3[i],
            ha="center", va="center", fontsize=6.9, fontweight="bold",
            color="white" if w > 0.48 else "#0f172a",
            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.78) if w <= 0.48 else None,
            zorder=5
        )
        ax_e3.text(
            0.0, yc - 0.14, sublabels_e3[i],
            ha="center", va="center", fontsize=6.7, fontweight="bold",
            color="#f8fafc" if w > 0.48 else "#1e293b",
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.78) if w <= 0.48 else None,
            zorder=5
        )

    ax_e3.set_xlim(-0.56, 0.56)
    ax_e3.set_ylim(-0.48, 2.52)
    ax_e3.axis("off")
    set_panel_header(
        ax_e3,
        "E3. Reference Sequence Funnel",
        "Assembly → CR → clean CR sequence",
        title_pts=14.5,
        sub_pts=2.5,
        title_size=11.2,
        sub_size=8.3,
    )

    # ======================================================================== #
    # PANEL C: Genome-Wide CR Coverage Map / Manhattan Plot (Row 1, Col 0 — Same Width as A & G)
    # ======================================================================== #
    ax_c = fig.add_subplot(master_gs[1, 0])
    ax_c.grid(True, alpha=0.6, zorder=0)

    # Visually indicate the linear [0, 1] band vs. log10 (> 1) region on the symlog Y-axis with stronger shading
    ax_c.axhspan(0.0, 1.0, color="#cbd5e1", alpha=0.55, zorder=0.5)
    ax_c.axhline(1.0, color="#475569", linestyle="--", linewidth=1.1, alpha=0.90, zorder=2.5)

    for r in bin_rows:
        if r["is_retained"].lower() != "true":
            ax_c.axvspan(
                float(r["bin_start_cum_mb"]), float(r["bin_end_cum_mb"]),
                color="#ffe4e6", alpha=0.82, zorder=1
            )

    cum_mb = cr_full["cum_mid_bp"].to_numpy() / 1e6
    rpkb = cr_full["reads_per_kb"].to_numpy()
    cr_cos = cr_full["cos_sim_516w"].to_numpy()

    max_mb_c = float(np.max(cum_mb)) if len(cum_mb) > 0 else 1.0
    pad_mb_c = max(0.02, 0.028 * max_mb_c)

    pos_mask = rpkb > 0
    zero_idx = np.where(~pos_mask)[0]
    if len(zero_idx) > 8000:
        zero_idx = zero_idx[:: max(1, len(zero_idx) // 8000)]
    ax_c.scatter(
        cum_mb[zero_idx], np.zeros_like(zero_idx, dtype=np.float64),
        s=3.5, c="#94a3b8", alpha=0.35, rasterized=True, zorder=2
    )

    if np.any(pos_mask):
        pos_idx = np.where(pos_mask)[0]
        order_pos = pos_idx[np.argsort(rpkb[pos_idx])]
        ax_c.scatter(
            cum_mb[order_pos], rpkb[order_pos],
            s=10.0,
            c=cr_cos[order_pos], cmap=cmap_cos, norm=norm_cos,
            alpha=0.88, edgecolors="none", rasterized=True, zorder=4
        )

    ax_c.set_yscale("symlog", linthresh=1.0)
    # Give equal horizontal padding on both the left and right edges so first/last CR dots are never cut in half
    ax_c.set_xlim(-pad_mb_c, max_mb_c + pad_mb_c)
    ax_c.set_ylim(bottom=-0.15, top=max(10.0, float(np.max(rpkb)) * 3.0 if len(rpkb) > 0 else 10.0))

    # Indicate the y=1 linear-to-log10 transition directly on the plot
    ax_c.text(
        max_mb_c + pad_mb_c * 0.85, 1.0,
        " y=1 (Linear 0–1 | Log₁₀ >1) ",
        ha="right", va="bottom", fontsize=7.0, fontweight="bold", color="#334155",
        bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="#94a3b8", alpha=0.90),
        zorder=5,
    )

    set_panel_header(
        ax_c,
        "C. Genome-Wide CR Coverage Map (Manhattan Plot)",
        f"Read density across {cr_full.height:,} CR segments (rose bands = pruned bins; gray band = linear 0–1 scale)",
    )
    ax_c.set_xlabel("Cumulative Classifiable Region Coordinate (Mb, bin_001 → bin_100)", fontsize=9.0, color="#1e293b")
    ax_c.set_ylabel("Read Density (Reads / kb)", fontsize=9.2, color="#1e293b")

    # ======================================================================== #
    # PANEL G: 10-Year Longitudinal eDNA Time Series (Row 2, Col 0 — Same Width as A & C)
    # ======================================================================== #
    ax_g = fig.add_subplot(master_gs[2, 0])
    ax_g.grid(True, alpha=0.6, zorder=0)

    w_idx = np.arange(len(weekly_rows))
    years_seq = [int(r["year"]) for r in weekly_rows]
    year_ticks: List[int] = []
    year_labels: List[str] = []
    for y_i, yr in enumerate(sorted(set(years_seq))):
        pos = [i for i, val in enumerate(years_seq) if val == yr]
        if pos:
            year_ticks.append(pos[0])
            year_labels.append(str(yr))
            if y_i % 2 == 1:
                ax_g.axvspan(pos[0] - 0.5, pos[-1] + 0.5, color="#f8fafc", zorder=0)

    ax_g.fill_between(
        w_idx, 0, raw_w, color="#94a3b8", alpha=0.30,
        label=f"Raw Kraken Clade Pairs (sum = {int(np.sum(raw_w)):,})", zorder=1
    )
    ax_g.plot(w_idx, raw_w, color="#64748b", linewidth=0.85, alpha=0.75, zorder=2)
    ax_g.plot(
        w_idx, all_cr_w, color="#0284c7", linewidth=1.05, alpha=0.85,
        label=f"All-CR Mapped Pairs (sum = {int(np.sum(all_cr_w)):,})", zorder=3
    )
    ax_g.plot(
        w_idx, clean_w, color=palette["primary"], linewidth=1.55,
        label=f"Validated Clean-Bin Pairs (sum = {int(np.sum(clean_w)):,})", zorder=4
    )

    set_panel_header(
        ax_g,
        "G. 10-Year Longitudinal eDNA Time Series (2013–2022)",
        "Weekly read pairs across 516 air-filter samples before and after validation",
    )
    ax_g.set_xlabel("Monitoring Year (Weekly Air-Filter Samples, Ljungbyhed)", fontsize=9.0, color="#1e293b")
    ax_g.set_ylabel("Weekly Read Pairs", fontsize=9.2, color="#1e293b")
    if year_ticks:
        ax_g.set_xticks(year_ticks)
        ax_g.set_xticklabels(year_labels, fontsize=8.8)
    ax_g.set_xlim(0, max(1, len(weekly_rows) - 1))
    ax_g.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt_compact(v)))
    ax_g.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.90, fontsize=8.0)

    # ======================================================================== #
    # PANEL F: Weekly Signal Retention (Row 1, Col 1 — Same Width as Panel B)
    # ======================================================================== #
    ax_f = fig.add_subplot(master_gs[1, 1])
    ax_f.grid(True, alpha=0.7, zorder=0)

    log_raw_w = np.log10(raw_w + 1.0)
    log_clean_w = np.log10(clean_w + 1.0)
    r_wk_log = float(stats.pearsonr(log_raw_w, log_clean_w)[0]) if (np.std(log_raw_w) > 0 and np.std(log_clean_w) > 0) else 0.0
    rho_wk = float(stats.spearmanr(raw_w, clean_w).statistic) if (np.std(raw_w) > 0 and np.std(clean_w) > 0) else 0.0

    max_log_f = max(2.0, float(np.max(log_raw_w)), float(np.max(log_clean_w))) * 1.08
    ax_f.plot([0, max_log_f], [0, max_log_f], linestyle="--", linewidth=1.3, color="#94a3b8", label="1:1 Identity (Y = X)", zorder=2)

    if np.sum(raw_w) > 0:
        ret_ratio = float(np.sum(clean_w)) / float(np.sum(raw_w))
        x_curve = np.linspace(0, max_log_f, 200)
        y_curve = np.log10(ret_ratio * (10.0 ** x_curve - 1.0) + 1.0)
        ax_f.plot(
            x_curve, y_curve, color="#0f172a", linestyle="-", linewidth=1.6,
            label=f"Y = {ret_ratio:.2f}·X (R_log={r_wk_log:.3f}, ρ={rho_wk:.3f})", zorder=3
        )

    ax_f.scatter(
        log_raw_w, log_clean_w, s=28, c=palette["primary"], edgecolors="white", linewidths=0.45,
        alpha=0.78, label=f"516 Wks (Lin R={r_wk_lin:.3f})", zorder=4
    )
    ax_f.set_xlim(-0.08, max_log_f)
    ax_f.set_ylim(-0.08, max_log_f)
    tick_powers = list(range(0, int(math.ceil(max_log_f)) + 1))
    tick_Text = ["0" if p == 0 else f"10^{p}" for p in tick_powers]
    ax_f.set_xticks(tick_powers)
    ax_f.set_xticklabels(tick_Text, fontsize=7.6)
    ax_f.set_yticks(tick_powers)
    ax_f.set_yticklabels(tick_Text, fontsize=7.6)
    set_panel_header(
        ax_f,
        "F. Weekly Signal Retention",
        "Clean vs. raw counts across 516 weeks",
    )
    ax_f.set_xlabel("Raw Kraken Pairs [log₁₀(x+1)]", fontsize=8.8, color="#1e293b")
    ax_f.set_ylabel("Clean Validated Pairs [log₁₀(y+1)]", fontsize=8.8, color="#1e293b")
    ax_f.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.88, fontsize=6.8)

    # ======================================================================== #
    # PANEL I: Bin Temporal Occupancy (Row 1, Col 2 — Same Width as Panel D)
    # ======================================================================== #
    gs_i = gridspec.GridSpecFromSubplotSpec(
        2, 1,
        subplot_spec=master_gs[1, 2],
        height_ratios=[0.24, 0.76],
        hspace=0.07,
    )
    ax_i_top = fig.add_subplot(gs_i[0])
    ax_i = fig.add_subplot(gs_i[1], sharex=ax_i_top)

    # Effective Number of Active Weeks (Inverse Simpson N_eff) for the clean time series
    sum_clean_w = float(np.sum(clean_w))
    sum_sq_clean_w = float(np.sum(clean_w ** 2))
    n_eff_weeks = (sum_clean_w ** 2) / max(1e-9, sum_sq_clean_w) if sum_clean_w > 0 else 0.0
    med_occ_ret = float(np.median(occ_pct_all[is_ret])) if np.any(is_ret) else 0.0
    med_wks_ret = int(round(med_occ_ret * len(weekly_rows) / 100.0))

    max_occ_x = max(5.0, float(np.max(occ_pct_all)) * 1.12) if len(occ_pct_all) > 0 else 10.0
    bins_occ = np.linspace(0.0, max_occ_x, 22)

    ax_i_top.grid(True, alpha=0.45, zorder=0)
    if np.any(is_ret):
        ax_i_top.hist(occ_pct_all[is_ret], bins=bins_occ, color=palette["primary"], alpha=0.65, edgecolor="white", linewidth=0.5, zorder=2)
    if np.any(~is_ret):
        ax_i_top.hist(occ_pct_all[~is_ret], bins=bins_occ, color="#e11d48", alpha=0.65, edgecolor="white", linewidth=0.5, zorder=3)
    ax_i_top.axvline(med_occ_ret, color="#0f172a", linestyle="--", linewidth=1.1, zorder=4)
    ax_i_top.set_ylabel("Bins", fontsize=7.4, color="#475569")
    ax_i_top.tick_params(axis="x", labelbottom=False, bottom=False)
    ax_i_top.tick_params(axis="y", labelsize=6.8)
    set_panel_header(
        ax_i_top,
        "I. Bin Temporal Occupancy",
        f"Clean median: {med_occ_ret:.1f}% ({med_wks_ret}/516 wks) | Effective active weeks (N_eff) = {n_eff_weeks:.1f}w",
    )

    ax_i.grid(True, alpha=0.7, zorder=0)
    ax_i.axvline(med_occ_ret, color="#0f172a", linestyle="--", linewidth=1.1, label=f"Clean Med ({med_occ_ret:.1f}%, {med_wks_ret}/516w)", zorder=2)
    ax_i.scatter(
        occ_pct_all[is_ret], cos_all[is_ret], s=38, c=palette["primary"],
        edgecolors="white", linewidths=0.5, alpha=0.84, marker="o",
        label=f"Retained (n={int(np.sum(is_ret))})", zorder=3
    )
    if np.any(~is_ret):
        ax_i.scatter(
            occ_pct_all[~is_ret], cos_all[~is_ret], s=58, c="#e11d48",
            edgecolors="#881337", linewidths=0.8, alpha=0.90, marker="^",
            label=f"Pruned (n={int(np.sum(~is_ret))})", zorder=4
        )
    ax_i.set_xlim(-0.03 * max_occ_x, max_occ_x)
    ax_i_top.set_xlim(-0.03 * max_occ_x, max_occ_x)
    ax_i.set_ylim(-0.04, 1.08)
    ax_i.set_xlabel("Weekly Sample Occupancy (% of 516 weeks)", fontsize=8.8, color="#1e293b")
    ax_i.set_ylabel("Temporal Concordance (Cosine Sim)", fontsize=8.8, color="#1e293b")
    ax_i.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.88, fontsize=6.8)

    # ======================================================================== #
    # PANEL H: 10-Year Seasonal Phenology Profile (Row 2, Cols 1..2 — Spans B+D / F+I)
    # ======================================================================== #
    ax_h = fig.add_subplot(master_gs[2, 1:3])
    ax_h.grid(True, alpha=0.65, zorder=0)

    wks = np.array(pheno["weeks"], dtype=np.float64)
    raw_share = np.array(pheno["raw_share"], dtype=np.float64)
    clean_share = np.array(pheno["clean_share"], dtype=np.float64)
    clean_mean_reads = np.array(pheno["clean_mean_reads"], dtype=np.float64)

    ax_h.plot(
        wks, raw_share, color="#475569", linestyle="--", linewidth=1.7,
        label="Pre-Validation Raw (3w Circular Smooth)", zorder=3
    )
    ax_h.plot(
        wks, clean_share, color=palette["primary"], linestyle="-", linewidth=2.3,
        label="Post-Validation Clean (3w Circular Smooth)", zorder=4
    )
    ax_h.fill_between(
        wks, 0, clean_share, color=palette["fill"], alpha=0.40, zorder=2
    )
    ax_h.fill_between(
        wks, clean_share, np.maximum(raw_share, clean_share),
        where=(raw_share > clean_share), color="#fda4af", alpha=0.55,
        label="Purged Off-Season / Cross-Mapping Noise", zorder=2.5
    )

    # Right Y-axis: Mean Weekly Clean Read Pairs (10-year average reads/week)
    ax_h_right = ax_h.twinx()
    max_share = max(1.0, float(np.max(np.maximum(raw_share, clean_share))) * 1.18)
    scale_reads_per_pct = float(np.max(clean_mean_reads)) / max(1e-6, float(np.max(clean_share))) if np.max(clean_share) > 0 else 1.0
    ax_h.set_ylim(0, max_share)
    ax_h_right.set_ylim(0, max_share * scale_reads_per_pct)
    ax_h_right.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt_compact(v)))
    ax_h_right.set_ylabel("Mean Clean Pairs / Week", fontsize=8.8, color="#047857")
    ax_h_right.tick_params(axis="y", labelsize=7.6, colors="#047857")

    ax_h.set_xlim(1, 52)
    month_starts = [1, 5, 9, 14, 18, 22, 27, 31, 35, 40, 44, 48]
    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    ax_h.set_xticks(month_starts)
    ax_h.set_xticklabels(month_names, fontsize=8.4)
    set_panel_header(
        ax_h,
        "H. 10-Year Seasonal Phenology Profile",
        f"Multi-year average cycle across weeks 1–52 (pre- vs. post-validation concordance r = {pheno['r_pheno']:.3f})",
    )
    ax_h.set_xlabel("Calendar Month (Weeks 1–52, 3-Week Circular Moving Average Across Dec–Jan)", fontsize=9.0, color="#1e293b")
    ax_h.set_ylabel("Share of Annual Signal (% / week)", fontsize=9.0, color="#1e293b")
    ax_h.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.88, fontsize=7.5)

    # ======================================================================== #
    # PANEL J: Bin x Week Hierarchical Clustermap (Rows 1..2, Col 3 — Same Width as Panel E)
    #          with Strip Labels at Bottom so Title is Left-Adjusted Like All Other Panels!
    # ======================================================================== #
    gs_j_outer = gridspec.GridSpecFromSubplotSpec(
        2, 1,
        subplot_spec=master_gs[1:3, 3],
        height_ratios=[0.980, 0.020],
        hspace=0.165,
    )
    gs_j_top = gridspec.GridSpecFromSubplotSpec(
        1, 5,
        subplot_spec=gs_j_outer[0],
        width_ratios=[0.16, 0.042, 0.042, 0.042, 0.714],
        wspace=0.035,
    )
    gs_j_bot = gridspec.GridSpecFromSubplotSpec(
        1, 2,
        subplot_spec=gs_j_outer[1],
        width_ratios=[0.38, 0.62],
        wspace=0.22,
    )

    ax_dendro = fig.add_subplot(gs_j_top[0, 0])
    ax_strip_cos = fig.add_subplot(gs_j_top[0, 1])
    ax_strip_occ = fig.add_subplot(gs_j_top[0, 2])
    ax_strip_prn = fig.add_subplot(gs_j_top[0, 3])
    ax_hm = fig.add_subplot(gs_j_top[0, 4])
    ax_occ_cbar = fig.add_subplot(gs_j_bot[0, 0])
    ax_hm_cbar = fig.add_subplot(gs_j_bot[0, 1])

    # Normalize each bin's weekly counts to Reads / kb of CR in that bin, then log10(Reads/kb + 1)
    bin_kb = np.maximum(0.001, x_all / 1000.0)[:, None]
    bin_rpkb_ts = bin_ts_matrix / bin_kb

    # Select top active weeks (up to 65 weeks with highest total reads, sorted chronologically)
    wk_totals = np.sum(bin_ts_matrix, axis=0)
    pos_wks = np.where(wk_totals > 0)[0]
    if len(pos_wks) > 65:
        top_w = pos_wks[np.argsort(wk_totals[pos_wks])[::-1][:65]]
        active_w_idx = np.sort(top_w)
    elif len(pos_wks) >= 5:
        active_w_idx = pos_wks
    else:
        active_w_idx = np.arange(min(52, bin_ts_matrix.shape[1]))

    hm_raw = np.log10(bin_rpkb_ts[:, active_w_idx] + 1.0)

    # Hierarchical clustering of the 100 bins using Ward minimum-variance linkage on Euclidean
    # distance of the log10(Reads/kb + 1) profiles, with Bar-Joseph optimal leaf ordering
    if hm_raw.shape[0] >= 3 and np.std(hm_raw) > 0:
        Z = linkage(hm_raw, method="ward", metric="euclidean", optimal_ordering=True)
        dendro_info = dendrogram(
            Z, orientation="left", ax=ax_dendro, no_labels=True,
            color_threshold=0, above_threshold_color="#475569"
        )
        leaf_order = dendro_info["leaves"][::-1]
    else:
        leaf_order = list(range(len(bin_rows)))

    ax_dendro.invert_yaxis()
    ax_dendro.axis("off")

    hm_ordered = hm_raw[leaf_order, :]
    cos_ordered = cos_all[leaf_order]
    occ_ordered = occ_pct_all[leaf_order]
    ret_ordered = is_ret[leaf_order]

    # Strip 1: Cosine Similarity (Viridis 0..1) — Label at BOTTOM
    ax_strip_cos.imshow(
        cos_ordered[:, None], aspect="auto", cmap=cmap_cos, norm=norm_cos, interpolation="nearest"
    )
    ax_strip_cos.xaxis.tick_bottom()
    ax_strip_cos.set_xticks([0])
    ax_strip_cos.set_xticklabels(["CosSim"], fontsize=6.8, fontweight="bold", rotation=40, ha="right")
    ax_strip_cos.tick_params(axis="x", pad=1.5)
    ax_strip_cos.set_yticks([])

    # Strip 2: Sample Occupancy % across 516 weeks (PuBu 0..max_occ) — Label at BOTTOM
    norm_occ_strip = mcolors.Normalize(vmin=0.0, vmax=max(1.0, float(np.max(occ_pct_all))))
    im_occ = ax_strip_occ.imshow(
        occ_ordered[:, None], aspect="auto", cmap="PuBu", norm=norm_occ_strip, interpolation="nearest"
    )
    ax_strip_occ.xaxis.tick_bottom()
    ax_strip_occ.set_xticks([0])
    ax_strip_occ.set_xticklabels(["Occ%"], fontsize=6.8, fontweight="bold", rotation=40, ha="right")
    ax_strip_occ.tick_params(axis="x", pad=1.5)
    ax_strip_occ.set_yticks([])

    # Strip 3: Pruned Status (Emerald = Retained, Crimson = Pruned) — Label at BOTTOM
    status_rgb = np.array([
        mcolors.to_rgb(palette["primary"]) if flag else mcolors.to_rgb("#e11d48")
        for flag in ret_ordered
    ], dtype=np.float64)[:, None, :]
    ax_strip_prn.imshow(status_rgb, aspect="auto", interpolation="nearest")
    ax_strip_prn.xaxis.tick_bottom()
    ax_strip_prn.set_xticks([0])
    ax_strip_prn.set_xticklabels(["Status"], fontsize=6.8, fontweight="bold", rotation=40, ha="right")
    ax_strip_prn.tick_params(axis="x", pad=1.5)
    ax_strip_prn.set_yticks([])

    # Main Heatmap: 100 Clustered Bins x Active Monitoring Weeks
    im_hm = ax_hm.imshow(
        hm_ordered, aspect="auto", cmap="YlGnBu", interpolation="nearest"
    )
    ax_hm.set_yticks([])
    n_act_w = len(active_w_idx)
    if n_act_w > 0:
        xt_pos = np.linspace(0, n_act_w - 1, min(8, n_act_w), dtype=int)
        xt_lbls = [weekly_rows[active_w_idx[p]]["year_week"] for p in xt_pos]
        ax_hm.set_xticks(xt_pos)
        ax_hm.set_xticklabels(xt_lbls, fontsize=7.0, rotation=35, ha="right")

    # Left-adjusted title & subtitle at the left edge of Panel J (ax_dendro), matching all other panels!
    set_panel_header(
        ax_dendro,
        "J. Bin × Week Hierarchical Clustermap",
        f"100 CR bins clustered by temporal profile across {n_act_w} active weeks (Ward linkage)",
    )

    # Legend for Pruned Status strip at top-right of Panel J
    status_patches = [
        Patch(facecolor=palette["primary"], edgecolor="none", label="Retained"),
        Patch(facecolor="#e11d48", edgecolor="none", label="Pruned"),
    ]
    ax_hm.legend(
        handles=status_patches, loc="upper right", frameon=True,
        facecolor="white", framealpha=0.88, fontsize=7.0, title="Status Strip", title_fontsize=7.0
    )

    # Slim Horizontal Colorbar 1: Occupancy Strip (% of 516 weeks)
    cbar_occ = fig.colorbar(im_occ, cax=ax_occ_cbar, orientation="horizontal")
    cbar_occ.set_label("Occupancy (% of 516 Wks)", fontsize=7.8, color="#1e293b", labelpad=2)
    cbar_occ.ax.tick_params(labelsize=6.9, pad=1)

    # Slim Horizontal Colorbar 2: Heatmap Read Density
    cbar_hm = fig.colorbar(im_hm, cax=ax_hm_cbar, orientation="horizontal")
    cbar_hm.set_label("Read Density [log₁₀(Reads / kb + 1)]", fontsize=7.8, color="#1e293b", labelpad=2)
    cbar_hm.ax.tick_params(labelsize=6.9, pad=1)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate widescreen 10-panel static publication diagnostic dashboards (PNG + PDF)."
    )
    parser.add_argument("--taxid", type=int, default=None, help="Single species TaxID (for Snakemake rule)")
    parser.add_argument("--species-slug", default=None, help="Single species slug (for Snakemake rule)")
    parser.add_argument("--out-png", default=None, help="Explicit output PNG path (for Snakemake rule)")
    parser.add_argument("--out-pdf", default=None, help="Explicit output PDF path (for Snakemake rule)")
    parser.add_argument("--target-r", type=float, default=0.90, help="Target Pearson R threshold")
    parser.add_argument("--min-retained-cr-frac", type=float, default=0.50, help="Min retained CR fraction")
    parser.add_argument("--log", default=None, help="Optional log file path")
    parser.add_argument(
        "--species",
        nargs="+",
        default=[
            "78633:Pieris_napi",
            "270461:Pieris_melete",
            "34373:Blumeria_graminis",
            "30423:Coloeus_monedula",
        ],
        help="List of taxid:species_slug pairs to process in batch mode",
    )
    parser.add_argument("--project-root", default=".", help="Project root directory")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    logger = setup_logger("widescreen_dashboard", args.log)

    if args.taxid is not None and args.species_slug is not None:
        species_list = [(int(args.taxid), str(args.species_slug))]
    else:
        species_list = []
        for item in args.species:
            tid_str, slug = item.split(":", 1)
            species_list.append((int(tid_str), slug))

    for taxid, slug in species_list:
        logger.info(f"=== Processing 10-Panel Diagnostic Suite for {slug} (TaxID: {taxid}) ===")
        data = compute_species_enhanced_data(taxid, slug, project_root, logger)
        val_dir: Path = data["val_dir"]

        if args.out_png and args.out_pdf and len(species_list) == 1:
            out_png = Path(args.out_png)
            out_pdf = Path(args.out_pdf)
            render_10panel_figure(
                data, out_png, out_pdf,
                target_r=args.target_r,
                min_retained_cr_frac=args.min_retained_cr_frac,
            )
        else:
            out_png_10panel = val_dir / f"{taxid}.diagnostic_10panel.png"
            out_pdf_10panel = val_dir / f"{taxid}.diagnostic_10panel.pdf"
            logger.info(f"[{taxid}] Rendering 10-panel figure ({out_png_10panel.name}, {out_pdf_10panel.name})...")
            render_10panel_figure(
                data, out_png_10panel, out_pdf_10panel,
                target_r=args.target_r,
                min_retained_cr_frac=args.min_retained_cr_frac,
            )


if __name__ == "__main__":
    main()
