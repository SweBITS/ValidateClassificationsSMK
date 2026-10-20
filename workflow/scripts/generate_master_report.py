#!/usr/bin/env python3
"""
generate_master_report.py — Stage 5: Master Summary Table, Unified Timeseries Parquet & HTML Dashboard
======================================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
This final aggregation step unifies all 6 stages of `ValidateClassificationsSMK` into three
publication- and downstream-ready deliverables:

1. **`validation_master_summary.tsv` & `validation_master_summary.parquet` (via Polars)**:
   Every target species from `config/targets.tsv` is guaranteed to appear here with its full
   lifecycle metrics across Stage 0 (clade read counts, assembly selection), Stage 1 (>= 1 kb
   contig filtering audit), Stage 2 (2x and 4x Classifiable Region lengths and % of genome),
   Stage 4 (CR read alignment rates), and Stage 5 (initial vs. final Pearson R and Spearman rho,
   outlier bins pruned, % CR and mapped reads retained, and terminal `validation_status` +
   human-readable `status_reason`).

2. **`validated_timeseries_long.parquet`**:
   Unified Polars Parquet table containing the 516-week (2013–2022) longitudinal time series
   for every target taxon (`raw_kraken_reads`, `all_cr_mapped_reads`, `validated_clean_bin_reads`,
   and `validation_status`), ready for immediate ecological modeling in Python or R.

3. **`validation_report_4panel.html` (and `validation_report.html`) & `validation_report_10panel.html`**:
   Self-contained, offline-portable HTML5 visual dashboards with executive summary cards,
   interactive search/filtering across all validated, failed, and excluded species, embedded
   base64 4-panel (`{t_id}.diagnostic_4panel.png`) and 10-panel (`{t_id}.diagnostic_10panel.png`)
   diagnostic figures, pruned-bin inspection tables, and SLURM resource benchmark summaries.
"""

import argparse
import base64
import csv
import datetime
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import jinja2
import polars as pl
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import ValidationStatus, setup_logger


def encode_png_base64(png_path: Path) -> str:
    """Read a PNG file and return a base64 data URI string for self-contained HTML embedding."""
    if not png_path.is_file():
        return ""
    raw = png_path.read_bytes()
    b64 = base64.b64encode(raw).decode("ascii")
    return f"data:image/png;base64,{b64}"


def collect_benchmarks_summary(benchmarks_dir: Path) -> List[Dict[str, Any]]:
    """
    Scan Snakemake `benchmarks/**/*.tsv` files and summarize job count, mean/max walltime (s),
    and peak RSS memory (MB) by pipeline stage.
    """
    if not benchmarks_dir.is_dir():
        return []

    by_stage: Dict[str, List[Dict[str, float]]] = {}
    for tsv_file in sorted(benchmarks_dir.rglob("*.tsv")):
        stage_name = tsv_file.parent.name
        try:
            with open(tsv_file, "r", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
                if rows:
                    r = rows[0]
                    sec = float(r.get("s") or 0.0)
                    rss = float(r.get("max_rss") or 0.0) if r.get("max_rss") not in (None, "-", "NA", "") else 0.0
                    io_in = float(r.get("io_in") or 0.0) if r.get("io_in") not in (None, "-", "NA", "") else 0.0
                    io_out = float(r.get("io_out") or 0.0) if r.get("io_out") not in (None, "-", "NA", "") else 0.0
                    by_stage.setdefault(stage_name, []).append({
                        "s": sec, "max_rss": rss, "io_in": io_in, "io_out": io_out
                    })
        except Exception:
            continue

    summary: List[Dict[str, Any]] = []
    for stage, items in sorted(by_stage.items()):
        n = len(items)
        secs = [x["s"] for x in items]
        rss_vals = [x["max_rss"] for x in items]
        summary.append({
            "stage": stage,
            "job_count": n,
            "total_walltime_min": round(sum(secs) / 60.0, 2),
            "max_walltime_min": round(max(secs) / 60.0, 2),
            "mean_walltime_sec": round(sum(secs) / n, 1),
            "max_rss_mb": round(max(rss_vals), 1),
            "mean_rss_mb": round(sum(rss_vals) / n, 1),
        })
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 5: Generate master validation summary TSV/Parquet, unified 516-week Parquet, and self-contained HTML report."
    )
    parser.add_argument("--wf-info", required=True, help="Path to workflow_run_info.json")
    parser.add_argument("--config-yaml", required=True, help="Path to config/config.yaml")
    parser.add_argument("--stage0-ledger", required=True, help="Path to stage0_preflight_and_assemblies.tsv")
    parser.add_argument("--stage1-audit", required=True, help="Path to assembly_filtering_summary.tsv")
    parser.add_argument("--cr-summary", required=True, help="Path to cr_streaming_summary.tsv")
    parser.add_argument("--raw-weekly-long", required=True, help="Path to raw_weekly_counts_long.tsv")
    parser.add_argument("--validation-dir", required=True, help="Base directory results/04_validation")
    parser.add_argument("--benchmarks-dir", required=True, help="Base directory benchmarks/")
    parser.add_argument("--template", required=True, help="Path to validation_report.html.j2")
    parser.add_argument("--out-master-tsv", required=True, help="Output TSV: validation_master_summary.tsv")
    parser.add_argument("--out-master-parquet", required=True, help="Output Parquet: validation_master_summary.parquet")
    parser.add_argument("--out-timeseries-parquet", required=True, help="Output Parquet: validated_timeseries_long.parquet")
    parser.add_argument("--out-html", required=True, help="Output HTML (4-panel report): validation_report.html")
    parser.add_argument("--out-html-4panel", "--out-html-2x2", dest="out_html_4panel", default=None, help="Explicit output HTML for 4-panel report: validation_report_4panel.html")
    parser.add_argument("--out-html-10panel", "--out-html-widescreen", dest="out_html_10panel", default=None, help="Output HTML for 10-panel report: validation_report_10panel.html")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger("generate_master_report", args.log)

    # 1. Load workflow run info & configuration
    wf_info = json.loads(Path(args.wf_info).read_text(encoding="utf-8")) if Path(args.wf_info).is_file() else {}
    wf_info["end_timestamp"] = datetime.datetime.now().astimezone().isoformat()
    with open(args.config_yaml, "r", encoding="utf-8") as cfh:
        cfg = yaml.safe_load(cfh) or {}

    # 2. Load Stage 0 ledger, Stage 1 audit, and Stage 2 CR summary
    with open(args.stage0_ledger, "r", encoding="utf-8") as fh:
        stage0_rows = list(csv.DictReader(fh, delimiter="\t"))

    stage1_by_tid: Dict[int, Dict[str, str]] = {}
    with open(args.stage1_audit, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            stage1_by_tid[int(row["t_id"])] = row

    cr_by_tid: Dict[int, Dict[str, str]] = {}
    if Path(args.cr_summary).is_file():
        with open(args.cr_summary, "r", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                cr_by_tid[int(row["t_id"])] = row

    val_dir = Path(args.validation_dir)

    # 3. Build Master Summary rows and HTML species cards for EVERY target in config/targets.tsv
    master_rows: List[Dict[str, Any]] = []
    species_cards_4panel: List[Dict[str, Any]] = []
    species_cards_10panel: List[Dict[str, Any]] = []
    weekly_frames: List[pl.DataFrame] = []

    # Load baseline raw weekly long table for fallback on Stage 0 / Stage 1 excluded species
    raw_weekly_df = pl.read_csv(args.raw_weekly_long, separator="\t")

    for s0 in stage0_rows:
        tid = int(s0["t_id"])
        slug = s0["species_slug"]
        s1 = stage1_by_tid.get(tid, {})
        cr = cr_by_tid.get(tid, {})

        sp_val_dir = val_dir / f"{tid}_{slug}"
        metrics_json_path = sp_val_dir / f"{tid}.validation_metrics.json"
        png_4panel_path = sp_val_dir / f"{tid}.diagnostic_4panel.png"
        if not png_4panel_path.is_file():
            png_4panel_path = sp_val_dir / f"{tid}.diagnostic_2x2.png"
        if not png_4panel_path.is_file():
            png_4panel_path = sp_val_dir / f"{tid}.diagnostic_3panel.png"
        png_10panel_path = sp_val_dir / f"{tid}.diagnostic_10panel.png"
        if not png_10panel_path.is_file():
            png_10panel_path = sp_val_dir / f"{tid}.diagnostic_widescreen.png"
        if not png_10panel_path.is_file():
            png_10panel_path = sp_val_dir / f"{tid}.diagnostic_8panel.png"
        if not png_10panel_path.is_file():
            png_10panel_path = png_4panel_path
        traj_path = sp_val_dir / f"{tid}.pruning_trajectory.tsv"
        weekly_path = sp_val_dir / f"{tid}.weekly_timeseries_516w.tsv"

        if metrics_json_path.is_file():
            m = json.loads(metrics_json_path.read_text(encoding="utf-8"))
            final_status = m["validation_status"]
            final_reason = m["status_reason"]
        else:
            m = {}
            final_status = s1.get("stage1_status") or s0.get("stage0_status") or ValidationStatus.EXCLUDED_NO_ASSEMBLY_IN_DB
            final_reason = s1.get("stage1_reason") or s0.get("stage0_reason") or "Excluded during pre-flight checks"

        # Collect 516-week time series for unified Parquet
        if weekly_path.is_file():
            sp_w_df = pl.read_csv(weekly_path, separator="\t")
            weekly_frames.append(sp_w_df)
        else:
            # Construct zero-mapped weekly frame from raw_weekly_df for excluded species
            sub_raw = raw_weekly_df.filter(pl.col("t_id") == tid)
            if len(sub_raw) > 0:
                sp_w_df = sub_raw.with_columns([
                    pl.lit(0).cast(pl.Int64).alias("all_cr_mapped_reads"),
                    pl.lit(0).cast(pl.Int64).alias("validated_clean_bin_reads"),
                    pl.lit(final_status).alias("validation_status"),
                ]).select([
                    "t_id", "t_species", "year_week", "year", "week",
                    "raw_kraken_reads", "all_cr_mapped_reads", "validated_clean_bin_reads",
                    "validation_status",
                ])
                weekly_frames.append(sp_w_df)

        # Parse pruned bins list from trajectory TSV if available
        pruned_steps_list: List[Dict[str, Any]] = []
        if traj_path.is_file():
            with open(traj_path, "r", encoding="utf-8") as tfh:
                for trow in csv.DictReader(tfh, delimiter="\t"):
                    if int(trow["iteration"]) > 0:
                        pruned_steps_list.append(trow)

        row_dict: Dict[str, Any] = {
            "t_id": tid,
            "t_species": s0["t_species"],
            "scientific_name": s0["scientific_name"],
            "common_name": s0.get("common_name", ""),
            "t_domain": s0.get("t_domain", ""),
            "t_phylum": s0.get("t_phylum", ""),
            "t_class": s0.get("t_class", ""),
            "t_order": s0.get("t_order", ""),
            "t_family": s0.get("t_family", ""),
            "t_genus": s0.get("t_genus", ""),
            "clade_size": int(s0.get("clade_size") or 0),
            "total_clade_reads": int(s0.get("total_clade_reads") or 0),
            "weeks_present": int(s0.get("weeks_present") or 0),
            "chosen_accession": s0.get("chosen_accession", "NA"),
            "clean_ncbi_accession": s0.get("clean_ncbi_accession", "NA"),
            "chosen_assembly_taxid": s0.get("chosen_assembly_taxid", "NA"),
            "chosen_level": s0.get("chosen_level", "NA"),
            "chosen_refseq_category": s0.get("chosen_refseq_category", "NA"),
            "selection_reason": s0.get("selection_reason", "NA"),
            "raw_contig_count": int(s1.get("raw_contig_count") or 0),
            "kept_contig_count": int(s1.get("kept_contig_count") or 0),
            "discarded_contig_pct": float(s1.get("discarded_contig_pct") or 0.0),
            "raw_assembly_bp": int(s1.get("raw_assembly_bp") or 0),
            "kept_assembly_bp": int(s1.get("kept_assembly_bp") or 0),
            "discarded_bp": int(s1.get("discarded_bp") or 0),
            "discarded_bp_pct": float(s1.get("discarded_bp_pct") or 0.0),
            "cr_bp_2x": int(cr.get("cr_bp_2x") or 0),
            "cr_pct_of_filtered_assembly_2x": float(cr.get("cr_pct_of_filtered_assembly_2x") or 0.0),
            "cr_segments_2x": int(cr.get("cr_segments_2x") or 0),
            "cr_bp_4x": int(cr.get("cr_bp_4x") or 0),
            "cr_pct_of_filtered_assembly_4x": float(cr.get("cr_pct_of_filtered_assembly_4x") or 0.0),
            "cr_segments_4x": int(cr.get("cr_segments_4x") or 0),
            "n_bins": int(m.get("n_bins") or 0),
            "occupied_bins": int(m.get("occupied_bins") or 0),
            "occupied_bins_frac": float(m.get("occupied_bins_frac") or 0.0),
            "initial_pearson_r": float(m.get("initial_pearson_r") or 0.0),
            "initial_pearson_p": float(m.get("initial_pearson_p") or 1.0),
            "initial_spearman_rho": float(m.get("initial_spearman_rho") or 0.0),
            "initial_ols_slope": float(m.get("initial_ols_slope") or 0.0),
            "final_pearson_r": float(m.get("final_pearson_r") or 0.0),
            "final_pearson_p": float(m.get("final_pearson_p") or 1.0),
            "final_spearman_rho": float(m.get("final_spearman_rho") or 0.0),
            "final_ols_slope": float(m.get("final_ols_slope") or 0.0),
            "bins_pruned_count": int(m.get("bins_pruned_count") or 0),
            "retained_cr_bp": int(m.get("retained_cr_bp") or 0),
            "retained_cr_pct": float(m.get("retained_cr_pct") or 0.0),
            "total_mapped_read_pairs": int(m.get("total_mapped_read_pairs") or 0),
            "retained_mapped_read_pairs": int(m.get("retained_mapped_read_pairs") or 0),
            "retained_mapped_reads_pct": float(m.get("retained_mapped_reads_pct") or 0.0),
            "validation_status": final_status,
            "status_reason": final_reason,
        }
        master_rows.append(row_dict)

        # Category group for UI filtering
        if final_status.startswith("VALIDATED_"):
            status_group = "VALIDATED"
        elif final_status.startswith("FAILED_"):
            status_group = "FAILED"
        else:
            status_group = "EXCLUDED"

        card_4panel = dict(row_dict)
        card_4panel["species_slug"] = slug
        card_4panel["status_group"] = status_group
        card_4panel["plot_data_uri"] = encode_png_base64(png_4panel_path)
        card_4panel["pruned_steps"] = pruned_steps_list
        species_cards_4panel.append(card_4panel)

        card_10panel = dict(row_dict)
        card_10panel["species_slug"] = slug
        card_10panel["status_group"] = status_group
        card_10panel["plot_data_uri"] = encode_png_base64(png_10panel_path)
        card_10panel["pruned_steps"] = pruned_steps_list
        species_cards_10panel.append(card_10panel)

    # 4. Write validation_master_summary.tsv & validation_master_summary.parquet via Polars
    master_df = pl.DataFrame(master_rows)
    out_tsv = Path(args.out_master_tsv)
    out_pq = Path(args.out_master_parquet)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_pq.parent.mkdir(parents=True, exist_ok=True)
    master_df.write_csv(out_tsv, separator="\t")
    master_df.write_parquet(out_pq)

    # 5. Write validated_timeseries_long.parquet
    out_ts_pq = Path(args.out_timeseries_parquet)
    out_ts_pq.parent.mkdir(parents=True, exist_ok=True)
    if weekly_frames:
        all_ts_df = pl.concat(weekly_frames, how="vertical_relaxed")
    else:
        all_ts_df = pl.DataFrame()
    all_ts_df.write_parquet(out_ts_pq)

    # 6. Compute Executive KPI counts & render both 4-panel and 10-panel HTML reports
    kpi = {
        "total_targets": len(master_rows),
        "validated_clean": sum(1 for r in master_rows if r["validation_status"] == ValidationStatus.VALIDATED_CLEAN),
        "validated_pruned": sum(1 for r in master_rows if r["validation_status"] == ValidationStatus.VALIDATED_PRUNED),
        "failed_total": sum(1 for r in master_rows if r["validation_status"].startswith("FAILED_")),
        "excluded_total": sum(1 for r in master_rows if r["validation_status"].startswith("EXCLUDED_")),
    }
    kpi["validated_total"] = kpi["validated_clean"] + kpi["validated_pruned"]

    benchmarks_summary = collect_benchmarks_summary(Path(args.benchmarks_dir))

    template_path = Path(args.template)
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(template_path.parent)),
        autoescape=jinja2.select_autoescape(["html", "xml"]),
    )
    template = env.get_template(template_path.name)
    html_4panel = template.render(
        wf_info=wf_info,
        config=cfg,
        kpi=kpi,
        species_cards=species_cards_4panel,
        benchmarks=benchmarks_summary,
    )

    out_html = Path(args.out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html_4panel, encoding="utf-8")

    out_html_4panel = Path(args.out_html_4panel) if args.out_html_4panel else out_html.parent / "validation_report_4panel.html"
    if out_html_4panel.resolve() != out_html.resolve():
        out_html_4panel.parent.mkdir(parents=True, exist_ok=True)
        out_html_4panel.write_text(html_4panel, encoding="utf-8")

    out_html_10panel = Path(args.out_html_10panel) if args.out_html_10panel else out_html.parent / "validation_report_10panel.html"
    html_10panel = template.render(
        wf_info=wf_info,
        config=cfg,
        kpi=kpi,
        species_cards=species_cards_10panel,
        benchmarks=benchmarks_summary,
    ).replace("max-width: 1520px;", "max-width: 1880px;")
    out_html_10panel.parent.mkdir(parents=True, exist_ok=True)
    out_html_10panel.write_text(html_10panel, encoding="utf-8")

    logger.info(
        f"Master validation reports generated ({out_html_4panel.name} & {out_html_10panel.name}): "
        f"{kpi['validated_total']}/{kpi['total_targets']} validated "
        f"({kpi['validated_clean']} clean, {kpi['validated_pruned']} pruned), "
        f"{kpi['failed_total']} failed, {kpi['excluded_total']} excluded."
    )


if __name__ == "__main__":
    main()
