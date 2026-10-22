#!/usr/bin/env python3
"""
compare_stage2_cr_and_benchmarks.py — Compare Baseline (4x/72c) vs. New (2x/36c + Thread Sweep) Stage 2 & Downstream Results
=============================================================================================================================
"""

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple

KRAKEN_STATS_RE = re.compile(
    r"(\d+)\s+sequences\s+\(([\d.]+)\s+Mbp\)\s+processed\s+in\s+([\d.]+)s\s+\(([\d.]+)\s+Kseq/m,\s+([\d.]+)\s+Mbp/m\)"
)


def load_bed_intervals(bed_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    if not bed_path.is_file():
        return by_contig
    with open(bed_path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) >= 3:
                by_contig.setdefault(parts[0], []).append((int(parts[1]), int(parts[2])))
    return by_contig


def interval_intersection_bp(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> int:
    """Two-pointer exact base-pair intersection between two sorted disjoint interval lists."""
    i, j = 0, 0
    inter = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if e > s:
            inter += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return inter


def compare_two_beds(old_bed: Path, new_bed: Path) -> Dict[str, object]:
    old_map = load_bed_intervals(old_bed)
    new_map = load_bed_intervals(new_bed)
    all_contigs = sorted(set(old_map.keys()) | set(new_map.keys()))

    old_bp = sum(e - s for ivs in old_map.values() for s, e in ivs)
    new_bp = sum(e - s for ivs in new_map.values() for s, e in ivs)
    old_segs = sum(len(ivs) for ivs in old_map.values())
    new_segs = sum(len(ivs) for ivs in new_map.values())

    inter_bp = 0
    for c in all_contigs:
        inter_bp += interval_intersection_bp(old_map.get(c, []), new_map.get(c, []))

    union_bp = old_bp + new_bp - inter_bp
    jaccard = (inter_bp / union_bp) if union_bp > 0 else 1.0
    return {
        "old_bp": old_bp,
        "new_bp": new_bp,
        "diff_bp": new_bp - old_bp,
        "old_segs": old_segs,
        "new_segs": new_segs,
        "inter_bp": inter_bp,
        "union_bp": union_bp,
        "only_old_bp": old_bp - inter_bp,
        "only_new_bp": new_bp - inter_bp,
        "jaccard": jaccard,
    }


def parse_kraken_log(log_path: Path) -> Dict[str, float]:
    if not log_path.is_file():
        return {}
    txt = log_path.read_text(encoding="utf-8", errors="replace")
    m = KRAKEN_STATS_RE.search(txt)
    if not m:
        return {}
    return {
        "sequences": float(m.group(1)),
        "mbp": float(m.group(2)),
        "classify_sec": float(m.group(3)),
        "kseq_per_min": float(m.group(4)),
        "mbp_per_min": float(m.group(5)),
    }


def parse_snakemake_bench(bench_path: Path) -> Dict[str, str]:
    if not bench_path.is_file():
        return {}
    with open(bench_path, "r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else {}


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare baseline vs. new Stage 2 benchmarks, thread sweep, and CR BEDs.")
    parser.add_argument("--baseline-dir", default="results/02_classifiable_regions/_baseline_4x_72t", help="Baseline Stage 2 backup directory")
    parser.add_argument("--current-cr-dir", default="results/02_classifiable_regions", help="Current Stage 2 directory")
    parser.add_argument("--current-bench", default="benchmarks/03_stream_kraken_cr/stream_kraken_cr.tsv", help="Current Snakemake benchmark TSV")
    parser.add_argument("--current-kraken-log", default="logs/03_stream_kraken_cr/kraken2_stream.log", help="Current Kraken 2 stream log")
    parser.add_argument("--thread-sweep-tsv", default="benchmarks/03_stream_kraken_cr/kraken2_thread_scaling.tsv", help="Thread scaling TSV")
    parser.add_argument("--out-md", default="results/02_classifiable_regions/stage2_comparison_report.md", help="Output markdown report")
    parser.add_argument("--out-tsv", default="results/02_classifiable_regions/stage2_cr_bed_comparison.tsv", help="Output per-species CR diff TSV")
    args = parser.parse_args()

    base_dir = Path(args.baseline_dir)
    curr_dir = Path(args.current_cr_dir)

    old_klog = parse_kraken_log(base_dir / "logs" / "kraken2_stream.log")
    new_klog = parse_kraken_log(Path(args.current_kraken_log))
    old_bench = parse_snakemake_bench(base_dir / "benchmarks" / "stream_kraken_cr.tsv")
    new_bench = parse_snakemake_bench(Path(args.current_bench))

    lines: List[str] = []
    lines.append("# Stage 2 Benchmark, Core-Scaling & Classifiable Region (CR) Comparison Report\n")

    # 1. Production Stream Benchmark Table
    lines.append("## 1. Production Stage 2 Stream Benchmark (`Baseline 4x/72t` vs. `New 2x/36t`)\n")
    lines.append("| Metric | Baseline (`100bp` step, `72` threads, `seqkit\\|awk`) | New (`200bp` step, `36` threads, `stream_tiled_fasta.py`, `intel/2025a`) | Speedup / Change |")
    lines.append("| :--- | ---: | ---: | ---: |")
    if old_klog and new_klog:
        spd_kseq = new_klog["kseq_per_min"] / max(1e-9, old_klog["kseq_per_min"])
        spd_mbp = new_klog["mbp_per_min"] / max(1e-9, old_klog["mbp_per_min"])
        spd_time = old_klog["classify_sec"] / max(1e-9, new_klog["classify_sec"])
        lines.append(f"| **Simulated Windows Classified** | {int(old_klog['sequences']):,} | {int(new_klog['sequences']):,} | {new_klog['sequences']/old_klog['sequences']:.2f}x |")
        lines.append(f"| **Total Bases Classified (Mbp)** | {old_klog['mbp']:,.2f} | {new_klog['mbp']:,.2f} | {new_klog['mbp']/old_klog['mbp']:.2f}x |")
        lines.append(f"| **Kraken 2 Classify Time (s)** | {old_klog['classify_sec']:,.1f} s ({old_klog['classify_sec']/3600:.2f} h) | {new_klog['classify_sec']:,.1f} s ({new_klog['classify_sec']/3600:.2f} h) | **{spd_time:.2f}x faster** |")
        lines.append(f"| **Throughput (`Kseq/min`)** | {old_klog['kseq_per_min']:,.1f} | {new_klog['kseq_per_min']:,.1f} | **{spd_kseq:.2f}x** |")
        lines.append(f"| **Throughput (`Mbp/min`)** | {old_klog['mbp_per_min']:,.2f} | {new_klog['mbp_per_min']:,.2f} | **{spd_mbp:.2f}x** |")
    if old_bench and new_bench:
        lines.append(f"| **Total Job Walltime (`h:m:s`)** | `{old_bench.get('h:m:s', 'N/A')}` ({float(old_bench.get('s', 0)):,.1f} s) | `{new_bench.get('h:m:s', 'N/A')}` ({float(new_bench.get('s', 0)):,.1f} s) | — |")
        lines.append(f"| **Mean CPU Load (`mean_load`)** | {old_bench.get('mean_load', 'N/A')}% | {new_bench.get('mean_load', 'N/A')}% | — |")
        lines.append(f"| **Peak RSS (`max_rss` MB)** | {float(old_bench.get('max_rss', 0)):,.1f} | {float(new_bench.get('max_rss', 0)):,.1f} | — |")
    lines.append("")

    # 2. Thread-Scaling Sweep Table
    sweep_path = Path(args.thread_sweep_tsv)
    if sweep_path.is_file():
        lines.append("## 2. Controlled Kraken 2 Core-Scaling Sweep (`6, 12, 18, 24, 36, 72` Cores)\n")
        lines.append("| Run Type | Threads | Sequences | Classify Time (s) | Throughput (`Kseq/min`) | Throughput (`Mbp/min`) | Speedup vs. 6c | Parallel Efficiency | CPU % | Voluntary Context Switches |")
        lines.append("| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        with open(sweep_path, "r", encoding="utf-8") as sfh:
            for r in csv.DictReader(sfh, delimiter="\t"):
                lines.append(
                    f"| `{r['run_type']}` | **{r['threads']}** | {int(r['sequences']):,} | "
                    f"{float(r['kraken_classify_sec']):,.2f} s | **{float(r['kseq_per_min']):,.1f}** | "
                    f"{float(r['mbp_per_min']):,.2f} | {float(r['speedup_vs_min_threads']):.2f}x | "
                    f"{float(r['parallel_efficiency_pct']):.1f}% | {float(r['cpu_percent']):.0f}% | "
                    f"{int(r['voluntary_ctxt_switches']):,} |"
                )
        lines.append("")

    # 3. Per-Species CR BED Comparison
    old_sum_path = base_dir / "cr_streaming_summary.tsv"
    cr_rows: List[Dict[str, object]] = []
    if old_sum_path.is_file():
        with open(old_sum_path, "r", encoding="utf-8") as fh:
            old_species = list(csv.DictReader(fh, delimiter="\t"))
        for sp in old_species:
            tid = sp["t_id"]
            old_bed = base_dir / "beds" / f"{tid}.cr_2x.bed"
            new_bed = curr_dir / str(tid) / f"{tid}.cr_2x.bed"
            cmp_res = compare_two_beds(old_bed, new_bed)
            cr_rows.append({
                "t_id": tid,
                "scientific_name": sp["scientific_name"],
                "kept_assembly_bp": int(sp["kept_assembly_bp"]),
                **cmp_res,
            })

    if cr_rows:
        tot_old = sum(int(r["old_bp"]) for r in cr_rows)
        tot_new = sum(int(r["new_bp"]) for r in cr_rows)
        tot_inter = sum(int(r["inter_bp"]) for r in cr_rows)
        tot_union = sum(int(r["union_bp"]) for r in cr_rows)
        tot_jacc = (tot_inter / tot_union) if tot_union > 0 else 1.0
        n_exact = sum(1 for r in cr_rows if int(r["only_old_bp"]) == 0 and int(r["only_new_bp"]) == 0)

        lines.append("## 3. Classifiable Region (`cr_2x.bed`) Concordance Across All 44 Species\n")
        lines.append(f"- **Total Old `cr_2x.bed` Base Pairs**: `{tot_old:,} bp`")
        lines.append(f"- **Total New `cr_2x.bed` Base Pairs**: `{tot_new:,} bp` (`{tot_new - tot_old:+,} bp`, `{(tot_new - tot_old) * 100.0 / max(1, tot_old):+.4f}%`)")
        lines.append(f"- **Genome-Wide Base-Pair Jaccard Similarity**: **`{tot_jacc * 100.0:.6f}%`** (`{tot_inter:,} / {tot_union:,} bp`)")
        lines.append(f"- **Species with 100.0000% Exact Base-Pair Identity**: `{n_exact} / {len(cr_rows)}`\n")

        lines.append("| TaxID | Species | Kept Assembly (bp) | Old CR (bp) | New CR (bp) | Net Diff (bp) | Old Segs | New Segs | Jaccard (%) |")
        lines.append("| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for r in cr_rows:
            lines.append(
                f"| `{r['t_id']}` | *{r['scientific_name']}* | {int(r['kept_assembly_bp']):,} | "
                f"{int(r['old_bp']):,} | {int(r['new_bp']):,} | {int(r['diff_bp']):+,} | "
                f"{int(r['old_segs']):,} | {int(r['new_segs']):,} | {float(r['jaccard']) * 100.0:.4f}% |"
            )
        lines.append("")

        out_tsv_path = Path(args.out_tsv)
        out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_tsv_path, "w", encoding="utf-8", newline="") as tfh:
            w = csv.DictWriter(tfh, fieldnames=list(cr_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(cr_rows)

    report_txt = "\n".join(lines)
    out_md_path = Path(args.out_md)
    out_md_path.parent.mkdir(parents=True, exist_ok=True)
    out_md_path.write_text(report_txt, encoding="utf-8")
    print(report_txt)


if __name__ == "__main__":
    main()
