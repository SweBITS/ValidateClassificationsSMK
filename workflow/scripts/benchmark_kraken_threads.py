#!/usr/bin/env python3
"""
benchmark_kraken_threads.py — In-Job Kraken 2 Core-Scaling Benchmark (6 to 72 Cores)
====================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
Once the main Stage 2 `stream_kraken_cr` pass finishes on a `cpu_largemem` node, the entire
2.6 TB Kraken 2 database (`hash.k2d`, `taxo.k2d`, `opts.k2d`) is warm in the Linux kernel
page cache in RAM.

By running a controlled thread-scaling sweep immediately afterward inside the same SLURM
allocation with `--memory-mapping` (which maps the already-cached RAM pages directly without
re-reading 2.6 TB from Lustre), we can measure pure Kraken 2 classification throughput across
`6, 12, 18, 24, 36, 72` OpenMP threads on Kebnekaise's 4-socket x 18-core Broadwell `largemem`
architecture in just a few minutes.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import setup_logger

KRAKEN_STATS_RE = re.compile(
    r"(\d+)\s+sequences\s+\(([\d.]+)\s+Mbp\)\s+processed\s+in\s+([\d.]+)s\s+\(([\d.]+)\s+Kseq/m,\s+([\d.]+)\s+Mbp/m\)"
)
KRAKEN_CLASS_RE = re.compile(r"(\d+)\s+sequences\s+classified\s+\(([\d.]+)%\)")


def parse_gnu_time_and_kraken_stderr(stderr_text: str) -> Dict[str, object]:
    """Extract Kraken 2 throughput metrics and GNU `/usr/bin/time -v` resource metrics."""
    res: Dict[str, object] = {
        "sequences": 0,
        "total_mbp": 0.0,
        "kraken_classify_sec": 0.0,
        "kseq_per_min": 0.0,
        "mbp_per_min": 0.0,
        "classified_sequences": 0,
        "classified_pct": 0.0,
        "user_time_sec": 0.0,
        "system_time_sec": 0.0,
        "cpu_percent": 0.0,
        "max_rss_kb": 0,
        "voluntary_ctxt_switches": 0,
        "involuntary_ctxt_switches": 0,
    }
    m_stats = KRAKEN_STATS_RE.search(stderr_text)
    if m_stats:
        res["sequences"] = int(m_stats.group(1))
        res["total_mbp"] = float(m_stats.group(2))
        res["kraken_classify_sec"] = float(m_stats.group(3))
        res["kseq_per_min"] = float(m_stats.group(4))
        res["mbp_per_min"] = float(m_stats.group(5))

    m_class = KRAKEN_CLASS_RE.search(stderr_text)
    if m_class:
        res["classified_sequences"] = int(m_class.group(1))
        res["classified_pct"] = float(m_class.group(2))

    for line in stderr_text.splitlines():
        s = line.strip()
        if s.startswith("User time (seconds):"):
            res["user_time_sec"] = float(s.split(":", 1)[1].strip())
        elif s.startswith("System time (seconds):"):
            res["system_time_sec"] = float(s.split(":", 1)[1].strip())
        elif s.startswith("Percent of CPU this job got:"):
            val = s.split(":", 1)[1].strip().rstrip("%")
            res["cpu_percent"] = float(val) if val else 0.0
        elif s.startswith("Maximum resident set size (kbytes):"):
            res["max_rss_kb"] = int(s.split(":", 1)[1].strip())
        elif s.startswith("Voluntary context switches:"):
            res["voluntary_ctxt_switches"] = int(s.split(":", 1)[1].strip())
        elif s.startswith("Involuntary context switches:"):
            res["involuntary_ctxt_switches"] = int(s.split(":", 1)[1].strip())

    return res


def build_benchmark_fasta_slice(
    tiler_script: Path,
    active_species_tsv: Path,
    out_fasta: Path,
    target_windows: int,
    window_bp: int,
    step_bp: int,
    workers: int,
    logger,
) -> int:
    """
    Generate a deterministic FASTA slice of `target_windows` (2 lines per FASTA record)
    in `$SNIC_TMP` so every thread count in the sweep classifies the exact same input.
    """
    target_lines = target_windows * 2
    logger.info(
        f"Building {target_windows:,}-window ({target_windows * window_bp / 1e6:.1f} Mbp) "
        f"benchmark slice at {out_fasta}..."
    )
    cmd = [
        sys.executable,
        str(tiler_script),
        "--active-species",
        str(active_species_tsv),
        "--window-bp",
        str(window_bp),
        "--step-bp",
        str(step_bp),
        "--workers",
        str(workers),
    ]
    lines_written = 0
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert proc.stdout is not None
    with open(out_fasta, "wb") as out_fh:
        for raw_line in proc.stdout:
            out_fh.write(raw_line)
            lines_written += 1
            if lines_written >= target_lines:
                break
        out_fh.flush()
    proc.terminate()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()
    n_seqs = lines_written // 2
    logger.info(
        f"Benchmark slice ready: {n_seqs:,} windows ({out_fasta.stat().st_size / 1e6:.1f} MB)."
    )
    return n_seqs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a controlled Kraken 2 thread-scaling benchmark (e.g. 6,12,18,24,36,72 cores) against the warm page-cached DB."
    )
    parser.add_argument("--kraken2-bin", required=True, help="Path to kraken2 executable")
    parser.add_argument("--db-path", required=True, help="Path to Kraken 2 database")
    parser.add_argument("--confidence", type=float, default=0.1, help="Kraken 2 confidence threshold")
    parser.add_argument("--mhg", type=int, default=10, help="Kraken 2 minimum hit groups")
    parser.add_argument("--tiler-script", required=True, help="Path to stream_tiled_fasta.py")
    parser.add_argument("--active-species", required=True, help="Path to active_species.tsv")
    parser.add_argument("--window-bp", type=int, default=400, help="Tile window size (bp)")
    parser.add_argument("--step-bp", type=int, default=200, help="Tile step size (bp)")
    parser.add_argument("--tiler-workers", type=int, default=8, help="Tiler workers")
    parser.add_argument("--threads-list", default="6,12,18,24,36,72", help="Comma-separated thread counts to test")
    parser.add_argument("--bench-windows", type=int, default=5_000_000, help="Number of 400bp windows for each sweep trial (default: 5,000,000 = 2.0 Gbp)")
    parser.add_argument("--main-kraken-log", default=None, help="Path to main production kraken2_stream.log to include as reference row")
    parser.add_argument("--main-threads", type=int, default=36, help="Thread count used by main production run")
    parser.add_argument("--out-tsv", required=True, help="Output TSV path for thread scaling results")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger("benchmark_kraken_threads", args.log)
    thread_counts = [int(x.strip()) for x in args.threads_list.split(",") if x.strip()]

    rows: List[Dict[str, object]] = []

    # 1. Record main production run metrics from logs/03_stream_kraken_cr/kraken2_stream.log if present
    if args.main_kraken_log and Path(args.main_kraken_log).is_file():
        main_txt = Path(args.main_kraken_log).read_text(encoding="utf-8", errors="replace")
        m_parsed = parse_gnu_time_and_kraken_stderr(main_txt)
        if int(m_parsed["sequences"]) > 0:
            rows.append({
                "run_type": "production_stream_fifo_all_species",
                "threads": args.main_threads,
                "input_mode": "FIFO_stream_2x",
                "wall_sec": round(float(m_parsed["kraken_classify_sec"]), 3),
                **m_parsed,
                "speedup_vs_min_threads": 0.0,
                "parallel_efficiency_pct": 0.0,
            })

    # 2. Build fixed benchmark slice in $SNIC_TMP (or /tmp) and run thread sweep
    tmp_root = os.environ.get("SNIC_TMP") or "/tmp"
    with tempfile.TemporaryDirectory(dir=tmp_root, prefix="kraken_bench_") as tmpdir:
        bench_fa = Path(tmpdir) / "bench_slice.fa"
        n_seqs = build_benchmark_fasta_slice(
            Path(args.tiler_script),
            Path(args.active_species),
            bench_fa,
            args.bench_windows,
            args.window_bp,
            args.step_bp,
            args.tiler_workers,
            logger,
        )
        if n_seqs == 0:
            logger.warning("No sequences generated for benchmark slice; skipping sweep.")
            return

        baseline_kseq: Optional[float] = None
        base_threads: int = thread_counts[0] if thread_counts else 1

        for t in thread_counts:
            logger.info(f"Running Kraken 2 core-scaling trial with threads={t} (OMP_NUM_THREADS={t})...")
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(t)

            cmd = [
                "/usr/bin/time",
                "-v",
                str(args.kraken2_bin),
                "--db",
                str(args.db_path),
                "--threads",
                str(t),
                "--memory-mapping",
                "--confidence",
                str(args.confidence),
                "--minimum-hit-groups",
                str(args.mhg),
                "--output",
                "-",
                str(bench_fa),
            ]
            t0 = time.monotonic()
            proc = subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                env=env,
                text=True,
                check=False,
            )
            wall_s = time.monotonic() - t0
            parsed = parse_gnu_time_and_kraken_stderr(proc.stderr or "")
            kseq_m = float(parsed["kseq_per_min"])
            if baseline_kseq is None and kseq_m > 0:
                baseline_kseq = kseq_m
                base_threads = t

            speedup = (kseq_m / baseline_kseq) if (baseline_kseq and baseline_kseq > 0) else 1.0
            ideal_scale = t / max(1, base_threads)
            eff_pct = (100.0 * speedup / ideal_scale) if ideal_scale > 0 else 100.0

            logger.info(
                f"  [threads={t:2d}] {kseq_m:,.1f} Kseq/m ({parsed['mbp_per_min']:,.2f} Mbp/m) | "
                f"classify_time={parsed['kraken_classify_sec']:.2f}s (wall={wall_s:.2f}s) | "
                f"CPU={parsed['cpu_percent']:.0f}% | speedup vs {base_threads}c = {speedup:.2f}x "
                f"(efficiency={eff_pct:.1f}%)"
            )

            rows.append({
                "run_type": "warm_cache_thread_sweep",
                "threads": t,
                "input_mode": "RAM_scratch_FASTA",
                "wall_sec": round(wall_s, 3),
                **parsed,
                "speedup_vs_min_threads": round(speedup, 4),
                "parallel_efficiency_pct": round(eff_pct, 2),
            })

    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "run_type",
        "threads",
        "input_mode",
        "sequences",
        "total_mbp",
        "kraken_classify_sec",
        "wall_sec",
        "kseq_per_min",
        "mbp_per_min",
        "speedup_vs_min_threads",
        "parallel_efficiency_pct",
        "cpu_percent",
        "user_time_sec",
        "system_time_sec",
        "voluntary_ctxt_switches",
        "involuntary_ctxt_switches",
        "max_rss_kb",
        "classified_sequences",
        "classified_pct",
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        fh.flush()
        os.fsync(fh.fileno())

    logger.info(f"Saved Kraken 2 thread-scaling benchmark table to {out_path}")


if __name__ == "__main__":
    main()
