#!/usr/bin/env python3
"""
consume_kraken_stream_to_cr.py — Stage 2: Single-Pass Kraken 2 Stream Consumer & Flanked CR Builder
===================================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
To discover the exact **Classifiable Regions (CRs)** of each target species' genome—i.e., the
genomic intervals capable of yielding a Kraken 2 classification to that species clade given our
2.6 TB database (`hash.k2d`) and parameters (`--confidence 0.1 --minimum-hit-groups 10`)—we tile
all filtered reference genomes into 400 bp single-end windows stepped every 100 bp (4x coverage)
and stream them via a zero-disk POSIX FIFO into Kraken 2.

This consumer sits at the end of the pipeline:
  `seqkit sliding -> FIFO -> kraken2 (no --report) -> awk ($1=="C") -> consume_kraken_stream_to_cr.py`

Three key algorithmic optimizations happen here in a single streaming pass:

1. **The "2x Inside 4x" Subset Identity (`idx % 2 == 0`)**:
   In a 400 bp window / 100 bp step (4x coverage) tiling along each contig (`idx = 0, 1, 2, 3, ...`),
   the even-numbered windows (`idx % 2 == 0`, starting at 0, 200, 400, ...) ARE the exact 200 bp
   step (2x coverage) dataset! Thus Kraken 2 classifies the 4x stream once against the 2.6 TB
   database, and this consumer simultaneously constructs both the **2x CR BED** and the **4x CR BED**
   (plus 1x coverage statistics via `idx % 4 == 0`).

2. **Clade-Aware Classification Check via `JolTax`**:
   A simulated window from a species reference genome is considered "correctly classified" if
   Kraken 2 assigns it to the target species TaxID (`t_id`) OR any descendant subspecies/strain
   TaxID within its `JolTax` clade.

3. **Step-Size Flanking Buffers ($\\pm S$) with Contig-End Clamping & Interval Union**:
   Why pad each classified window $[s, e)$ by the step size $\\pm S$ (`+/- 200 bp` for 2x;
   `+/- 100 bp` for 4x)?
   - Because windows are stepped by $S$ bp, a diagnostic minimizer cluster near the boundary of a
     classified window can also be overlapped by real eDNA reads extending up to $S$ bp outside
     that window before the next adjacent window is sampled. Padding by $\\pm S$ prevents real
     eDNA reads at CR boundaries from being soft-clipped when aligned exclusively to excised CRs.
   - Every padded interval is strictly clamped to `[0, contig_length]` so coordinates never become
     negative or overrun the contig end.
   - Finally, all overlapping or adjacent padded intervals (`next_start <= curr_end`) on each contig
     are merged into disjoint intervals `[start, end)` before sequence extraction, guaranteeing
     **zero base-pair duplication or double-counting**.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import setup_logger


def compact_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Sort and merge overlapping or touching 0-based half-open intervals `[start, end)`.
    Uses strict `next_start <= curr_end` so adjacent intervals `[0, 600)` and `[600, 1000)`
    coalesce into a single continuous interval `[0, 1000)`.
    """
    if len(intervals) <= 1:
        return intervals
    intervals.sort(key=lambda x: (x[0], x[1]))
    merged: List[Tuple[int, int]] = [intervals[0]]
    for s, e in intervals[1:]:
        cs, ce = merged[-1]
        if s <= ce:
            if e > ce:
                merged[-1] = (cs, e)
        else:
            merged.append((s, e))
    return merged


def load_active_metadata(
    active_species_tsv: Path,
    clade_map_tsv: Path,
    logger,
) -> Tuple[Dict[int, Dict[str, str]], Dict[int, Set[int]], Dict[int, Dict[str, int]]]:
    """
    Load active species rows, per-species JolTax clade sets, and per-species contig length maps.
    """
    active_meta: Dict[int, Dict[str, str]] = {}
    contig_lengths: Dict[int, Dict[str, int]] = {}

    with open(active_species_tsv, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tid = int(row["t_id"])
            active_meta[tid] = row
            lengths_path = Path(row["contig_lengths_tsv"])
            clen_map: Dict[str, int] = {}
            if lengths_path.is_file():
                with open(lengths_path, "r", encoding="utf-8") as lfh:
                    for lrow in csv.DictReader(lfh, delimiter="\t"):
                        clen_map[lrow["contig_id"]] = int(lrow["length_bp"])
            contig_lengths[tid] = clen_map

    clade_sets: Dict[int, Set[int]] = {tid: {tid} for tid in active_meta}
    with open(clade_map_tsv, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            dtid = int(row["descendant_taxid"])
            rtid = int(row["target_t_id"])
            if rtid in clade_sets:
                clade_sets[rtid].add(dtid)

    logger.info(
        f"Initialized CR stream consumer for {len(active_meta)} active species "
        f"({sum(len(v) for v in contig_lengths.values()):,} total filtered contigs)."
    )
    return active_meta, clade_sets, contig_lengths


def write_and_fsync_bed(
    bed_path: Path,
    contig_intervals: Dict[str, List[Tuple[int, int]]],
) -> Tuple[int, int]:
    """
    Write merged 0-based BED intervals (`contig \t start \t end \t cr_length`) to disk,
    flushing and calling `os.fsync()` to guarantee Lustre/NFS persistence before SLURM job exit.
    Returns `(total_cr_bp, total_cr_segments)`.
    """
    bed_path.parent.mkdir(parents=True, exist_ok=True)
    total_bp = 0
    total_segs = 0

    with open(bed_path, "w", encoding="utf-8") as fh:
        for contig in sorted(contig_intervals.keys()):
            merged = compact_intervals(contig_intervals[contig])
            for s, e in merged:
                seg_len = e - s
                if seg_len > 0:
                    fh.write(f"{contig}\t{s}\t{e}\t{seg_len}\n")
                    total_bp += seg_len
                    total_segs += 1
        fh.flush()
        os.fsync(fh.fileno())

    return total_bp, total_segs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Consume awk-filtered Kraken 2 stdout stream to build 2x and 4x Classifiable Region (CR) BED files."
    )
    parser.add_argument("--active-species", required=True, help="Path to active_species.tsv")
    parser.add_argument("--clade-map", required=True, help="Path to stage0_clade_map.tsv")
    parser.add_argument("--out-cr-dir", required=True, help="Base output directory for per-species CR BEDs (results/02_classifiable_regions)")
    parser.add_argument("--out-summary-tsv", required=True, help="Output TSV summarizing 2x and 4x CR statistics per species")
    parser.add_argument("--out-sentinel", required=True, help="Output completion sentinel file verified after fsync")
    parser.add_argument("--pad-2x", type=int, default=200, help="Flanking buffer in bp for 2x tier (default: 200)")
    parser.add_argument("--pad-4x", type=int, default=100, help="Flanking buffer in bp for 4x tier (default: 100)")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger("consume_kraken_stream_to_cr", args.log)

    active_meta, clade_sets, contig_lengths = load_active_metadata(
        Path(args.active_species), Path(args.clade_map), logger
    )

    # Per-species, per-contig interval accumulators
    intervals_2x: Dict[int, Dict[str, List[Tuple[int, int]]]] = {
        tid: defaultdict(list) for tid in active_meta
    }
    intervals_4x: Dict[int, Dict[str, List[Tuple[int, int]]]] = {
        tid: defaultdict(list) for tid in active_meta
    }

    # Diagnostic counters
    windows_classified_total: Dict[int, int] = {tid: 0 for tid in active_meta}
    windows_clade_4x: Dict[int, int] = {tid: 0 for tid in active_meta}
    windows_clade_2x: Dict[int, int] = {tid: 0 for tid in active_meta}

    pad_2x = args.pad_2x
    pad_4x = args.pad_4x
    lines_processed = 0

    logger.info("Draining awk-filtered Kraken 2 classification stream from stdin...")
    for line in sys.stdin:
        line = line.rstrip("\r\n")
        if not line:
            continue
        lines_processed += 1
        if lines_processed % 5_000_000 == 0:
            logger.info(f"Processed {lines_processed:,} classified simulated windows from stream...")

        parts = line.split("\t")
        if len(parts) < 2:
            continue
        header_str, k2_taxid_str = parts[0], parts[1]

        # Header format injected during FIFO streaming:
        # `{t_id}|{contig}|{idx}|{start0}|{end0}`
        try:
            tid_part, rest = header_str.split("|", 1)
            contig_id, idx_str, s0_str, e0_str = rest.rsplit("|", 3)
            tid = int(tid_part)
            k2_taxid = int(k2_taxid_str)
            idx = int(idx_str)
            start0 = int(s0_str)
            end0 = int(e0_str)
        except ValueError:
            continue

        if tid not in active_meta:
            continue

        windows_classified_total[tid] += 1

        # Verify that Kraken 2 classified this window to the target species clade
        if k2_taxid not in clade_sets[tid]:
            continue

        clen = contig_lengths[tid].get(contig_id)
        if clen is None:
            continue

        # 1. Record 4x interval padded by +/- pad_4x (100 bp) and clamped to [0, clen]
        windows_clade_4x[tid] += 1
        s4 = start0 - pad_4x
        if s4 < 0:
            s4 = 0
        e4 = end0 + pad_4x
        if e4 > clen:
            e4 = clen
        list_4x = intervals_4x[tid][contig_id]
        list_4x.append((s4, e4))
        if len(list_4x) >= 32768:
            intervals_4x[tid][contig_id] = compact_intervals(list_4x)

        # 2. If even window index (idx % 2 == 0), record 2x interval padded by +/- pad_2x (200 bp)
        if (idx & 1) == 0:
            windows_clade_2x[tid] += 1
            s2 = start0 - pad_2x
            if s2 < 0:
                s2 = 0
            e2 = end0 + pad_2x
            if e2 > clen:
                e2 = clen
            list_2x = intervals_2x[tid][contig_id]
            list_2x.append((s2, e2))
            if len(list_2x) >= 32768:
                intervals_2x[tid][contig_id] = compact_intervals(list_2x)

    logger.info(
        f"Finished draining stdin ({lines_processed:,} total classified windows). "
        f"Merging intervals and writing 2x and 4x CR BED files..."
    )

    out_cr_dir = Path(args.out_cr_dir)
    summary_rows: List[Dict[str, object]] = []

    for tid, meta in active_meta.items():
        sp_dir = out_cr_dir / str(tid)
        bed_2x_path = sp_dir / f"{tid}.cr_2x.bed"
        bed_4x_path = sp_dir / f"{tid}.cr_4x.bed"

        cr_bp_2x, cr_segs_2x = write_and_fsync_bed(bed_2x_path, intervals_2x[tid])
        cr_bp_4x, cr_segs_4x = write_and_fsync_bed(bed_4x_path, intervals_4x[tid])

        kept_bp = int(meta.get("kept_assembly_bp") or 0)
        pct_2x = (100.0 * cr_bp_2x / kept_bp) if kept_bp > 0 else 0.0
        pct_4x = (100.0 * cr_bp_4x / kept_bp) if kept_bp > 0 else 0.0

        logger.info(
            f"TaxID {tid} ({meta['scientific_name']}): "
            f"2x CR = {cr_bp_2x:,} bp ({pct_2x:.2f}% of filtered assembly across {cr_segs_2x:,} segments); "
            f"4x CR = {cr_bp_4x:,} bp ({pct_4x:.2f}% across {cr_segs_4x:,} segments)"
        )

        summary_rows.append({
            "t_id": tid,
            "t_species": meta["t_species"],
            "scientific_name": meta["scientific_name"],
            "species_slug": meta["species_slug"],
            "kept_assembly_bp": kept_bp,
            "windows_classified_any_taxon_4x": windows_classified_total[tid],
            "windows_clade_classified_4x": windows_clade_4x[tid],
            "windows_clade_classified_2x": windows_clade_2x[tid],
            "cr_bp_2x": cr_bp_2x,
            "cr_pct_of_filtered_assembly_2x": round(pct_2x, 4),
            "cr_segments_2x": cr_segs_2x,
            "cr_bp_4x": cr_bp_4x,
            "cr_pct_of_filtered_assembly_4x": round(pct_4x, 4),
            "cr_segments_4x": cr_segs_4x,
        })

    # Write summary TSV with explicit flush + fsync
    summary_path = Path(args.out_summary_tsv)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "t_id", "t_species", "scientific_name", "species_slug", "kept_assembly_bp",
        "windows_classified_any_taxon_4x", "windows_clade_classified_4x",
        "windows_clade_classified_2x", "cr_bp_2x", "cr_pct_of_filtered_assembly_2x",
        "cr_segments_2x", "cr_bp_4x", "cr_pct_of_filtered_assembly_4x", "cr_segments_4x",
    ]
    with open(summary_path, "w", encoding="utf-8", newline="") as sfh:
        writer = csv.DictWriter(sfh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
        sfh.flush()
        os.fsync(sfh.fileno())

    # Write explicit completion token after all files are fsynced
    sentinel_path = Path(args.out_sentinel)
    sentinel_path.parent.mkdir(parents=True, exist_ok=True)
    with open(sentinel_path, "w", encoding="utf-8") as sent_fh:
        sent_fh.write("STREAM_DRAINED_AND_FSYNCED_OK\n")
        sent_fh.flush()
        os.fsync(sent_fh.fileno())

    logger.info("All CR BEDs and summary tables flushed and fsynced successfully.")


if __name__ == "__main__":
    main()
