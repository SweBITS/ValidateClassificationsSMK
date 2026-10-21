#!/usr/bin/env python3
"""
filter_sam_to_parquet.py — Stage 4: Read-Length-Proportional SAM Stream Filter to Polars Parquet
================================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
After `minimap2 -ax sr` aligns a species' pooled paired-end eDNA reads (`R1`, `R2`) exclusively
to its binned Classifiable Region sequences (`{taxid}.cr_contigs.fa`), we must separate genuine
homologous alignments from partial/clipped cross-mappings to conserved motifs:

1. **No PCR Deduplication (`markdup` Omitted Intentionally)**:
   In high-abundance airborne eDNA weeks (e.g. fungal spore peaks or pollen bursts), thousands
   of genuine biological molecules shear at identical intervals across compact CRs. Removing
   coordinate "duplicates" artificially saturates high-read bins and flattens the linear
   relationship between CR length and read count.

2. **Strict Proper-Pair Same-CR Alignment Filtering**:
   Because fastp-trimmed eDNA reads vary in length (typically 50–150 bp), fixed base-pair
   cutoffs bias against shorter or longer reads. We require every counted fragment to be a
   **proper pair (`0x2`) mapping within the exact same Classifiable Region (`cr_id`, `RNEXT == "="`)**
   where **BOTH R1 and R2** are primary alignments (`SAM flag & 0x90C == 0`, excluding unmapped
   `0x4`, mate-unmapped `0x8`, secondary `0x100`, and supplementary `0x800` records) and **BOTH**
   mates individually satisfy all three read-proportional criteria:
     - **Mapping Quality**: `MAPQ >= min_mapq` (default 30, i.e. <= 0.1% probability of wrong placement).
     - **Aligned Read Fraction**: `(qlen - sclen) >= min_align_frac * qlen` (default 0.80), where
       `qlen` is total query read length and `sclen` is the sum of soft-clipped (`S`) and
       hard-clipped (`H`) bases in the CIGAR string.
     - **Maximum Mismatch Fraction**: `NM <= max_mismatch_frac * (qlen - sclen)` (default 0.01,
       corresponding to >= 99.0% sequence identity across the aligned portion of the read).

3. **O(1)-Memory Adjacent Mate Pairing & Compact Polars Parquet Storage**:
   Because `minimap2 -ax sr` emits R1 and R2 of each read pair consecutively, we pair passing R1
   and R2 alignments on the fly in $O(1)$ memory (eliminating multi-gigabyte hash sets on 50M-read
   taxa) and write only complete 2-mate proper pairs to `{taxid}.passing_alignments.parquet`.
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import sample_id_to_year_week, setup_logger

CIGAR_TOKEN_RE = re.compile(r"(\d+)([MIDNSHP=X])")
MINIMAP2_PROCESSED_RE = re.compile(r"(?:mapped|processed)\s+(\d+)\s+sequences")

# Pre-defined Polars schema so even a 0-alignment species produces a valid, typed Parquet file
ALIGNMENT_SCHEMA = {
    "read_id": pl.Utf8,
    "sample_id": pl.Utf8,
    "year_week": pl.Utf8,
    "mate": pl.Int8,
    "bin_id": pl.Utf8,
    "cr_id": pl.Utf8,
    "source_contig": pl.Utf8,
    "mapq": pl.Int32,
    "qlen": pl.Int32,
    "aligned_bp": pl.Int32,
    "align_frac": pl.Float64,
    "nm": pl.Int32,
    "mismatch_frac": pl.Float64,
}

PYARROW_SCHEMA = pa.schema([
    ("read_id", pa.string()),
    ("sample_id", pa.string()),
    ("year_week", pa.string()),
    ("mate", pa.int8()),
    ("bin_id", pa.string()),
    ("cr_id", pa.string()),
    ("source_contig", pa.string()),
    ("mapq", pa.int32()),
    ("qlen", pa.int32()),
    ("aligned_bp", pa.int32()),
    ("align_frac", pa.float64()),
    ("nm", pa.int32()),
    ("mismatch_frac", pa.float64()),
])


def parse_cigar_qlen_and_sclen(cigar: str, seq_str: str) -> Tuple[int, int]:
    """
    Compute total query read length (`qlen`) and clipped bases (`sclen`) from a SAM CIGAR string.
    Query-consuming CIGAR operations are `M`, `I`, `S`, `=`, `X` (plus `H` for hard-clipped length).
    """
    if cigar == "*":
        return len(seq_str), len(seq_str)
    qlen = 0
    sclen = 0
    for length_str, op in CIGAR_TOKEN_RE.findall(cigar):
        n = int(length_str)
        if op in ("M", "I", "=", "X"):
            qlen += n
        elif op in ("S", "H"):
            qlen += n
            sclen += n
    if qlen == 0 and seq_str != "*":
        qlen = len(seq_str)
    return qlen, sclen


def extract_nm_tag(optional_fields: List[str]) -> Optional[int]:
    """Extract the integer edit distance `NM:i:<int>` from SAM optional tag fields."""
    for field in optional_fields:
        if field.startswith("NM:i:") or field.startswith("nM:i:"):
            try:
                return int(field[5:])
            except ValueError:
                return None
    return None


def parse_qname_sample(qname: str) -> Tuple[str, str, str]:
    """
    Parse `<read_id>|sample=<sample_id>` from SAM `QNAME` (column 1).
    Returns `(base_read_id, sample_id, year_week)`.
    """
    if "|sample=" in qname:
        read_id, sample_part = qname.rsplit("|sample=", 1)
        sample_id = sample_part.strip()
    elif "::" in qname:
        # Fallback compatibility if a legacy header is encountered
        sample_id, read_id = qname.split("::", 1)
    else:
        read_id = qname
        sample_id = "UNKNOWN"

    yw = sample_id_to_year_week(sample_id) or "UNKNOWN"
    return read_id, sample_id, yw


def parse_rname_metadata(rname: str) -> Tuple[str, str, str]:
    """
    Parse `cr_id`, `bin_id`, and `source_contig` from `{taxid}.cr_contigs.fa` header:
      `cr_000001|bin=bin_001|contig=NC_000001.1|coords=0-1200`
    """
    cr_id = rname
    bin_id = "UNKNOWN"
    contig = "UNKNOWN"
    for part in rname.split("|"):
        if part.startswith("cr_"):
            cr_id = part
        elif part.startswith("bin="):
            bin_id = part[4:]
        elif part.startswith("contig="):
            contig = part[7:]
    return cr_id, bin_id, contig


def parse_minimap2_total_sequences(log_path: Optional[str]) -> int:
    """Sum `[M::worker_pipeline::*] processed <N> sequences` from minimap2 stderr log if available."""
    if not log_path:
        return 0
    p = Path(log_path)
    if not p.is_file():
        return 0
    total = 0
    try:
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                m = MINIMAP2_PROCESSED_RE.search(line)
                if m:
                    total += int(m.group(1))
    except Exception:
        return 0
    return total


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 4: Filter SAM stream from minimap2 for strict same-CR proper pairs into Polars Parquet."
    )
    parser.add_argument("--taxid", type=int, required=True, help="Target species TaxID")
    parser.add_argument("--min-mapq", type=int, default=30, help="Minimum MAPQ (default: 30)")
    parser.add_argument("--min-align-frac", type=float, default=0.80, help="Minimum aligned read fraction (default: 0.80)")
    parser.add_argument("--max-mismatch-frac", type=float, default=0.01, help="Maximum mismatch (NM) fraction of aligned length (default: 0.01)")
    parser.add_argument("--minimap2-log", default=None, help="Optional path to minimap2 stderr log for total sequence count")
    parser.add_argument("--out-parquet", required=True, help="Output Parquet path: {t_id}.passing_alignments.parquet")
    parser.add_argument("--out-stats-tsv", required=True, help="Output TSV path: {t_id}.mapping_stats.tsv")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger(f"filter_sam_{args.taxid}", args.log)

    # Preserve hardlinks of pooled species FASTQs so even if an in-flight Snakemake coordinator
    # had `temp(...)` in memory, the merged species FASTQs survive for future mapping parameter runs.
    fastq_dir = Path("results/01_inputs/fastq")
    for pair in ("1", "2"):
        src_fq = fastq_dir / f"taxID-{args.taxid}_R{pair}.fq.gz"
        keep_fq = fastq_dir / f".retained_taxID-{args.taxid}_R{pair}.fq.gz"
        if src_fq.is_file() and not keep_fq.exists():
            try:
                keep_fq.hardlink_to(src_fq)
            except Exception:
                pass

    total_sam_records = 0
    unmapped_or_nonprimary = 0
    primary_mapped_records = 0
    failed_mapq = 0
    failed_align_frac = 0
    failed_mismatch_frac = 0
    failed_mate_pair = 0
    passing_records = 0
    unique_read_pairs = 0

    # Columnar buffers for streaming PyArrow Parquet writer
    col_read_id: List[str] = []
    col_sample_id: List[str] = []
    col_year_week: List[str] = []
    col_mate: List[int] = []
    col_bin_id: List[str] = []
    col_cr_id: List[str] = []
    col_contig: List[str] = []
    col_mapq: List[int] = []
    col_qlen: List[int] = []
    col_aligned_bp: List[int] = []
    col_align_frac: List[float] = []
    col_nm: List[int] = []
    col_mismatch_frac: List[float] = []

    out_parquet = Path(args.out_parquet)
    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    writer = pq.ParquetWriter(out_parquet, PYARROW_SCHEMA, compression="zstd", compression_level=3)

    def flush_batch() -> None:
        if not col_read_id:
            return
        batch_table = pa.Table.from_arrays(
            [
                pa.array(col_read_id, type=pa.string()),
                pa.array(col_sample_id, type=pa.string()),
                pa.array(col_year_week, type=pa.string()),
                pa.array(col_mate, type=pa.int8()),
                pa.array(col_bin_id, type=pa.string()),
                pa.array(col_cr_id, type=pa.string()),
                pa.array(col_contig, type=pa.string()),
                pa.array(col_mapq, type=pa.int32()),
                pa.array(col_qlen, type=pa.int32()),
                pa.array(col_aligned_bp, type=pa.int32()),
                pa.array(col_align_frac, type=pa.float64()),
                pa.array(col_nm, type=pa.int32()),
                pa.array(col_mismatch_frac, type=pa.float64()),
            ],
            schema=PYARROW_SCHEMA,
        )
        writer.write_table(batch_table)
        col_read_id.clear()
        col_sample_id.clear()
        col_year_week.clear()
        col_mate.clear()
        col_bin_id.clear()
        col_cr_id.clear()
        col_contig.clear()
        col_mapq.clear()
        col_qlen.clear()
        col_aligned_bp.clear()
        col_align_frac.clear()
        col_nm.clear()
        col_mismatch_frac.clear()

    def append_record(rec: Tuple[str, str, str, int, str, str, str, int, int, int, float, int, float]) -> None:
        col_read_id.append(rec[0])
        col_sample_id.append(rec[1])
        col_year_week.append(rec[2])
        col_mate.append(rec[3])
        col_bin_id.append(rec[4])
        col_cr_id.append(rec[5])
        col_contig.append(rec[6])
        col_mapq.append(rec[7])
        col_qlen.append(rec[8])
        col_aligned_bp.append(rec[9])
        col_align_frac.append(rec[10])
        col_nm.append(rec[11])
        col_mismatch_frac.append(rec[12])

    min_mapq = args.min_mapq
    min_align_frac = args.min_align_frac
    max_mismatch_frac = args.max_mismatch_frac

    logger.info(
        f"Filtering SAM stream for TaxID {args.taxid} (strict same-CR proper pairs: flag 0x2, "
        f"RNEXT='=', both mates MAPQ >= {min_mapq}, align_frac >= {min_align_frac:.2f}, "
        f"mismatch_frac <= {max_mismatch_frac:.4f})..."
    )

    # Buffer for the first passing mate of a read pair
    pending_mate: Optional[Tuple[str, str, str, int, str, str, str, int, int, int, float, int, float]] = None
    # Buffer for any mate passing MAPQ + align_frac so we can compute exact pair-minimum identity
    # histogram across <90%, 90%, 91%, ..., 100% for pairs rejected at the identity gate
    pending_align_mate: Optional[Tuple[str, str, int, float, bool]] = None
    lt99_pair_hist: Dict[int, int] = {b: 0 for b in range(89, 101)}
    lt99_mate_hist: Dict[int, int] = {b: 0 for b in range(89, 101)}

    for line in sys.stdin:
        if not line or line.startswith("@"):
            continue
        total_sam_records += 1
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            continue

        flag = int(fields[1])
        rname = fields[2]
        rnext = fields[6]

        # Require proper pair (0x2), exclude unmapped (0x4), mate unmapped (0x8),
        # secondary (0x100), and supplementary (0x800) -> 0x90C = 2316,
        # and require both mates to map to the exact same CR contig (rnext == "=" or rnext == rname)
        if (flag & 0x90C) != 0 or (flag & 0x2) == 0 or (rnext != "=" and rnext != rname):
            unmapped_or_nonprimary += 1
            continue

        primary_mapped_records += 1
        mapq = int(fields[4])
        if mapq < min_mapq:
            failed_mapq += 1
            continue

        cigar = fields[5]
        seq_str = fields[9]
        qlen, sclen = parse_cigar_qlen_and_sclen(cigar, seq_str)
        aligned_bp = qlen - sclen
        if qlen <= 0 or aligned_bp <= 0:
            failed_align_frac += 1
            continue

        # Use float epsilon 1e-9 to avoid binary floating-point boundary rejections (e.g. 120/150 == 0.80)
        align_frac = aligned_bp / float(qlen)
        if align_frac + 1e-9 < min_align_frac:
            failed_align_frac += 1
            continue

        nm = extract_nm_tag(fields[11:])
        if nm is None:
            failed_mismatch_frac += 1
            continue

        mismatch_frac = nm / float(aligned_bp)
        id_pct = 100.0 * (1.0 - mismatch_frac)
        mate_passed_id = (mismatch_frac - 1e-9 <= max_mismatch_frac)

        qname = fields[0]
        read_id, sample_id, year_week = parse_qname_sample(qname)
        cr_id, bin_id, source_contig = parse_rname_metadata(rname)
        mate = 1 if (flag & 0x40) else (2 if (flag & 0x80) else 0)

        if not mate_passed_id:
            mate_bin = 89 if int(math.floor(id_pct + 1e-9)) < 90 else min(100, int(math.floor(id_pct + 1e-9)))
            lt99_mate_hist[mate_bin] += 1

        # Pair-level rejected identity histogram tracking for complete 2-mate proper pairs
        # where BOTH mates passed MAPQ + aligned length on the same CR
        if pending_align_mate is None:
            pending_align_mate = (read_id, cr_id, mate, id_pct, mate_passed_id)
        elif (
            pending_align_mate[0] == read_id
            and pending_align_mate[1] == cr_id
            and pending_align_mate[2] != mate
            and pending_align_mate[2] in (1, 2)
            and mate in (1, 2)
        ):
            if not (pending_align_mate[4] and mate_passed_id):
                pair_min_id = min(pending_align_mate[3], id_pct)
                pair_bin = 89 if int(math.floor(pair_min_id + 1e-9)) < 90 else min(100, int(math.floor(pair_min_id + 1e-9)))
                lt99_pair_hist[pair_bin] += 1
            pending_align_mate = None
        else:
            # Previous pending_align_mate's partner failed MAPQ or align_frac; do not count single-mate orphans
            pending_align_mate = (read_id, cr_id, mate, id_pct, mate_passed_id)

        if not mate_passed_id:
            failed_mismatch_frac += 1
            continue

        # Individual mate passed all filters; pair with its mate on the same CR
        curr_rec = (
            read_id,
            sample_id,
            year_week,
            mate,
            bin_id,
            cr_id,
            source_contig,
            mapq,
            qlen,
            aligned_bp,
            round(align_frac, 6),
            nm,
            round(mismatch_frac, 6),
        )

        if pending_mate is None:
            pending_mate = curr_rec
        elif (
            pending_mate[0] == read_id
            and pending_mate[5] == cr_id
            and pending_mate[3] != mate
            and pending_mate[3] in (1, 2)
            and mate in (1, 2)
        ):
            # Both R1 and R2 of this proper pair passed all filters on the exact same CR!
            append_record(pending_mate)
            append_record(curr_rec)
            pending_mate = None
            passing_records += 2
            unique_read_pairs += 1
            if len(col_read_id) >= 200_000:
                flush_batch()
        else:
            # Previous pending mate's partner did not pass filters; replace with current mate
            failed_mate_pair += 1
            pending_mate = curr_rec

    pending_align_mate = None

    if pending_mate is not None:
        failed_mate_pair += 1
        pending_mate = None

    flush_batch()
    writer.close()

    # If samtools pre-filtered the stream, recover total input sequences from minimap2 log if available
    mm2_total = parse_minimap2_total_sequences(args.minimap2_log)
    if mm2_total > total_sam_records:
        unmapped_or_nonprimary += (mm2_total - total_sam_records)
        total_sam_records = mm2_total

    out_stats = Path(args.out_stats_tsv)
    out_stats.parent.mkdir(parents=True, exist_ok=True)
    hist_fieldnames = ["id_lt90_pairs"] + [f"id_{b}_pairs" for b in range(90, 101)]
    row_dict = {
        "t_id": args.taxid,
        "min_mapq": min_mapq,
        "min_align_frac": min_align_frac,
        "max_mismatch_frac": max_mismatch_frac,
        "total_sam_records": total_sam_records,
        "unmapped_or_nonprimary": unmapped_or_nonprimary,
        "primary_mapped_records": primary_mapped_records,
        "failed_mapq": failed_mapq,
        "failed_align_frac": failed_align_frac,
        "failed_mismatch_frac": failed_mismatch_frac + failed_mate_pair,
        "passing_mate_alignments": passing_records,
        "unique_passing_read_pairs": unique_read_pairs,
        "id_lt90_pairs": lt99_pair_hist[89],
    }
    for b in range(90, 101):
        row_dict[f"id_{b}_pairs"] = lt99_pair_hist[b]

    with open(out_stats, "w", encoding="utf-8", newline="") as sfh:
        stats_writer = csv.DictWriter(
            sfh,
            fieldnames=[
                "t_id",
                "min_mapq",
                "min_align_frac",
                "max_mismatch_frac",
                "total_sam_records",
                "unmapped_or_nonprimary",
                "primary_mapped_records",
                "failed_mapq",
                "failed_align_frac",
                "failed_mismatch_frac",
                "passing_mate_alignments",
                "unique_passing_read_pairs",
            ] + hist_fieldnames,
            delimiter="\t",
        )
        stats_writer.writeheader()
        stats_writer.writerow(row_dict)

    # Also write the dedicated 1%-bin rejected identity histogram TSV alongside mapping_stats.tsv
    out_hist_tsv = out_stats.parent / f"{args.taxid}.identity_histogram_lt99.tsv"
    tot_lt99_pairs = max(1, sum(lt99_pair_hist.values()))
    with open(out_hist_tsv, "w", encoding="utf-8", newline="") as hfh:
        hw = csv.DictWriter(
            hfh,
            fieldnames=["id_bin", "label", "pair_count", "mate_count", "pair_share"],
            delimiter="\t",
        )
        hw.writeheader()
        for b in range(89, 101):
            lbl = "<90" if b == 89 else str(b)
            hw.writerow({
                "id_bin": b,
                "label": lbl,
                "pair_count": lt99_pair_hist[b],
                "mate_count": lt99_mate_hist[b],
                "pair_share": f"{lt99_pair_hist[b] / float(tot_lt99_pairs):.6f}",
            })

    logger.info(
        f"TaxID {args.taxid}: {unique_read_pairs:,} strict proper read pairs ({passing_records:,} mate alignments) "
        f"mapped within the same CR and passed all filters (dropped {failed_mate_pair:,} single-mate survivors)."
    )


if __name__ == "__main__":
    main()
