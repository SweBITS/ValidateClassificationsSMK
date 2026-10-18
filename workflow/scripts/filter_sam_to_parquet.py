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

2. **Read-Length-Proportional Alignment Filtering**:
   Because fastp-trimmed eDNA reads vary in length (typically 50–150 bp), fixed base-pair
   cutoffs bias against shorter or longer reads. Instead, we filter every primary alignment
   (`SAM flag & 0x904 == 0`, excluding unmapped `0x4`, secondary `0x100`, and supplementary
   `0x800` records) using three read-proportional criteria:
     - **Mapping Quality**: `MAPQ >= min_mapq` (default 30, i.e. <= 0.1% probability of wrong placement).
     - **Aligned Read Fraction**: `(qlen - sclen) >= min_align_frac * qlen` (default 0.80), where
       `qlen` is total query read length and `sclen` is the sum of soft-clipped (`S`) and
       hard-clipped (`H`) bases in the CIGAR string.
     - **Maximum Mismatch Fraction**: `NM <= max_mismatch_frac * (qlen - sclen)` (default 0.01,
       corresponding to >= 99.0% sequence identity across the aligned portion of the read).

3. **Compact Polars Parquet Storage**:
   Passing alignments are written to `{taxid}.passing_alignments.parquet` via **Polars**,
   embedding the parsed `read_id`, `sample_id`, normalized `year_week` (`YYYY_WW`), `bin_id`,
   `cr_id`, `mapq`, `qlen`, `aligned_bp`, and `nm` for rapid Stage 5 aggregation.
"""

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import sample_id_to_year_week, setup_logger

CIGAR_TOKEN_RE = re.compile(r"(\d+)([MIDNSHP=X])")

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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 4: Filter SAM stream from minimap2 by MAPQ, aligned fraction, and mismatch fraction into Polars Parquet."
    )
    parser.add_argument("--taxid", type=int, required=True, help="Target species TaxID")
    parser.add_argument("--min-mapq", type=int, default=30, help="Minimum MAPQ (default: 30)")
    parser.add_argument("--min-align-frac", type=float, default=0.80, help="Minimum aligned read fraction (default: 0.80)")
    parser.add_argument("--max-mismatch-frac", type=float, default=0.01, help="Maximum mismatch (NM) fraction of aligned length (default: 0.01)")
    parser.add_argument("--out-parquet", required=True, help="Output Parquet path: {t_id}.passing_alignments.parquet")
    parser.add_argument("--out-stats-tsv", required=True, help="Output TSV path: {t_id}.mapping_stats.tsv")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger(f"filter_sam_{args.taxid}", args.log)

    total_sam_records = 0
    unmapped_or_nonprimary = 0
    primary_mapped_records = 0
    failed_mapq = 0
    failed_align_frac = 0
    failed_mismatch_frac = 0
    passing_records = 0

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
    seen_read_ids: Set[str] = set()

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

    min_mapq = args.min_mapq
    min_align_frac = args.min_align_frac
    max_mismatch_frac = args.max_mismatch_frac

    logger.info(
        f"Filtering SAM stream for TaxID {args.taxid} "
        f"(MAPQ >= {min_mapq}, align_frac >= {min_align_frac:.2f}, mismatch_frac <= {max_mismatch_frac:.4f})..."
    )

    for line in sys.stdin:
        if not line or line.startswith("@"):
            continue
        total_sam_records += 1
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            continue

        flag = int(fields[1])
        # Exclude unmapped (0x4), secondary (0x100), and supplementary (0x800) -> 0x904 = 2308
        if (flag & 0x904) != 0:
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
        if mismatch_frac - 1e-9 > max_mismatch_frac:
            failed_mismatch_frac += 1
            continue

        # Passed all Stage 4 alignment filters!
        passing_records += 1
        qname = fields[0]
        rname = fields[2]
        read_id, sample_id, year_week = parse_qname_sample(qname)
        cr_id, bin_id, source_contig = parse_rname_metadata(rname)
        mate = 1 if (flag & 0x40) else (2 if (flag & 0x80) else 0)

        col_read_id.append(read_id)
        col_sample_id.append(sample_id)
        col_year_week.append(year_week)
        col_mate.append(mate)
        col_bin_id.append(bin_id)
        col_cr_id.append(cr_id)
        col_contig.append(source_contig)
        col_mapq.append(mapq)
        col_qlen.append(qlen)
        col_aligned_bp.append(aligned_bp)
        col_align_frac.append(round(align_frac, 6))
        col_nm.append(nm)
        col_mismatch_frac.append(round(mismatch_frac, 6))
        seen_read_ids.add(read_id)

        if len(col_read_id) >= 200_000:
            flush_batch()

    flush_batch()
    writer.close()
    unique_read_pairs = len(seen_read_ids)

    out_stats = Path(args.out_stats_tsv)
    out_stats.parent.mkdir(parents=True, exist_ok=True)
    with open(out_stats, "w", encoding="utf-8", newline="") as sfh:
        writer = csv.DictWriter(
            sfh,
            fieldnames=[
                "t_id",
                "total_sam_records",
                "unmapped_or_nonprimary",
                "primary_mapped_records",
                "failed_mapq",
                "failed_align_frac",
                "failed_mismatch_frac",
                "passing_mate_alignments",
                "unique_passing_read_pairs",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({
            "t_id": args.taxid,
            "total_sam_records": total_sam_records,
            "unmapped_or_nonprimary": unmapped_or_nonprimary,
            "primary_mapped_records": primary_mapped_records,
            "failed_mapq": failed_mapq,
            "failed_align_frac": failed_align_frac,
            "failed_mismatch_frac": failed_mismatch_frac,
            "passing_mate_alignments": passing_records,
            "unique_passing_read_pairs": unique_read_pairs,
        })

    logger.info(
        f"TaxID {args.taxid}: {passing_records:,} mate alignments ({unique_read_pairs:,} unique read pairs) "
        f"passed Stage 4 filters out of {primary_mapped_records:,} primary mapped records."
    )


if __name__ == "__main__":
    main()
