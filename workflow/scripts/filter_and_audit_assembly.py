#!/usr/bin/env python3
"""
filter_and_audit_assembly.py — Stage 1: Upfront Reference Contig Filtering & Audit Gate
=======================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
Draft eukaryotic genome assemblies frequently contain thousands of short, unplaced
contigs (< 1,000 bp) that arise from unresolved repetitive elements, heterozygosity
bubbles, or trace bacterial/adapter contamination. Including < 1,000 bp contigs in our
400 bp sliding-window simulation creates two problems:
  1. Contigs shorter than 400 bp cannot generate even a single full-length 400 bp window,
     and contigs between 400 bp and 999 bp suffer severe edge-truncation bias.
  2. Short contaminant or high-copy repeat fragments distort the baseline assembly length
     and inflate noise during binning.

Therefore, Stage 1 filters every selected reference assembly upfront to retain only contigs
>= `min_contig_length` (default 1,000 bp) and records a transparent audit trail in
`assembly_filtering_summary.tsv`:
  - `raw_contig_count`, `kept_contig_count`, `discarded_contig_pct`
  - `raw_assembly_bp`, `kept_assembly_bp`, `discarded_bp`, `discarded_bp_pct`

Two automated quality gates protect downstream stages from degenerate draft assemblies:
  - **Insufficient Contig Length Gate (`EXCLUDED_INSUFFICIENT_CONTIG_LENGTH`)**:
    If `kept_assembly_bp < min_kept_assembly_bp` (default 100,000 bp, so each of the 100
    bins can hold at least ~1 kb of sequence), the species is excluded gracefully.
  - **Filtered Assembly Proportion Gate (`EXCLUDED_FILTERED_ASSEMBLY_PROPORTION`)**:
    If `discarded_bp_pct > max_discarded_bp_pct` (default 50.0% of total base pairs lost)
    or `discarded_contig_pct > max_discarded_contig_pct` (configurable, default 100.0%),
    the assembly is too fragmented to serve as a reliable genomic reference and is
    excluded gracefully.
"""

import argparse
import csv
import gzip
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import ValidationStatus, setup_logger


def filter_single_fasta_worker(
    task: Tuple[int, str, Path, Path, Path, int]
) -> Dict[str, Any]:
    """
    Worker function to stream-filter a single reference FASTA (`>= min_contig_len`),
    writing the filtered FASTA and a per-contig length table (`{t_id}.contig_lengths.tsv`).

    Returns raw and kept contig/base-pair statistics.
    """
    tid, acc, in_fasta, out_fasta, out_lengths_tsv, min_contig_len = task

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_lengths_tsv.parent.mkdir(parents=True, exist_ok=True)

    open_in = gzip.open if str(in_fasta).endswith(".gz") else open

    raw_contig_count = 0
    kept_contig_count = 0
    raw_assembly_bp = 0
    kept_assembly_bp = 0
    kept_contigs_meta: List[Tuple[str, int]] = []

    def flush_record(header_line: str, seq_chunks: List[str], out_fh) -> None:
        nonlocal raw_contig_count, kept_contig_count, raw_assembly_bp, kept_assembly_bp
        if not header_line:
            return
        seq = "".join(seq_chunks)
        seq_len = len(seq)
        raw_contig_count += 1
        raw_assembly_bp += seq_len
        if seq_len >= min_contig_len:
            kept_contig_count += 1
            kept_assembly_bp += seq_len
            # Extract clean contig ID (first whitespace-delimited token after '>')
            contig_id = header_line[1:].strip().split()[0]
            kept_contigs_meta.append((contig_id, seq_len))
            # Write clean header (`>contig_id`) and 80-char wrapped sequence so downstream
            # `seqkit sliding` produces predictable, whitespace-free FASTA headers.
            out_fh.write(f">{contig_id}\n")
            for i in range(0, seq_len, 80):
                out_fh.write(seq[i : i + 80] + "\n")

    with open_in(in_fasta, "rt", encoding="utf-8", errors="replace") as infh, open(
        out_fasta, "w", encoding="utf-8"
    ) as outfh:
        curr_header = ""
        curr_chunks: List[str] = []
        for line in infh:
            if line.startswith(">"):
                flush_record(curr_header, curr_chunks, outfh)
                curr_header = line.strip()
                curr_chunks = []
            else:
                curr_chunks.append(line.strip())
        flush_record(curr_header, curr_chunks, outfh)

    with open(out_lengths_tsv, "w", encoding="utf-8", newline="") as lfh:
        lwriter = csv.writer(lfh, delimiter="\t")
        lwriter.writerow(["contig_id", "length_bp"])
        for cid, clen in kept_contigs_meta:
            lwriter.writerow([cid, clen])

    discarded_contigs = raw_contig_count - kept_contig_count
    discarded_contig_pct = (
        (100.0 * discarded_contigs / raw_contig_count) if raw_contig_count > 0 else 100.0
    )
    discarded_bp = raw_assembly_bp - kept_assembly_bp
    discarded_bp_pct = (
        (100.0 * discarded_bp / raw_assembly_bp) if raw_assembly_bp > 0 else 100.0
    )

    return {
        "t_id": tid,
        "raw_contig_count": raw_contig_count,
        "kept_contig_count": kept_contig_count,
        "discarded_contig_pct": round(discarded_contig_pct, 4),
        "raw_assembly_bp": raw_assembly_bp,
        "kept_assembly_bp": kept_assembly_bp,
        "discarded_bp": discarded_bp,
        "discarded_bp_pct": round(discarded_bp_pct, 4),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 1: Filter contigs < min_contig_length, audit assembly retention, and emit active_species.tsv."
    )
    parser.add_argument("--stage0-ledger", required=True, help="Path to stage0_preflight_and_assemblies.tsv")
    parser.add_argument("--accession-manifest", required=True, help="Path to accession_to_fasta.tsv")
    parser.add_argument("--out-references-dir", required=True, help="Directory for filtered reference FASTAs")
    parser.add_argument("--out-audit-tsv", required=True, help="Output TSV: assembly_filtering_summary.tsv")
    parser.add_argument("--out-active-tsv", required=True, help="Output TSV: active_species.tsv (Stage 0 + 1 passing)")
    parser.add_argument("--min-contig-length", type=int, default=1000, help="Minimum contig length in bp (default: 1000)")
    parser.add_argument("--min-kept-bp", type=int, default=100000, help="Minimum retained assembly size in bp (default: 100000)")
    parser.add_argument("--max-discarded-bp-pct", type=float, default=50.0, help="Max discarded base-pair percentage (default: 50.0)")
    parser.add_argument("--max-discarded-contig-pct", type=float, default=100.0, help="Max discarded contig count percentage (default: 100.0)")
    parser.add_argument("--threads", type=int, default=4, help="Parallel worker processes")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger("filter_and_audit_assembly", args.log)

    # 1. Load Stage 0 ledger
    with open(args.stage0_ledger, "r", encoding="utf-8") as fh:
        stage0_rows = list(csv.DictReader(fh, delimiter="\t"))

    # 2. Load accession -> FASTA path mapping
    acc_to_fasta: Dict[str, Path] = {}
    with open(args.accession_manifest, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("clean_ncbi_accession") or "").strip()
            fp = (row.get("fasta_path") or "").strip()
            if acc and fp:
                acc_to_fasta[acc] = Path(fp)

    refs_dir = Path(args.out_references_dir)
    refs_dir.mkdir(parents=True, exist_ok=True)

    # 3. Build parallel tasks for all Stage 0 passing species
    tasks: List[Tuple[int, str, Path, Path, Path, int]] = []
    for row in stage0_rows:
        if row["stage0_status"] != ValidationStatus.PASSED_STAGE0:
            continue
        tid = int(row["t_id"])
        clean_acc = row["clean_ncbi_accession"]
        if clean_acc not in acc_to_fasta:
            raise FileNotFoundError(
                f"Rehydrated FASTA not found in {args.accession_manifest} for TaxID {tid} ({clean_acc})"
            )
        in_fasta = acc_to_fasta[clean_acc]
        out_fasta = refs_dir / f"{tid}.{clean_acc}.filtered.fna"
        out_lengths = refs_dir / f"{tid}.{clean_acc}.contig_lengths.tsv"
        tasks.append((tid, clean_acc, in_fasta, out_fasta, out_lengths, args.min_contig_length))

    logger.info(
        f"Filtering {len(tasks)} reference assemblies (min_contig_length >= {args.min_contig_length} bp) "
        f"using {args.threads} worker(s)..."
    )
    stats_by_tid: Dict[int, Dict[str, Any]] = {}
    if tasks:
        with ProcessPoolExecutor(max_workers=max(1, args.threads)) as pool:
            for res in pool.map(filter_single_fasta_worker, tasks):
                stats_by_tid[res["t_id"]] = res

    # 4. Evaluate Stage 1 exclusion gates and construct assembly_filtering_summary.tsv + active_species.tsv
    audit_rows: List[Dict[str, Any]] = []
    active_rows: List[Dict[str, Any]] = []

    for row in stage0_rows:
        tid = int(row["t_id"])
        s0_status = row["stage0_status"]
        s0_reason = row["stage0_reason"]
        clean_acc = row["clean_ncbi_accession"]

        if s0_status != ValidationStatus.PASSED_STAGE0:
            audit_rows.append({
                "t_id": tid,
                "t_species": row["t_species"],
                "scientific_name": row["scientific_name"],
                "species_slug": row["species_slug"],
                "chosen_accession": row["chosen_accession"],
                "clean_ncbi_accession": clean_acc,
                "chosen_assembly_taxid": row["chosen_assembly_taxid"],
                "selection_reason": row["selection_reason"],
                "raw_contig_count": 0,
                "kept_contig_count": 0,
                "discarded_contig_pct": 0.0,
                "raw_assembly_bp": 0,
                "kept_assembly_bp": 0,
                "discarded_bp": 0,
                "discarded_bp_pct": 0.0,
                "stage1_status": s0_status,
                "stage1_reason": s0_reason,
            })
            continue

        st = stats_by_tid[tid]
        kept_bp = st["kept_assembly_bp"]
        kept_contigs = st["kept_contig_count"]
        disc_bp_pct = st["discarded_bp_pct"]
        disc_contig_pct = st["discarded_contig_pct"]

        if kept_contigs == 0 or kept_bp < args.min_kept_bp:
            s1_status = ValidationStatus.EXCLUDED_INSUFFICIENT_CONTIG_LENGTH
            s1_reason = (
                f"Retained sequence after >= {args.min_contig_length} bp filter ({kept_bp:,} bp across "
                f"{kept_contigs:,} contigs) is below min_kept_assembly_bp ({args.min_kept_bp:,} bp)"
            )
            # Remove unused filtered FASTA content to save space while keeping empty placeholder if needed
        elif disc_bp_pct > args.max_discarded_bp_pct:
            s1_status = ValidationStatus.EXCLUDED_FILTERED_ASSEMBLY_PROPORTION
            s1_reason = (
                f"Discarded base-pair proportion ({disc_bp_pct:.2f}% of {st['raw_assembly_bp']:,} bp) "
                f"exceeds max_discarded_bp_pct threshold ({args.max_discarded_bp_pct:.2f}%)"
            )
        elif disc_contig_pct > args.max_discarded_contig_pct:
            s1_status = ValidationStatus.EXCLUDED_FILTERED_ASSEMBLY_PROPORTION
            s1_reason = (
                f"Discarded contig-count proportion ({disc_contig_pct:.2f}% of {st['raw_contig_count']:,} contigs) "
                f"exceeds max_discarded_contig_pct threshold ({args.max_discarded_contig_pct:.2f}%)"
            )
        else:
            s1_status = ValidationStatus.PASSED_STAGE1
            s1_reason = (
                f"Passed Stage 1: kept {kept_bp:,}/{st['raw_assembly_bp']:,} bp "
                f"({100.0 - disc_bp_pct:.2f}% bp retained; {kept_contigs:,}/{st['raw_contig_count']:,} contigs)"
            )
            active_rows.append({
                "t_id": tid,
                "t_species": row["t_species"],
                "scientific_name": row["scientific_name"],
                "species_slug": row["species_slug"],
                "chosen_accession": row["chosen_accession"],
                "clean_ncbi_accession": clean_acc,
                "filtered_fasta": str(refs_dir / f"{tid}.{clean_acc}.filtered.fna"),
                "contig_lengths_tsv": str(refs_dir / f"{tid}.{clean_acc}.contig_lengths.tsv"),
                "kept_assembly_bp": kept_bp,
                "kept_contig_count": kept_contigs,
            })

        logger.info(f"[{s1_status}] TaxID {tid} ({row['scientific_name']}): {s1_reason}")

        audit_rows.append({
            "t_id": tid,
            "t_species": row["t_species"],
            "scientific_name": row["scientific_name"],
            "species_slug": row["species_slug"],
            "chosen_accession": row["chosen_accession"],
            "clean_ncbi_accession": clean_acc,
            "chosen_assembly_taxid": row["chosen_assembly_taxid"],
            "selection_reason": row["selection_reason"],
            "raw_contig_count": st["raw_contig_count"],
            "kept_contig_count": kept_contigs,
            "discarded_contig_pct": disc_contig_pct,
            "raw_assembly_bp": st["raw_assembly_bp"],
            "kept_assembly_bp": kept_bp,
            "discarded_bp": st["discarded_bp"],
            "discarded_bp_pct": disc_bp_pct,
            "stage1_status": s1_status,
            "stage1_reason": s1_reason,
        })

    # Write assembly_filtering_summary.tsv
    Path(args.out_audit_tsv).parent.mkdir(parents=True, exist_ok=True)
    audit_fields = [
        "t_id", "t_species", "scientific_name", "species_slug",
        "chosen_accession", "clean_ncbi_accession", "chosen_assembly_taxid",
        "selection_reason", "raw_contig_count", "kept_contig_count",
        "discarded_contig_pct", "raw_assembly_bp", "kept_assembly_bp",
        "discarded_bp", "discarded_bp_pct", "stage1_status", "stage1_reason",
    ]
    with open(args.out_audit_tsv, "w", encoding="utf-8", newline="") as afh:
        writer = csv.DictWriter(afh, fieldnames=audit_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)

    # Write active_species.tsv (checkpoint output)
    Path(args.out_active_tsv).parent.mkdir(parents=True, exist_ok=True)
    active_fields = [
        "t_id", "t_species", "scientific_name", "species_slug",
        "chosen_accession", "clean_ncbi_accession", "filtered_fasta",
        "contig_lengths_tsv", "kept_assembly_bp", "kept_contig_count",
    ]
    with open(args.out_active_tsv, "w", encoding="utf-8", newline="") as act_fh:
        writer = csv.DictWriter(act_fh, fieldnames=active_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(active_rows)

    logger.info(
        f"Stage 1 audit complete: {len(active_rows)}/{len(stage0_rows)} target species passed "
        f"both Stage 0 and Stage 1 and are active for downstream validation."
    )


if __name__ == "__main__":
    main()
