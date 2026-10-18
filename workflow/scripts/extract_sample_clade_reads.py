#!/usr/bin/env python3
"""
extract_sample_clade_reads.py — Single-Pass Multi-Taxon eDNA Read Extraction per Sample
=======================================================================================

WHY THIS SCRIPT EXISTS & HOW IT AVOIDS I/O BOTTLENECKS:
-------------------------------------------------------
Our longitudinal metagenomic dataset comprises 516 weekly samples, each with ~150–200 million
paired-end reads (~10–15 GB compressed per R1/R2 pair) and a matching `.kraken2.gz` read-by-read
classification file.

A naive pipeline that extracts reads for $N$ target species by scanning each `.kraken2.gz` file
and running `seqtk subseq` $N$ separate times per sample would decompress > 250 TB of FASTQ data!
Instead, this script implements a **Single-Pass Multi-Taxon Demultiplexing Architecture**:

1. **C-Speed `.kraken2.gz` Clade Filtering (`pigz -dc | awk`)**:
   - We load the `JolTax`-expanded `descendant_taxid -> target_t_id` lookup table (from Stage 0,
     restricted to species that passed both Stage 0 and Stage 1).
   - An inline `awk` script pre-loads those target clade TaxIDs into a hash table and scans the
     decompressed `.kraken2.gz` stream at C speed (> 1 GB/s), emitting *only* rows where
     `Column 1 == "C"` and `Column 3 (classified_taxid)` belongs to an active target clade.
   - This remains 100% immune to the extra minimizer hit-group telemetry columns emitted by
     `kraken2_SweBITS`.

2. **Zero-Cost Short-Circuit for Empty Weeks**:
   - If a sample contains 0 reads across all active target species, the script writes empty
     chunk files and exits immediately in < 0.05 s without touching the 15 GB FASTQ files.

3. **Single-Pass `seqtk subseq` Union Extraction & Header Tagging**:
   - All matching `read_id`s in the sample (across all active target species combined) are written
     to a temporary `union_read_ids.txt` in `$SNIC_TMP`.
   - `seqtk subseq` streams `{sample}_R1.fastq.gz` **once** and `{sample}_R2.fastq.gz` **once**
     for the union of read IDs.
   - As `seqtk subseq` outputs the 4-line FASTQ records, Python demultiplexes each read to its
     assigned `target_t_id` chunk file while rewriting the FASTQ header strictly to:
       `@<read_id>|sample=<sample_id> [optional original comment]`
     This embeds exact weekly sample provenance directly into every read name so Stage 4 mapping
     can run once per species on the pooled FASTQ while Stage 5 reconstructs the exact 516-week
     time series.
"""

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Set

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import setup_logger


def load_active_clade_mapping(
    active_species_tsv: Path,
    clade_map_tsv: Path,
    logger,
) -> tuple[List[int], Dict[int, List[int]]]:
    """
    Load active target TaxIDs (passing Stage 0 & Stage 1) and their JolTax descendant clades.
    Returns `(active_tids, desc_to_targets)` where `desc_to_targets` maps each descendant
    TaxID to one or more active root `target_t_id`s.
    """
    active_tids: List[int] = []
    with open(active_species_tsv, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tid_str = (row.get("t_id") or "").strip()
            if tid_str:
                active_tids.append(int(tid_str))

    active_set: Set[int] = set(active_tids)
    desc_to_targets: Dict[int, List[int]] = {}
    with open(clade_map_tsv, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            dtid = int(row["descendant_taxid"])
            rtid = int(row["target_t_id"])
            if rtid in active_set:
                desc_to_targets.setdefault(dtid, []).append(rtid)

    logger.info(
        f"Loaded {len(active_tids)} active target species covering {len(desc_to_targets)} "
        f"total descendant clade TaxIDs."
    )
    return active_tids, desc_to_targets


def scan_kraken_for_sample(
    kraken_gz: Path,
    desc_to_targets: Dict[int, List[int]],
    tmp_dir: Path,
    threads: int,
    logger,
) -> Dict[str, List[int]]:
    """
    Use `pigz -dc` + `awk` to scan a sample's `.kraken2.gz` file at C speed for any read
    classified to a TaxID in `desc_to_targets`.

    Returns `read_to_targets: Dict[read_id, List[target_t_id]]`.
    """
    if not desc_to_targets:
        return {}

    # Write descendant -> comma-separated target_t_ids file for awk BEGIN block
    tax_lookup_file = tmp_dir / "active_descendant_taxids.tsv"
    with open(tax_lookup_file, "w", encoding="utf-8") as fh:
        for dtid, rtids in desc_to_targets.items():
            fh.write(f"{dtid}\t{','.join(str(x) for x in rtids)}\n")

    decomp_threads = max(1, min(threads, 4))
    awk_prog = (
        'BEGIN { FS="\\t"; OFS="\\t"; while ((getline < ARGV[1]) > 0) map[$1] = $2; ARGV[1] = "" } '
        '$1 == "C" && ($3 in map) { print $2, map[$3] }'
    )
    cmd = (
        f"pigz -dc -p {decomp_threads} {sh_quote(str(kraken_gz))} "
        f"| awk {sh_quote(awk_prog)} {sh_quote(str(tax_lookup_file))} -"
    )

    logger.info(f"Scanning {kraken_gz.name} with C-speed pigz+awk filter...")
    proc = subprocess.Popen(
        ["bash", "-o", "pipefail", "-c", cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )

    read_to_targets: Dict[str, List[int]] = {}
    assert proc.stdout is not None
    for line in proc.stdout:
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) == 2:
            read_id, rtid_csv = parts
            # Strip trailing /1 or /2 if present in Kraken read ID
            if read_id.endswith("/1") or read_id.endswith("/2"):
                read_id = read_id[:-2]
            rtids = [int(x) for x in rtid_csv.split(",")]
            if read_id not in read_to_targets:
                read_to_targets[read_id] = rtids
            else:
                for r in rtids:
                    if r not in read_to_targets[read_id]:
                        read_to_targets[read_id].append(r)

    stderr_out = proc.stderr.read() if proc.stderr else ""
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(
            f"Failed scanning {kraken_gz} (exit code {ret}): {stderr_out.strip()}"
        )

    logger.info(f"Found {len(read_to_targets):,} matching read pair(s) in {kraken_gz.name}")
    return read_to_targets


def sh_quote(s: str) -> str:
    """POSIX shell single-quote escaping."""
    return "'" + s.replace("'", "'\"'\"'") + "'"


def extract_and_demux_pair(
    fastq_gz: Path,
    union_ids_file: Path,
    read_to_targets: Dict[str, List[int]],
    active_tids: List[int],
    sample_id: str,
    pair_label: str,
    out_dir: Path,
    logger,
) -> Dict[int, int]:
    """
    Run `seqtk subseq` once on `fastq_gz` with `union_ids_file`, rewrite headers to
    `@<read_id>|sample=<sample_id> [comment]`, and demultiplex into per-species gzipped files.
    """
    counts_per_tid: Dict[int, int] = {tid: 0 for tid in active_tids}

    # Open gzipped output handles for only the species that actually have >= 1 read in this sample
    present_tids: Set[int] = set()
    for rtids in read_to_targets.values():
        present_tids.update(rtids)

    # Create empty .fq.gz files for species with 0 reads in this sample
    for tid in active_tids:
        out_fp = out_dir / f"{sample_id}_taxID-{tid}_{pair_label}.fq.gz"
        if tid not in present_tids or not read_to_targets:
            with gzip.open(out_fp, "wb", compresslevel=1):
                pass

    if not read_to_targets or not present_tids:
        return counts_per_tid

    handles: Dict[int, Any] = {}
    try:
        for tid in present_tids:
            out_fp = out_dir / f"{sample_id}_taxID-{tid}_{pair_label}.fq.gz"
            handles[tid] = gzip.open(out_fp, "wt", encoding="utf-8", compresslevel=3)

        cmd = ["seqtk", "subseq", str(fastq_gz), str(union_ids_file)]
        logger.info(
            f"Extracting {len(read_to_targets):,} {pair_label} reads from {fastq_gz.name} via seqtk subseq..."
        )
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1024 * 1024,
        )
        assert proc.stdout is not None

        while True:
            h_line = proc.stdout.readline()
            if not h_line:
                break
            seq_line = proc.stdout.readline()
            plus_line = proc.stdout.readline()
            qual_line = proc.stdout.readline()
            if not qual_line:
                raise RuntimeError(f"Truncated FASTQ record from seqtk subseq on {fastq_gz}")

            # Parse @read_id [optional comment]
            raw_header = h_line[1:].rstrip("\r\n")
            parts = raw_header.split(None, 1)
            raw_read_id = parts[0]
            comment = f" {parts[1]}" if len(parts) > 1 else ""

            lookup_id = raw_read_id
            if lookup_id.endswith("/1") or lookup_id.endswith("/2"):
                lookup_id = lookup_id[:-2]

            target_list = read_to_targets.get(lookup_id)
            if not target_list:
                continue

            # Strictly write `@<read_id>|sample=<sample_id>` as required by specification
            tagged_header = f"@{lookup_id}|sample={sample_id}{comment}\n"
            record_str = tagged_header + seq_line + plus_line + qual_line

            for tid in target_list:
                handles[tid].write(record_str)
                counts_per_tid[tid] += 1

        stderr_out = proc.stderr.read() if proc.stderr else ""
        ret = proc.wait()
        if ret != 0:
            raise RuntimeError(
                f"seqtk subseq failed on {fastq_gz} (exit code {ret}): {stderr_out.strip()}"
            )
    finally:
        for fh in handles.values():
            fh.close()

    return counts_per_tid


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract clade reads for all active species from a single sample's .kraken2.gz and R1/R2 FASTQs."
    )
    parser.add_argument("--sample-id", required=True, help="Sample identifier (e.g. Lj_2013_5_339)")
    parser.add_argument("--kraken-out", required=True, help="Path to sample .kraken2.gz file")
    parser.add_argument("--r1", required=True, help="Path to sample R1 .fastq.gz file")
    parser.add_argument("--r2", required=True, help="Path to sample R2 .fastq.gz file")
    parser.add_argument("--active-species", required=True, help="Path to active_species.tsv")
    parser.add_argument("--clade-map", required=True, help="Path to stage0_clade_map.tsv")
    parser.add_argument("--out-dir", required=True, help="Directory for per-sample species FASTQ chunks")
    parser.add_argument("--out-done", required=True, help="Output sentinel TSV recording extracted read counts per species")
    parser.add_argument("--threads", type=int, default=4, help="Threads for decompression")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger(f"extract_{args.sample_id}", args.log)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_done = Path(args.out_done)
    out_done.parent.mkdir(parents=True, exist_ok=True)

    active_tids, desc_to_targets = load_active_clade_mapping(
        Path(args.active_species), Path(args.clade_map), logger
    )

    scratch_base = os.environ.get("SNIC_TMP") or tempfile.gettempdir()
    with tempfile.TemporaryDirectory(prefix=f"edna_extract_{args.sample_id}_", dir=scratch_base) as tmp_str:
        tmp_dir = Path(tmp_str)
        read_to_targets = scan_kraken_for_sample(
            kraken_gz=Path(args.kraken_out),
            desc_to_targets=desc_to_targets,
            tmp_dir=tmp_dir,
            threads=args.threads,
            logger=logger,
        )

        union_ids_file = tmp_dir / "union_read_ids.txt"
        with open(union_ids_file, "w", encoding="utf-8") as ufh:
            for rid in sorted(read_to_targets.keys()):
                ufh.write(f"{rid}\n")

        r1_counts = extract_and_demux_pair(
            fastq_gz=Path(args.r1),
            union_ids_file=union_ids_file,
            read_to_targets=read_to_targets,
            active_tids=active_tids,
            sample_id=args.sample_id,
            pair_label="R1",
            out_dir=out_dir,
            logger=logger,
        )
        r2_counts = extract_and_demux_pair(
            fastq_gz=Path(args.r2),
            union_ids_file=union_ids_file,
            read_to_targets=read_to_targets,
            active_tids=active_tids,
            sample_id=args.sample_id,
            pair_label="R2",
            out_dir=out_dir,
            logger=logger,
        )

    with open(out_done, "w", encoding="utf-8", newline="") as dfh:
        writer = csv.writer(dfh, delimiter="\t")
        writer.writerow(["sample_id", "t_id", "extracted_r1_reads", "extracted_r2_reads"])
        for tid in active_tids:
            writer.writerow([args.sample_id, tid, r1_counts.get(tid, 0), r2_counts.get(tid, 0)])

    logger.info(
        f"Sample {args.sample_id} extraction complete (total R1 extracted across active species: "
        f"{sum(r1_counts.values()):,})."
    )


if __name__ == "__main__":
    main()
