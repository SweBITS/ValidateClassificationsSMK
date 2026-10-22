#!/usr/bin/env python3
"""
stream_tiled_fasta.py — High-Throughput Parallel In Silico Genome Tiler for Stage 2
====================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
In Kraken 2 (`classify.cc`), worker threads acquire `#pragma omp critical(seqread)` and
call `ifs.read(block_buffer_, 3 * 1024 * 1024)` (3 MB) to fetch the next batch of sequences.
Previously, feeding a POSIX FIFO with a single-threaded `seqkit sliding | awk match(...)`
loop meant that regex parsing on hundreds of millions of FASTA headers throttled the pipe,
causing a single Kraken 2 thread to block inside `#pragma omp critical(seqread)` while all
other Kraken 2 threads sat idle (resulting in 1.15 average CPU core utilization).

This script eliminates that bottleneck:
1. **Direct C-Memory Slicing at 200 bp Step (2x Coverage)**:
   Slices each contig directly in memory (`seq[s0:e0]`) with step `S = 200 bp` and window
   `W = 400 bp` (plus the exact contig-end window `[L - W, L)` when `(L - W) % S != 0` and
   `L >= W`, or `[0, L)` when `L < W`), emitting `>{t_id}|{contig}|{s0}|{e0}` headers
   directly without intermediate string formatting or `awk` regex parsing.
2. **Multi-Core Producer Pool**:
   Uses a pool of worker processes (default: 8 workers, running on the non-Kraken cores of
   the 72-core `largemem` node) to slice contigs in parallel and stream large multi-megabyte
   buffers into the FIFO/pipe so Kraken 2's 3 MB `LoadBlock()` reads return instantaneously.
"""

import argparse
import csv
import fcntl
import multiprocessing as mp
import os
import sys
from pathlib import Path
from typing import Iterator, List, Tuple

# Increase Linux pipe/FIFO buffer size where supported (up to 1 MB default unprivileged max)
F_SETPIPE_SZ = 1031


def iter_fasta_records(fasta_path: Path) -> Iterator[Tuple[str, str]]:
    """Yield `(contig_id, sequence_str)` from an uncompressed FASTA file."""
    header = None
    chunks: List[str] = []
    with open(fasta_path, "r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == ">":
                if header is not None:
                    yield header, "".join(chunks)
                # First whitespace-delimited token after '>'
                header = line[1:].strip().split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def tile_contig_chunk(args: Tuple[str, str, str, int, int]) -> bytes:
    """
    Slice a single contig `(tid_str, contig_id, seq, window_bp, step_bp)` into
    FASTA records `>{tid}|{contig}|{s0}|{e0}\\n{subseq}\\n` and return encoded bytes.
    """
    tid_str, contig_id, seq, window_bp, step_bp = args
    seq_len = len(seq)
    if seq_len == 0:
        return b""

    prefix = f">{tid_str}|{contig_id}|"
    out_lines: List[str] = []

    if seq_len <= window_bp:
        out_lines.append(f"{prefix}0|{seq_len}\n{seq}\n")
        return "".join(out_lines).encode("ascii")

    last_full_start = seq_len - window_bp
    for s0 in range(0, last_full_start + 1, step_bp):
        e0 = s0 + window_bp
        out_lines.append(f"{prefix}{s0}|{e0}\n{seq[s0:e0]}\n")

    # Match `seqkit sliding -g` greedy terminal window if contig end wasn't hit exactly
    if (last_full_start % step_bp) != 0:
        out_lines.append(
            f"{prefix}{last_full_start}|{seq_len}\n{seq[last_full_start:seq_len]}\n"
        )

    return "".join(out_lines).encode("ascii")


def task_generator(
    active_species_tsv: Path,
    window_bp: int,
    step_bp: int,
    max_contig_slice_bp: int = 20_000_000,
) -> Iterator[Tuple[str, str, str, int, int]]:
    """
    Yield contig (or sub-contig) slicing tasks across all active species.
    Large chromosomes (> 20 Mbp) are split along exact `step_bp`-aligned boundaries
    so worker processes stay memory-light and load-balanced.
    """
    with open(active_species_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tid_str = str(row["t_id"]).strip()
            fasta_str = str(row["filtered_fasta"]).strip()
            if not tid_str or not fasta_str:
                continue
            fasta_path = Path(fasta_str)
            if not fasta_path.is_file():
                raise FileNotFoundError(
                    f"Filtered reference FASTA not found for TaxID {tid_str}: {fasta_path}"
                )
            for contig_id, seq in iter_fasta_records(fasta_path):
                seq_len = len(seq)
                if seq_len <= max_contig_slice_bp:
                    yield (tid_str, contig_id, seq, window_bp, step_bp)
                else:
                    # Split chromosome into step-aligned chunks with window overlap
                    # Wait: to keep global coordinates s0/e0 exact, pass sub-sequence with offset!
                    yield (tid_str, contig_id, seq, window_bp, step_bp)


def tile_species_or_batch(
    args: Tuple[str, List[Tuple[str, str, int, int, bool]], int, int]
) -> bytes:
    """
    Slice a batch of contigs/segments `(tid_str, [(contig_id, sub_seq, offset0, orig_len, is_last), ...], window_bp, step_bp)`
    into exact global 0-based coordinates `[s0, e0)` in a single ~10 MB buffer.
    """
    tid_str, segments, window_bp, step_bp = args
    out_lines: List[str] = []
    append = out_lines.append

    for contig_id, sub_seq, offset0, orig_len, is_last in segments:
        sub_len = len(sub_seq)
        if sub_len == 0:
            continue
        prefix = f">{tid_str}|{contig_id}|"
        if orig_len <= window_bp:
            if offset0 == 0:
                append(f"{prefix}0|{orig_len}\n{sub_seq}\n")
            continue

        last_local_start = sub_len - window_bp
        for local_s0 in range(0, last_local_start + 1, step_bp):
            s0 = offset0 + local_s0
            append(f"{prefix}{s0}|{s0 + window_bp}\n{sub_seq[local_s0:local_s0 + window_bp]}\n")

        if is_last:
            last_global_start = orig_len - window_bp
            if (last_global_start % step_bp) != 0:
                local_tail = last_global_start - offset0
                if 0 <= local_tail < sub_len:
                    append(f"{prefix}{last_global_start}|{orig_len}\n{sub_seq[local_tail:]}\n")

    return "".join(out_lines).encode("ascii")


def batched_task_generator(
    active_species_tsv: Path,
    window_bp: int,
    step_bp: int,
    target_batch_bp: int = 5_000_000,
) -> Iterator[Tuple[str, List[Tuple[str, str, int, int, bool]], int, int]]:
    """
    Yield ~5 Mbp batches of contigs/segments across all active species.
    Groups thousands of short draft contigs into a single 5 Mbp IPC batch while
    splitting multi-hundred-megabyte chromosomes into 5 Mbp step-aligned slices.
    """
    stride_bp = (target_batch_bp // step_bp) * step_bp
    with open(active_species_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tid_str = str(row["t_id"]).strip()
            fasta_str = str(row["filtered_fasta"]).strip()
            if not tid_str or not fasta_str:
                continue
            fasta_path = Path(fasta_str)
            if not fasta_path.is_file():
                raise FileNotFoundError(
                    f"Filtered reference FASTA not found for TaxID {tid_str}: {fasta_path}"
                )
            batch: List[Tuple[str, str, int, int, bool]] = []
            batch_bp = 0
            for contig_id, seq in iter_fasta_records(fasta_path):
                orig_len = len(seq)
                if orig_len == 0:
                    continue
                if orig_len <= stride_bp + window_bp:
                    batch.append((contig_id, seq, 0, orig_len, True))
                    batch_bp += orig_len
                    if batch_bp >= stride_bp:
                        yield (tid_str, batch, window_bp, step_bp)
                        batch = []
                        batch_bp = 0
                else:
                    if batch:
                        yield (tid_str, batch, window_bp, step_bp)
                        batch = []
                        batch_bp = 0
                    last_full_start = orig_len - window_bp
                    offset0 = 0
                    while offset0 <= last_full_start:
                        next_offset = offset0 + stride_bp
                        if next_offset > last_full_start:
                            yield (
                                tid_str,
                                [(contig_id, seq[offset0:], offset0, orig_len, True)],
                                window_bp,
                                step_bp,
                            )
                            break
                        else:
                            end_slice = (next_offset - step_bp) + window_bp
                            yield (
                                tid_str,
                                [(contig_id, seq[offset0:end_slice], offset0, orig_len, False)],
                                window_bp,
                                step_bp,
                            )
                            offset0 = next_offset
            if batch:
                yield (tid_str, batch, window_bp, step_bp)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stream 2x in silico tiled FASTA windows across all active species in parallel."
    )
    parser.add_argument("--active-species", required=True, help="Path to active_species.tsv")
    parser.add_argument("--window-bp", type=int, default=400, help="Sliding window size in bp (default: 400)")
    parser.add_argument("--step-bp", type=int, default=200, help="Sliding step size in bp (default: 200)")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel producer worker processes (default: 8)")
    args = parser.parse_args()

    out_fd = sys.stdout.fileno()
    try:
        fcntl.fcntl(out_fd, F_SETPIPE_SZ, 1048576)
    except OSError:
        pass

    out_buf = sys.stdout.buffer
    tasks = batched_task_generator(Path(args.active_species), args.window_bp, args.step_bp)

    if args.workers <= 1:
        for t in tasks:
            data = tile_species_or_batch(t)
            if data:
                out_buf.write(data)
    else:
        with mp.Pool(processes=args.workers) as pool:
            for data in pool.imap(tile_species_or_batch, tasks, chunksize=1):
                if data:
                    out_buf.write(data)

    out_buf.flush()


if __name__ == "__main__":
    main()
