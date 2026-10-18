#!/usr/bin/env python3
"""
build_cr_bins.py — Stage 3: Shortest-to-Longest Contig Sorting & Linear Variable-Size CR Binning
================================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
Testing whether mapped eDNA reads distribute uniformly across a species' Classifiable Regions
(CRs) requires regressing Mapped Read Count ($Y$) against Genomic Bin CR Length ($X$). Two
naive binning approaches fail in practice:
  1. **Whole-contig binning** fails on chromosome-scale assemblies (only a few points; a single
     macrochromosome dominates OLS leverage) and on highly fragmented draft assemblies (thousands
     of tiny contigs with extreme Poisson variance).
  2. **Equal-sized bins** ($X_1 = X_2 = \\dots = X_{100} = \\bar{L}$) collapse all $X$ values onto
     a single vertical line ($\text{Var}(X) = 0$), making linear correlation undefined.

OUR ALGORITHM (DETERMINISTIC LINEAR VARIABLE-SIZE BINNING):
-----------------------------------------------------------
1. **Shortest-to-Longest Contig Ordering**:
   We take all contigs from the filtered reference assembly that contain at least one merged CR
   and sort those contigs by **total contig length ascending** (`(contig_length_bp, contig_id)`).
   Within each contig, CR intervals are ordered linearly from $5' \\to 3'$ (`start0 < end0`).

2. **Linearly Increasing Bin Target Capacities ($0.3\\bar{L}$ to $1.7\\bar{L}$)**:
   Let $L_{\\text{CR}}$ be the total classifiable sequence length across the genome and
   $\\bar{L} = L_{\\text{CR}} / B$ be the mean CR length per bin (default $B = 100$).
   Each sequential bin $k \\in \\{1, 2, \\dots, B\\}$ receives target length:
     $$W_k = \\bar{L} \\times \\left(0.3 + 1.4 \\times \\frac{k - 1}{B - 1}\\right)$$
   Because $(0.3 + 1.7) / 2 = 1.0$, $\\sum_{k=1}^B W_k = L_{\\text{CR}}$. Using cumulative integer
   rounding $C_k = \\text{round}(\\sum_{j=1}^k W_j)$ with $\\text{target\\_bp}_k = C_k - C_{k-1}$
   guarantees that integer bin capacities sum to $L_{\\text{CR}}$ to the **exact single base pair**.

3. **Linear Walk & Exact Boundary Splitting**:
   Walking along the sorted CR intervals, we fill Bin 1 ($0.3\\bar{L}$), then Bin 2, ..., up to
   Bin $B$ ($1.7\\bar{L}$). Whenever a continuous CR interval $[s, e)$ crosses the remaining
   capacity boundary of Bin $k$, we split that interval cleanly at coordinate $s + \\text{rem}$
   so Bin $k$ reaches its exact target capacity and the remaining bases start Bin $k+1$.

BIOLOGICAL & STATISTICAL BENEFITS:
  - Small draft contigs naturally coalesce into the smaller initial bins ($0.3\\bar{L} \\dots$),
    while large chromosomal CRs remain in contiguous linear neighborhoods inside the larger bins
    ($\\dots 1.7\\bar{L}$).
  - Produces a clean $5.7\\times$ horizontal spread on the X-axis with **zero random numbers**.
  - When Stage 5 prunes an outlier bin (e.g., an rRNA tandem array, integrated viral element, or
    mitochondrial contig), it excises a contiguous genomic neighborhood rather than scattered loci.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import setup_logger


def load_contig_sequences(filtered_fasta: Path) -> Dict[str, str]:
    """Load contig sequences from the Stage 1 filtered FASTA into memory."""
    seqs: Dict[str, str] = {}
    if not filtered_fasta.is_file():
        return seqs

    curr_id = ""
    curr_chunks: List[str] = []
    with open(filtered_fasta, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if curr_id:
                    seqs[curr_id] = "".join(curr_chunks)
                curr_id = line[1:].strip().split()[0]
                curr_chunks = []
            else:
                curr_chunks.append(line.strip())
        if curr_id:
            seqs[curr_id] = "".join(curr_chunks)
    return seqs


def load_cr_bed(bed_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Load merged 0-based half-open CR intervals `[start, end)` per contig from BED."""
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    if not bed_path.is_file():
        return by_contig

    with open(bed_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                contig = parts[0]
                s = int(parts[1])
                e = int(parts[2])
                if e > s:
                    by_contig.setdefault(contig, []).append((s, e))

    for contig in by_contig:
        by_contig[contig].sort(key=lambda x: (x[0], x[1]))
    return by_contig


def compute_exact_bin_capacities(
    total_cr_bp: int,
    n_bins: int,
    min_scale: float = 0.3,
    max_scale: float = 1.7,
) -> List[int]:
    """
    Compute exact integer target lengths $W_1 \\dots W_B$ summing to `total_cr_bp` via
    cumulative rounding of the linear progression from `min_scale * L_bar` to `max_scale * L_bar`.
    """
    if n_bins <= 1:
        return [total_cr_bp]

    mean_len = total_cr_bp / float(n_bins)
    slope = max_scale - min_scale
    raw_weights = [
        mean_len * (min_scale + slope * ((k - 1) / float(n_bins - 1)))
        for k in range(1, n_bins + 1)
    ]
    # Normalize raw_weights just in case (min_scale + max_scale) / 2 != 1.0 in custom user config
    weight_sum = sum(raw_weights)
    scale_factor = (total_cr_bp / weight_sum) if weight_sum > 0 else 0.0

    capacities: List[int] = []
    prev_cum = 0
    running_float = 0.0
    for k, w in enumerate(raw_weights, start=1):
        running_float += w * scale_factor
        curr_cum = total_cr_bp if k == n_bins else int(round(running_float))
        cap = max(0, curr_cum - prev_cum)
        capacities.append(cap)
        prev_cum = curr_cum

    return capacities


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 3: Sort contigs shortest-to-longest and partition CRs into 100 linearly increasing bins."
    )
    parser.add_argument("--taxid", type=int, required=True, help="Target species TaxID")
    parser.add_argument("--filtered-fasta", required=True, help="Path to Stage 1 filtered reference FASTA")
    parser.add_argument("--contig-lengths", required=True, help="Path to Stage 1 contig_lengths.tsv")
    parser.add_argument("--cr-bed", required=True, help="Path to active-tier CR BED file (e.g. {t_id}.cr_2x.bed)")
    parser.add_argument("--out-cr-fasta", required=True, help="Output FASTA: {t_id}.cr_contigs.fa")
    parser.add_argument("--out-bins-tsv", required=True, help="Output TSV: {t_id}.bins_metadata.tsv")
    parser.add_argument("--n-bins", type=int, default=100, help="Number of bins B (default: 100)")
    parser.add_argument("--min-scale", type=float, default=0.3, help="Smallest bin scale factor (default: 0.3)")
    parser.add_argument("--max-scale", type=float, default=1.7, help="Largest bin scale factor (default: 1.7)")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger(f"build_cr_bins_{args.taxid}", args.log)

    # 1. Load contig lengths and CR BED intervals
    contig_lengths: Dict[str, int] = {}
    with open(args.contig_lengths, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            contig_lengths[row["contig_id"]] = int(row["length_bp"])

    cr_by_contig = load_cr_bed(Path(args.cr_bed))
    total_cr_bp = sum(e - s for segs in cr_by_contig.values() for s, e in segs)

    out_fasta_path = Path(args.out_cr_fasta)
    out_bins_path = Path(args.out_bins_tsv)
    out_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    out_bins_path.parent.mkdir(parents=True, exist_ok=True)

    # Graceful handling if a species has insufficient CR sequence to partition into n_bins
    if total_cr_bp < args.n_bins:
        logger.warning(
            f"TaxID {args.taxid} has only {total_cr_bp} bp of CR sequence (< {args.n_bins} bins). "
            f"Writing empty CR FASTA and zero-capacity bin table for graceful downstream exclusion."
        )
        out_fasta_path.write_text("", encoding="utf-8")
        with open(out_bins_path, "w", encoding="utf-8", newline="") as bfh:
            writer = csv.writer(bfh, delimiter="\t")
            writer.writerow([
                "bin_id", "bin_index", "target_cr_bp", "actual_cr_bp",
                "num_cr_segments", "num_distinct_contigs", "source_contigs_summary",
            ])
            for k in range(1, args.n_bins + 1):
                writer.writerow([f"bin_{k:03d}", k, 0, 0, 0, 0, "NONE"])
        return

    # 2. Sort contigs containing >= 1 CR from shortest total contig length to longest,
    #    breaking ties deterministically by contig_id
    contigs_with_cr = sorted(
        cr_by_contig.keys(),
        key=lambda cid: (contig_lengths.get(cid, 0), cid),
    )
    logger.info(
        f"TaxID {args.taxid}: {len(contigs_with_cr):,} contigs contain CRs "
        f"(total L_CR = {total_cr_bp:,} bp; shortest contig = {contig_lengths.get(contigs_with_cr[0], 0):,} bp, "
        f"longest contig = {contig_lengths.get(contigs_with_cr[-1], 0):,} bp)."
    )

    # Flatten into an ordered queue of (contig_id, start0, end0)
    ordered_intervals: List[Tuple[str, int, int]] = []
    for cid in contigs_with_cr:
        for s, e in cr_by_contig[cid]:
            ordered_intervals.append((cid, s, e))

    # 3. Compute exact integer target capacities W_1 .. W_B
    capacities = compute_exact_bin_capacities(
        total_cr_bp=total_cr_bp,
        n_bins=args.n_bins,
        min_scale=args.min_scale,
        max_scale=args.max_scale,
    )
    assert sum(capacities) == total_cr_bp, f"Capacity sum {sum(capacities)} != {total_cr_bp}"

    # 4. Walk linearly along ordered CR intervals, splitting cleanly at bin capacity boundaries
    # Each bin stores a list of (contig_id, seg_start0, seg_end0)
    binned_segments: List[List[Tuple[str, int, int]]] = [[] for _ in range(args.n_bins)]
    bin_idx = 0
    rem_cap = capacities[0]

    for cid, s, e in ordered_intervals:
        curr_s = s
        while curr_s < e:
            # Advance to next bin with non-zero remaining capacity
            while bin_idx < args.n_bins and rem_cap == 0:
                bin_idx += 1
                if bin_idx < args.n_bins:
                    rem_cap = capacities[bin_idx]

            if bin_idx >= args.n_bins:
                break

            avail = e - curr_s
            take = min(avail, rem_cap)
            seg_end = curr_s + take
            binned_segments[bin_idx].append((cid, curr_s, seg_end))
            curr_s = seg_end
            rem_cap -= take

    # 5. Load reference sequences and write {t_id}.cr_contigs.fa + {t_id}.bins_metadata.tsv
    seqs = load_contig_sequences(Path(args.filtered_fasta))

    global_seg_counter = 0
    bins_meta_rows: List[Dict[str, object]] = []

    with open(out_fasta_path, "w", encoding="utf-8") as fa_fh:
        for k_zero, seg_list in enumerate(binned_segments):
            k_one = k_zero + 1
            bin_id = f"bin_{k_one:03d}"
            target_bp = capacities[k_zero]
            actual_bp = 0
            contig_bp_tally: Dict[str, int] = {}

            for cid, s, e in seg_list:
                global_seg_counter += 1
                cr_id = f"cr_{global_seg_counter:06d}"
                seg_len = e - s
                actual_bp += seg_len
                contig_bp_tally[cid] = contig_bp_tally.get(cid, 0) + seg_len

                subseq = seqs[cid][s:e]
                if len(subseq) != seg_len:
                    raise RuntimeError(
                        f"Sequence slice length mismatch on {cid}:{s}-{e}: expected {seg_len}, got {len(subseq)}"
                    )

                # Encode cr_id, bin_id, source_contig, and 0-based start-end into the FASTA header
                header = f">{cr_id}|bin={bin_id}|contig={cid}|coords={s}-{e}"
                fa_fh.write(header + "\n")
                for i in range(0, seg_len, 80):
                    fa_fh.write(subseq[i : i + 80] + "\n")

            # Summarize source contigs for tabular metadata (up to first 5 contigs + count if more)
            distinct_contigs = list(contig_bp_tally.keys())
            summary_parts = [f"{c}:{contig_bp_tally[c]}bp" for c in distinct_contigs[:5]]
            if len(distinct_contigs) > 5:
                summary_parts.append(f"...(+{len(distinct_contigs) - 5}_more_contigs)")
            contigs_summary = ";".join(summary_parts) if summary_parts else "NONE"

            bins_meta_rows.append({
                "bin_id": bin_id,
                "bin_index": k_one,
                "target_cr_bp": target_bp,
                "actual_cr_bp": actual_bp,
                "num_cr_segments": len(seg_list),
                "num_distinct_contigs": len(distinct_contigs),
                "source_contigs_summary": contigs_summary,
            })

    with open(out_bins_path, "w", encoding="utf-8", newline="") as bfh:
        fields = [
            "bin_id", "bin_index", "target_cr_bp", "actual_cr_bp",
            "num_cr_segments", "num_distinct_contigs", "source_contigs_summary",
        ]
        writer = csv.DictWriter(bfh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(bins_meta_rows)

    logger.info(
        f"TaxID {args.taxid}: Wrote {global_seg_counter:,} CR segments across {args.n_bins} bins "
        f"(bin_001 = {capacities[0]:,} bp .. bin_{args.n_bins:03d} = {capacities[-1]:,} bp) "
        f"to {out_fasta_path.name}."
    )


if __name__ == "__main__":
    main()
