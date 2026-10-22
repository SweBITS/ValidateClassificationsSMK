# Stage 2 (`stream_kraken_cr`) Optimization, Core-Scaling Benchmark & CR Comparison Plan

## 1. Objectives

When re-running the pipeline from Stage 2 (`stream_kraken_cr`), we want to evaluate three things empirically:
1. **End-to-End Stage 2 Speedup (`Baseline 4x / 72t` vs. `New 2x / 36t`)**:
   - Measure the walltime and throughput (`Kseq/min`, `Mbp/min`, CPU %, context switches) gained by:
     - Switching from `100 bp` step ($4\times$ tiling coverage, `~242.7M` windows) to **`200 bp` step ($2\times$ tiling coverage, `~121.3M` windows)**.
     - Replacing the single-threaded `seqkit sliding | awk match(...)` FIFO loop (which blocked Kraken 2 inside `#pragma omp critical(seqread)` and limited CPU utilization to `115.11%`) with the multi-worker parallel byte slicer ([`workflow/scripts/stream_tiled_fasta.py`](workflow/scripts/stream_tiled_fasta.py), `8` workers, `1 MB` pipe buffer).
     - Loading `intel/2025a`, exporting `OMP_NUM_THREADS=36`, and running Kraken 2 with **`36` OpenMP threads** (2 sockets on the 4-socket × 18-core `largemem` node).
2. **Controlled Kraken 2 Core-Scaling Curve (`6, 12, 18, 24, 36, 72` Cores)**:
   - Immediately after the main production stream finishes inside `stream_kraken_cr`, the entire 2.6 TB database (`hash.k2d`) is warm in the Linux kernel page cache on the `largemem` node.
   - [`workflow/scripts/benchmark_kraken_threads.py`](workflow/scripts/benchmark_kraken_threads.py) automatically generates a deterministic `5,000,000`-window (`2.0 Gbp`) slice in `$SNIC_TMP` and benchmarks Kraken 2 with `--memory-mapping` across **`6, 12, 18, 24, 36, 72` cores**:
     - `6` cores (1/3 socket)
     - `12` cores (2/3 socket)
     - `18` cores (1 full Broadwell socket / NUMA node)
     - `24` cores (1.33 sockets — matching earlier `01-CLASSIFICATION` runs)
     - `36` cores (2 full Broadwell sockets)
     - `72` cores (all 4 Broadwell sockets — testing cross-socket NUMA & OpenMP lock contention)
   - Outputs are saved to `benchmarks/03_stream_kraken_cr/kraken2_thread_scaling.tsv` and logged in `logs/03_stream_kraken_cr/kraken2_thread_scaling.log`.
3. **Classifiable Region (`cr_2x.bed`) & Downstream Stage 3–5 Concordance**:
   - Compare every species' new `{t_id}.cr_2x.bed` against the backed-up baseline (`results/02_classifiable_regions/_baseline_4x_72t/beds/{t_id}.cr_2x.bed`):
     - Interior `400 bp` windows stepped every `200 bp` (`0, 200, 400, ...`) are mathematically identical to the even-indexed (`idx % 2 == 0`) windows of the `100 bp` step.
     - Only the greedy terminal window at the right tip of a contig (`[L - 400, L)` when `(L - 400) % 200 != 0`) can differ by `< 200 bp`.
   - Compare post-Stage-2 validation summaries (`results/05_reports/` vs. `results/05_reports/_baseline_4x_72t/`) to verify whether any species' Pearson $R$, pruned bin count, or validation verdict changes.

---

## 2. Archived Baseline Snapshot Locations

All pre-rerun baseline artifacts (from the `100 bp` step / `72`-thread run) have been backed up on disk:

| Baseline Artifact | Backup Path | Key Baseline Value |
| :--- | :--- | :--- |
| **Snakemake Benchmark** | `results/02_classifiable_regions/_baseline_4x_72t/benchmarks/stream_kraken_cr.tsv` | `40,518.75 s` (`11:15:18`), `mean_load = 115.11%`, `max_rss = 2,678,414 MB` |
| **Kraken 2 Stream Log** | `results/02_classifiable_regions/_baseline_4x_72t/logs/kraken2_stream.log` | `242,682,268` seqs (`97,065.21 Mbp`) in `36,929.33 s` (`394.3 Kseq/m`, `157.70 Mbp/m`) |
| **CR Consumer Log** | `results/02_classifiable_regions/_baseline_4x_72t/logs/consume_kraken_stream_to_cr.log` | `240,101,293` classified windows drained |
| **CR Summary Table** | `results/02_classifiable_regions/_baseline_4x_72t/cr_streaming_summary.tsv` | 44 active species (`cr_bp_2x`, `cr_segments_2x`) |
| **Per-Species CR BEDs** | `results/02_classifiable_regions/_baseline_4x_72t/beds/{t_id}.cr_2x.bed` | All 44 baseline `cr_2x.bed` (and `cr_4x.bed`) files |
| **Stage 5 Reports** | `results/05_reports/_baseline_4x_72t/bins{50,100}_*/` | All 6 completed `bins100_*` and `bins50_*` report directories |

---

## 3. How to Launch the Rerun

To force Snakemake to re-run Stage 2 (`stream_kraken_cr`, including the `6, 12, 18, 24, 36, 72`-core benchmark sweep) and automatically propagate through Stages 3–5 for the active config profile:

```bash
# Force Stage 2 (`stream_kraken_cr`) and all downstream Stage 3–5 rules to re-run
# (Reuses Stage 0, Stage 1 filtered references, and Stage 4 extracted sample FASTQs)
snakemake --profile slurm/hpc2n --use-conda --forcerun stream_kraken_cr
```

*(Note: After this one-time benchmark run is complete, if you don't want future Stage 2 runs to perform the 6-point thread sweep, set `benchmark_threads_list: ""` under `resources.stream_kraken_cr` in `config/config.yaml`.)*

---

## 4. How to Evaluate & Compare Results Once Finished

As soon as the Snakemake run completes, run:

```bash
python3 workflow/scripts/compare_stage2_cr_and_benchmarks.py
```

This script automatically reads the baseline backup (`results/02_classifiable_regions/_baseline_4x_72t/`), the new Stage 2 logs/benchmarks, the core-scaling table (`benchmarks/03_stream_kraken_cr/kraken2_thread_scaling.tsv`), and all 44 `{t_id}.cr_2x.bed` files, and generates:
1. **`results/02_classifiable_regions/stage2_comparison_report.md`**:
   - **Section 1**: Side-by-side table of the production Stage 2 run (`Baseline 4x/72t` vs. `New 2x/36t`).
   - **Section 2**: Complete `6, 12, 18, 24, 36, 72`-core scaling table (`Kseq/min`, `Mbp/min`, speedup vs. 6 cores, parallel efficiency %, CPU %, and context switches) to pick the permanent optimal `kraken_threads` value.
   - **Section 3**: Genome-wide and per-species `cr_2x.bed` concordance table (`Old CR bp`, `New CR bp`, `Net Diff bp`, `Old Segs`, `New Segs`, and exact base-pair `Jaccard %`).
2. **`results/02_classifiable_regions/stage2_cr_bed_comparison.tsv`**:
   - Machine-readable TSV of the 44-species CR coordinate comparison.
