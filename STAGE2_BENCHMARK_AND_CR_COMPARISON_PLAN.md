# Stage 2 (`stream_kraken_cr`) Optimization, Core-Scaling Benchmark & CR Concordance Report

## 1. Executive Summary

Stage 2 (`stream_kraken_cr` in [`workflow/rules/03_stream_kraken_cr.smk`](workflow/rules/03_stream_kraken_cr.smk)) was refactored and benchmarked on Kebnekaise (`cpu_largemem`, 4-socket × 18-core Intel Xeon E7-8860v4 Broadwell, 3.0 TB RAM) with three improvements:
1. **Direct $2\times$ Tiling (`200 bp` step, `400 bp` window)**:
   - Replaced `100 bp` step ($4\times$ coverage, `242.7M` windows) with **`200 bp` step ($2\times$ coverage, `121.5M` windows)**, halving the number of windows while preserving identical interior window coordinates (`0, 200, 400, ...`).
2. **Parallel Multi-Worker Python Tiler ([`workflow/scripts/stream_tiled_fasta.py`](workflow/scripts/stream_tiled_fasta.py))**:
   - Replaced the single-threaded `seqkit sliding | awk match(...)` FIFO loop (which bottlenecked on regex header parsing and starved Kraken 2 inside `#pragma omp critical(seqread)`) with an 8-worker Python C-memory byte slicer feeding a `1 MB` POSIX pipe buffer.
3. **Optimal Toolchain & OpenMP Thread Count (`intel/2025a`, `kraken_threads: 36`)**:
   - Loaded `intel/2025a` (`GCCcore/14.2.0` runtime matching `kraken2_SweBiTS_MHG_intel`) and set `OMP_NUM_THREADS=36` (`2` of `4` NUMA sockets), leaving the other `36` cores on the 72-core node free for the parallel Python tiler, `awk`, and CR interval consumer ([`workflow/scripts/consume_kraken_stream_to_cr.py`](workflow/scripts/consume_kraken_stream_to_cr.py)).

---

## 2. Production Stage 2 Stream Benchmark (`Baseline 4x/72t` vs. `New 2x/36t`)

| Metric | Baseline (`100bp` step, `72` threads, `seqkit\|awk`) | New (`200bp` step, `36` threads, `stream_tiled_fasta.py`, `intel/2025a`) | Speedup / Change |
| :--- | ---: | ---: | ---: |
| **Simulated Windows Classified** | `242,682,268` | `121,458,409` | `0.50x` |
| **Total Bases Classified (`Mbp`)** | `97,065.21` | `48,583.36` | `0.50x` |
| **Kraken 2 Classify Time (`s`)** | `36,929.3 s` (`10.26 h`) | `6,327.2 s` (`1.76 h`) | **`5.84x` faster** |
| **Per-Window Throughput (`Kseq/min`)** | `394.3` | `1,151.8` | **`2.92x` higher** |
| **Base-Pair Throughput (`Mbp/min`)** | `157.70` | `460.71` | **`2.92x` higher** |
| **Kraken 2 Walltime (incl. ~63m Lustre DB Load)** | `11:15:18` (`40,518.8 s`) | `2:53:11` (`10,391.0 s`) | **`3.90x` faster** |
| **Mean CPU Load (`mean_load`)** | `115.11%` | `186.83%` | `+62.3%` |
| **Peak RSS (`max_rss` MB)** | `2,678,414.9` | `2,704,419.4` | `+0.97%` |

---

## 3. Controlled 44-Species Kraken 2 Core-Scaling Sweep (`6, 12, 18, 24, 36, 72` Cores)

To determine the optimal `kraken_threads` setting on Kebnekaise's 72-core (`4` sockets × `18` cores/socket) Broadwell `largemem` nodes (`SLURM Job 48083443` on `b-cn1050`), the `2.68 TB` `hash.k2d` database was loaded into heap RAM once and tested across **`6, 12, 18, 24, 36, 72` OpenMP threads** on the exact same **`2,500,002`-window (`1,000.00 Mbp`)** dataset stratified and interleaved across all **44 active species** (`~56,818` windows = `22.7 Mbp` per species) stored in local `/scratch` RAM disk:

| Threads | NUMA Footprint (Broadwell 4×18c) | Classify Wall Time (`s`) | Throughput (`Kseq/min`) | Throughput (`Mbp/min`) | Speedup vs. 6c | Speedup vs. Warm 12c | CPU Utilization (`%` / Active Cores) | Per-Core Efficiency (`CPU % / Threads`) | User / Sys CPU (`s`) | Voluntary / Involuntary Context Switches |
| ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **6** *(cold 1st pass)* | 0.33 socket | `23.64 s` | `6,345.9` | `2,538.34` | `1.00x` | `0.43x` | `590.8%` (`5.91c`) | **98.5%** | `115.3s` / `24.3s` | `288` / `299` |
| **12** | 0.67 socket | `10.28 s` | `14,595.7` | `5,838.27` | `2.30x` | `1.00x` | `1,179.1%` (`11.79c`) | **98.3%** | `118.6s` / `2.6s` | `205` / `266` |
| **18** | **1 full socket** | `7.18 s` | `20,882.8` | `8,353.13` | `3.29x` | `1.43x` *(95.4% eff.)* | `1,737.5%` (`17.38c`) | **96.5%** | `122.6s` / `2.2s` | `321` / `263` |
| **24** | 1.33 sockets | `5.60 s` | `26,768.6` | `10,707.45` | `4.22x` | `1.83x` *(91.7% eff.)* | `2,283.8%` (`22.84c`) | **95.2%** | `125.9s` / `2.1s` | `466` / `255` |
| **36** | **2 full sockets (50% node)** | **`3.89 s`** | **`38,605.0`** | **`15,441.99`** | **`6.08x`** | **`2.64x`** *(88.2% eff.)* | **`3,318.9%` (`33.19c`)** | **92.2%** | `126.9s` / `2.1s` | `912` / `285` |
| **72** | **All 4 sockets (100% node)** | `3.08 s` | `48,690.5` | `19,476.18` | `7.67x` | `3.34x` *(55.6% eff.)* | `4,488.1%` (`44.88c`) | **62.3%** | `134.7s` / `3.6s` | **`4,963` / `1,005`** |

*(Raw TSV saved in [`benchmarks/03_stream_kraken_cr/kraken2_thread_scaling.tsv`](benchmarks/03_stream_kraken_cr/kraken2_thread_scaling.tsv).)*

### Why `36` Cores Is the Optimal Setting (`config/config.yaml`)
1. **High Parallel Efficiency Up to 2 Sockets (`36` Cores)**:
   - Scaling from `12` to `36` cores (`3.0x` cores) yields a **`2.64x` speedup (`88.2%` scaling efficiency)** and maintains **`92.2%` per-core CPU utilization (`33.19` active cores out of `36`)**.
2. **Severe Diminishing Returns & Lock Contention at 4 Sockets (`72` Cores)**:
   - Doubling threads from `36` to `72` (`+100%` cores) increases throughput by only **`+26.1%`** (`3.89 s` $\rightarrow$ `3.08 s` per Gbp), while per-core CPU efficiency drops from **`92.2%` to `62.3%`** (`27+` cores idle/spinning) and voluntary context switches jump **`5.4x`** (`912` $\rightarrow$ `4,963`) due to cross-socket QPI memory-controller contention across the 2.68 TB hash table and OpenMP `#pragma omp critical` serialization (`seqread` and `output_lock`).
   - Allocating `36` threads to Kraken 2 leaves the remaining `36` cores on the 72-core node dedicated to the 8-worker Python FIFO producer (`stream_tiled_fasta.py`), `awk`, and `consume_kraken_stream_to_cr.py`.
3. **Why Different Genomes Classify at Different Speeds (`ResolveTree` Hotspot)**:
   - In Kraken 2 (`classify.cc:442-483`), `ResolveTree()` runs a nested ancestor-walk loop over all distinct TaxIDs hit within a 400 bp window. Species-unique windows hitting 1–3 taxa resolve in nanoseconds, whereas repeat-rich eukaryotic scaffolds (*Orobanche minor*, *Salmo trutta*, *Columba*) hit dozens of distinct taxa across the tree per window, making `ResolveTree()` 50×–200× more CPU-intensive on repeat-heavy regions.
4. **Why `--memory-mapping` Cannot Be Used Against Lustre on Kebnekaise Without `large_shm`**:
   - Kebnekaise's Lustre client (`llite.*.max_cached_mb`) caps per-node Lustre page caching at `~24.9 GiB` even on `3.0 TB` `largemem` nodes. Passing `--memory-mapping` (`mmap(MAP_PRIVATE)` on a Lustre file descriptor) causes constant Lustre page-cache thrashing (`majflt > 60,000`) when random hash lookups span `2.6 TB`. Standard heap loading (without `--memory-mapping`) reads `hash.k2d` sequentially at `~1.5 GB/s` into `2.68 TB` of anonymous RAM once and runs with zero disk I/O thereafter.

---

## 4. Classifiable Region (`cr_2x.bed`) & Downstream Stage 3–5 Concordance

Full per-species comparison tables are stored in [`results/02_classifiable_regions/stage2_comparison_report.md`](results/02_classifiable_regions/stage2_comparison_report.md) and [`results/02_classifiable_regions/stage2_cr_bed_comparison.tsv`](results/02_classifiable_regions/stage2_cr_bed_comparison.tsv):
- **Old Total `cr_2x.bed` Base Pairs**: `12,450,502,358 bp`
- **New Total `cr_2x.bed` Base Pairs**: `12,451,033,634 bp` (`+531,276 bp`, `+0.0043%`, arising solely from using full `400 bp` terminal windows `[L - 400, L)` at contig right tips where `(L - 400) % 200 != 0`)
- **Genome-Wide Base-Pair Jaccard Identity Across All 44 Species**: **`99.989797%`** (`12,450,132,801 / 12,451,403,191 bp`)
- **Downstream Stage 3–5 Concordance (`bins50_s0p3-1p7_mapq20_aln80_id95`)**: `49 / 50` targets (`43 / 44` active species) retained the exact same `validation_status` (only borderline moss `127167` *Orthodicranum fulvum*, which previously pruned `27 / 50` bins from initial `R = 0.28`, shifted across the `50%` retention floor from `VALIDATED_PRUNED` to `FAILED_LOW_CORRELATION`).

---

## 5. Archived Baseline Snapshot Locations

All pre-rerun baseline artifacts (from the `100 bp` step / `72`-thread run) are preserved on disk for reference:

| Baseline Artifact | Backup Path | Key Baseline Value |
| :--- | :--- | :--- |
| **Snakemake Benchmark** | `results/02_classifiable_regions/_baseline_4x_72t/benchmarks/stream_kraken_cr.tsv` | `40,518.75 s` (`11:15:18`), `mean_load = 115.11%`, `max_rss = 2,678,414 MB` |
| **Kraken 2 Stream Log** | `results/02_classifiable_regions/_baseline_4x_72t/logs/kraken2_stream.log` | `242,682,268` seqs (`97,065.21 Mbp`) in `36,929.33 s` (`394.3 Kseq/m`, `157.70 Mbp/m`) |
| **CR Consumer Log** | `results/02_classifiable_regions/_baseline_4x_72t/logs/consume_kraken_stream_to_cr.log` | `240,101,293` classified windows drained |
| **CR Summary Table** | `results/02_classifiable_regions/_baseline_4x_72t/cr_streaming_summary.tsv` | 44 active species (`cr_bp_2x`, `cr_segments_2x`) |
| **Per-Species CR BEDs** | `results/02_classifiable_regions/_baseline_4x_72t/beds/{t_id}.cr_2x.bed` | All 44 baseline `cr_2x.bed` (and `cr_4x.bed`) files |
| **Stage 5 Reports** | `results/05_reports/_baseline_4x_72t/bins{50,100}_*/` | All 6 completed `bins100_*` and `bins50_*` baseline report directories |
