# Mission: Architect & Implement `ValidateClassificationsSMK` (v2) — High-Throughput Metagenomic eDNA Validation Pipeline

## 1. Role, Authority & Collaborative Workflow
You are an expert Bioinformatics Architect and Snakemake Workflow Engineer working with a research team in the **SweBITS** project (Umeå University). Your mission is to completely refactor and implement **`ValidateClassificationsSMK`**, a production-grade Snakemake pipeline running on the **Kebnekaise HPC cluster (HPC2N, SLURM)** that validates putative taxonomic detections (across all domains of life, with an urgent focus on eukaryotes) from a 10-year airborne environmental DNA (eDNA) shotgun metagenomic timeseries (**516 weekly air-filter samples, ~150–200M paired-end reads per sample**).

### Your Authority & Operating Rules
1. **Full Architectural Authority (With One Sacred Exception)**:
   - You have full authority to restructure the project directory layout, rewrite the `Snakefile`, modularize `workflow/rules/*.smk`, and write clean Python/Bash/Awk scripts in `workflow/scripts/`.
   - **THE SACRED EXCEPTION**: Do  **NOT**  modify or delete files outside this project directory (`/pfs/proj/nobackup/fs/projnb10/snic2019-35-73/private/daniel/02-LJUNGBYHED/01-SERIES_2013/03-ANALYSIS/validations/ValidateClassificationsSMK`).
   -  The core mechanics of my existing  **SLURM profile and SLURM submit scripts** are battle-tested on Kebnekaise after extensive troubleshooting, you may edit them on your own risk ;)
2. **Ask Before Assuming**:
   - Whenever an input file format (e.g., Kraken 2 sample output format, DB manifest column names, FASTQ directory structure), edge-case behavior, or aesthetic/style choice is ambiguous, **stop and ask me clarifying questions** before generating large blocks of dependent code.
   - If you need sample lines from any input file (`manifest.tsv`, `.kraken.out`, sample sheet, existing SLURM profile, etc.), ask me and I will provide them if they are not found.
3. **Code Documentation Philosophy ("The WHY Over The HOW")**:
   - This repository will be hosted on the **SweBITS GitHub organization** for other researchers in our group to clone and run easily.
   - Every Snakemake rule and Python/Bash script must have clear, human-readable comments that explain **WHY** each biological, statistical, or algorithmic choice exists (e.g., why we pad CRs by the step size, why contigs are sorted shortest-to-longest before binning, why `--report` is omitted in Kraken 2). Make complex streaming and mathematical steps intuitive for a biologist or colleague reading the code a year from now.

---

## 2. Scientific & Algorithmic Background

### The Core Hypothesis & Challenge
- **Biological Premise**: If a species is genuinely present in an eDNA sample, shotgun reads classified to that taxon by Kraken 2 should map across its genome in proportion to the length of its **Classifiable Regions (CRs)**—the specific genomic intervals capable of yielding a Kraken 2 classification to that species clade given our exact database (2.6 TB `hash.k2d`) and classification parameters (`--confidence 0.1 --minimum-hit-groups 10`). False-positive detections (conserved sequence cross-mapping, PCR artifacts, or unsequenced sister-species noise) pile into a few localized spikes while leaving most CRs empty.
- **Architectural Inversion (Crucial Design Shift)**:
  - We do **NOT** map eDNA reads to the full reference genome, and we do **NOT** back-map simulated reads with `minimap2`.
  - Instead, we:
    1. Filter reference genomes ($\ge 1,000\text{ bp}$ contigs).
    2. Slice references into deterministic **`400 bp` single-end in silico reads** (`100 bp` step = $4\times$ coverage), stream them through Kraken 2 once, and merge the known genomic coordinates of reads classified to the target species clade (resolved via **`JolTax`**) directly into **Classifiable Regions (CRs)**.
    3. Partition the merged CRs in linear genomic order into **100 bins of linearly increasing target size** ($0.3\times$ to $1.7\times$ mean bin CR length) and extract only those CR sequences as a compact reference (`{taxid}.cr_contigs.fa`).
    4. Extract the real eDNA reads classified to that species clade across all 516 samples and map them **exclusively to `{taxid}.cr_contigs.fa`** (without PCR deduplication).
    5. Run **Recursive Linear Regression with Outlier Pruning** ($R \ge 0.90$) on Mapped eDNA Reads (Y) vs. Bin CR Length (X), reconstruct the cleaned 516-week time series from the validated bins, and generate an overhauled, publication-quality visual report and master summary table.

### Initial Target Taxa (Manuscript Section 2 Benchmark Set)
```tsv
t_id	t_species
325452	Phytophthora kernoviae
40932	Rhopalosiphum padi
43146	Rhopalosiphum maidi
116150	Pararge aegeria
111881	Nymphalis urticae
64459	Pieris rapae
9031	Gallus gallus
30423	Coloeus monedula
8934	Columba palumbus
8932	Columba livia
8032	Salmo trutta
4298	Ilex aquifolium
36056	Boletus edulis
13563	Heterobasidion annosum
746836	Hymenoscyphus fraxineus
27923	Mnemiopsis leidyi
47425	Armillaria borealis
64369	Armillaria calvescens
47426	Armillaria cepistipes
47427	Armillaria gallica
47428	Armillaria ostoyae
47429	Armillaria mellea
```

---

## 3. Detailed Pipeline Architecture & Step-by-Step Specifications

### User Experience & Repository Structure (Snakemake Best Practices)
Design a clean, standard Snakemake directory structure:
- `config/config.yaml`: Algorithmic parameters (with sensible defaults pre-filled) and static cluster/database paths.
- `config/targets.tsv`: Simple user-facing table (`t_id`, `t_species`) of taxa to validate. Users should **not** have to manually enter assembly accessions or FASTA paths into a complex YAML file.
- `workflow/Snakefile`, `workflow/rules/*.smk`, `workflow/scripts/*.py`, `workflow/envs/*.yaml`.
- `profiles/slurm/`: Preserve my existing SLURM profile and submission scripts.
- `results/` & `logs/` & `benchmarks/`: Cleanly organized outputs, per-rule logs, and Snakemake `benchmark:` TSVs.

### Graceful Failure & Ledger Pattern (CRITICAL REQUIREMENT)
If a user inputs 25 target species and 4 of them fail at different pre-checks (e.g., `< 1,000` total eDNA reads across the timeseries, organellar-only genome in the database manifest, or `< 500 kb` of sequence remaining after removing `< 1,000 bp` contigs), **the Snakemake pipeline must NEVER crash or halt the remaining 21 valid species**.
- Implement a clean **Checkpoint / Status Ledger pattern**: every target species from `config/targets.tsv` is tracked through the pipeline and guaranteed to appear in the final master report table with an explicit `validation_status` (e.g., `VALIDATED_CLEAN`, `VALIDATED_PRUNED`, `FAILED_LOW_CORRELATION`, `FAILED_LOW_OCCUPANCY`, `EXCLUDED_LOW_READ_COUNT`, `EXCLUDED_ORGANELLAR_ONLY`, `EXCLUDED_INSUFFICIENT_CONTIG_LENGTH`) and a human-readable `status_reason`.

---

### Stage 0: Clade Resolution (`JolTax`), Pre-Flight Gates, Assembly Selection & eDNA Read Extraction
1. **Taxonomic Clade Resolution with `JolTax`**:
   - Use **`JolTax`** (from the SweBITS GitHub repository) for all taxonomic operations.
   - For each target species `t_id` in `config/targets.tsv`, use `JolTax` to retrieve the full set of descendant TaxIDs in the clade rooted at `t_id` (including the species node itself and all subordinate subspecies, strains, or formae speciales).
2. **Pre-Flight Gate 1 — Minimum Read Count Across Timeseries (`min_total_reads: 1000`)**:
   - Count (or extract) all eDNA read pairs classified by Kraken 2 to any TaxID within the species' `JolTax` clade across the 516 samples.
   - **Statistical Power Threshold**: If a species clade has fewer than `min_total_reads` (default `1,000` read pairs, configurable in `config.yaml`, with an optional `force_include_taxids` override list), mark the species as `EXCLUDED_LOW_READ_COUNT` in the ledger and skip downstream heavy compute. (With 100 bins averaging $0.3\bar{L}$ to $1.7\bar{L}$, 1,000 reads ensures $\lambda \ge 3$ reads in the smallest bin so $\ge 95\%$ of bins are non-zero under Poisson statistics).
   - **eDNA Read Extraction**: Extract the fastp-trimmed paired-end eDNA reads (`R1` and `R2`) assigned to the target species clade across the samples, tagging or preserving the `sample_id` in the read header (`|sample={sample_id}`) so downstream mapping runs in a single pooled pass per species while retaining exact per-sample (weekly) provenance.
3. **Pre-Flight Gate 2 — Organellar-Only Clade Exclusion via DB Manifest**:
   - Inspect the Kraken 2 database manifest file for all assemblies belonging to the species' `JolTax` clade.
   - If the clade contains **only organellar genomes** (mitochondrial, chloroplast, plastid, apicoplast) and no nuclear/chromosomal assembly, exclude the species (`EXCLUDED_ORGANELLAR_ONLY`) and log the reason in the final summary table.
4. **Deterministic Best-Assembly Selection (`select_assemblies.py` Logic)**:
   - Our Kraken 2 database was constructed with at most one assembly per taxonomic node, meaning a species clade may have one assembly at the species node (`t_id`) plus assemblies at subspecies nodes.
   - Among all `a_status == "current"` non-organellar assemblies in the target species' `JolTax` clade within the runtime manifest, select the single reference assembly using this deterministic cascade:
     1. **`a_refseq_category`**: Filter to `'reference genome'` or `'representative genome'` if present.
     2. **`a_level`**: Filter to the highest available assembly level (`Complete Genome` > `Chromosome` > `Scaffold` > `Contig`).
     3. **Species-Node Preference**: If multiple candidates remain tied at the highest RefSeq category and assembly level, prefer the assembly directly assigned to the species node (`t_id`).
     4. **Quality & Date Tie-Breakers**: Highest N50 (`a_contig_n50` for Contig level, else `a_scaffold_n50`) $\rightarrow$ Lowest sequence count (`a_num_contigs` / `a_num_scaffolds`) $\rightarrow$ Most recent `a_release_date` $\rightarrow$ **Deterministic seeded random choice (`random_state=42`)**.
   - Log the selected accession, TaxID, total candidate assemblies in clade, and selection reason (`a_selection_reason`).

---

### Stage 1: Upfront Reference Contig Filtering (`>= 1,000 bp`) & Audit Logging
1. **Filter Short Contigs**:
   - Upon obtaining/locating the selected reference FASTA, filter out all contigs shorter than `min_contig_length: 1000` bp (short contigs in draft assemblies represent unresolved repeats or assembly artifacts).
   - For this first version of the pipeline, treat organellar contigs present inside a full nuclear assembly like any other sequence (do not build special organellar splitters yet).
2. **Mandatory Assembly Filtering Audit Report (`assembly_filtering_summary.tsv`)**:
   - Record for every species: `t_id`, `t_species`, `chosen_accession`, `chosen_assembly_taxid`, `selection_reason`, `raw_contig_count`, `kept_contig_count`, `discarded_contig_pct`, `raw_assembly_bp`, `kept_assembly_bp`, `discarded_bp`, and `discarded_bp_pct`.
   - If `kept_assembly_bp` is below a minimum viable nuclear threshold (e.g., `< 100,000 bp` or fewer than `100` possible windows), mark the species as `EXCLUDED_INSUFFICIENT_CONTIG_LENGTH` and handle gracefully.
   - If  `discarded_contig_pct`  is above a viable threshold (e.g.,  `50%`), mark the species as  `EXCLUDED_FILTERED_ASSEMBLY_PROPORTION`  and handle gracefully.
  
---

### Stage 2: Zero-Disk Streaming `400 bp` Tiling, Single-Pass Kraken 2 & Flanked CR Construction
1. **Hardcoded `400 bp` Single-End In Silico Tiling at $4\times$ Coverage (`100 bp` Step)**:
   - Slice all valid filtered reference genomes using `seqkit sliding -W 400 -s 100 -g` ($4\times$ coverage).
   - Encode the contig name, 0-based/1-based window index (`idx`), exact `start-end` genomic coordinates, and target `t_id` into each FASTA header.
   - **Crucial Optimization ($2\times$ Inside $4\times$)**: Do **NOT** generate separate $2\times$ and $4\times$ FASTA streams for Kraken 2! In a $400\text{ bp}$ window / $100\text{ bp}$ step ($4\times$) tiling, the **even-numbered windows (`idx % 2 == 0`) ARE the exact $2\times$ coverage dataset (`200 bp` step)**! Thus, Kraken 2 only classifies the $4\times$ stream once, and the downstream Python stream consumer simultaneously builds both the **$2\times$ CR BED** and the **$4\times$ CR BED** (and optionally $1\times$ via `idx % 4 == 0`).
2. **Zero-Disk POSIX FIFO Streaming into Kraken 2**:
   - Stream the tiled reads from all target genomes into a POSIX named pipe (`mkfifo "$SNIC_TMP/sim_reads.fa"`) so that 0 GB of simulated reads are written to disk.
   - Run Kraken 2 on a full 72-core `cpu_largemem` node (`--confidence 0.1 --minimum-hit-groups 10`, configurable).
   - **DO NOT generate a Kraken 2 `--report` file** (omit the `--report` flag completely; generating a report on a 2.6 TB DB wastes substantial time at the end of classification).
   - **Modified Kraken 2 Compatibility & C-Speed Stream Filtering**: We use a custom-compiled Kraken 2 fork (`kraken2_SweBITS`) that outputs extra minimizer hit-group telemetry columns alongside the massive $k$-mer hit list string. Pipe `kraken2` stdout directly into `awk -F'\t' '$1 == "C" {print $2 "\t" $3}'` before Python. This strips all unclassified (`U`) reads and drops columns $4\dots N$ at C speed (shrinking stdout volume $>300\times$) while remaining 100% immune to extra columns.
   - **Bulletproof Stream Draining**: Use `set -euo pipefail`, explicitly `wait` on background FIFO producers, flush/fsync all output files in the Python consumer, and verify an explicit completion marker so SLURM never terminates the job before `awk` and Python have completely drained the Kraken 2 stream.
3. **Clade-Aware Coordinate Merging & Step-Size Flanking Buffers (Zero Overlap Guarantee)**:
   - In the Python stream consumer, use **`JolTax`** to check if the classified TaxID (Column 3) belongs to the target species' clade (species node or any descendant node).
   - **Why We Pad CRs by the Step Size ($S$)**: To prevent real eDNA reads that straddle the edge of a diagnostic minimizer cluster from being soft-clipped when mapped exclusively to excised CRs, pad each raw classified window $[s, e)$ by $\pm S$ on both ends (`+200 bp` for the $2\times$ analysis tier where $S = 200\text{ bp}$; `+100 bp` for the $4\times$ tier where $S = 100\text{ bp}$).
   - **Contig-End Clamping & Zero Double-Counting**:
     - Clamp every padded interval strictly to `[0, contig_length]` so no coordinate ever goes negative or exceeds the contig end.
     - Sort and **merge all overlapping or touching padded intervals (`next_start <= curr_end`)** per contig *before* extracting sequences. This guarantees that no base pair of the reference assembly is ever duplicated or double-counted across CRs.
   - Default downstream analysis uses the **$2\times$ CR dataset** (configurable via `coverage_tier: "2x"` or `"4x"` in `config.yaml`).

---

### Stage 3: Small-to-Large Contig Sorting & Linear Variable-Size Binning (`build_cr_bins.py`)
1. **Why & How We Bin (`n_bins: 100`, configurable)**:
   - Correlating at the whole-contig level fails for chromosome-scale assemblies (too few points, extreme macrochromosome leverage) and draft assemblies (thousands of tiny points).
   - Dividing CRs into 100 equal-sized bins (`CR Length = constant`) collapses all X-axis values into a single vertical line, making linear correlation impossible.
   - **Our Algorithm (Easy to Explain in Materials & Methods)**:
     1. Take all contigs that contain at least one merged CR and **sort the contigs from shortest total contig length to longest total contig length** (breaking ties deterministically by contig name).
     2. Compute total classifiable sequence length $L_{\text{CR}}$ across the genome and mean bin length $\bar{L} = L_{\text{CR}} / B$ (where $B = \text{config}[\text{"n\_bins"}]$, default `100`).
     3. Assign each sequential bin $k \in \{1, 2, \dots, B\}$ a linearly increasing target CR length from $0.3\bar{L}$ to $1.7\bar{L}$:
        $$W_k = \bar{L} \times \left(0.3 + 1.4 \times \frac{k - 1}{B - 1}\right)$$
        Because $\frac{0.3 + 1.7}{2} = 1.0$, the sum of target lengths across all $B$ bins equals $L_{\text{CR}}$ ($\sum_{k=1}^B W_k = L_{\text{CR}}$).
     4. Walk linearly along the CR intervals of the sorted contigs (from shortest contig to longest chromosome), filling Bin 1 ($0.3\bar{L}$), then Bin 2, ..., up to Bin $B$ ($1.7\bar{L}$). If a continuous CR interval crosses the capacity boundary of Bin $k$, split that CR interval cleanly at the boundary coordinate so Bin $k$ reaches $W_k$ and the remaining bases start Bin $k+1$.
     5. **Biological & Statistical Benefit**: CRs from small contigs naturally group together into the smaller initial bins, while CRs along large chromosomes stay in contiguous linear neighborhoods inside the larger bins—giving a clean $5.7\times$ horizontal spread on the X-axis with zero random numbers, while ensuring that dropping an outlier bin removes a contiguous genomic neighborhood.
2. **Outputs**:
   - `{taxid}.cr_contigs.fa`: Extracted FASTA of all CR segments, with headers encoding `cr_id`, `bin_id` (`bin_001`..`bin_100`), `source_contig`, and `start-end`.
   - `{taxid}.bins_metadata.tsv`: Tabular metadata per bin (`bin_id`, `target_cr_bp`, `actual_cr_bp`, `num_cr_segments`, `source_contigs_summary`).

---

### Stage 4: Mapping eDNA Reads Exclusively to CR Contigs (No Deduplication)
1. **Index & Align (`minimap2 -ax sr`)**:
   - Index `{taxid}.cr_contigs.fa` with `minimap2` and align the extracted paired-end eDNA reads (`R1`, `R2`) for that species using default `minimap2 -ax sr` settings.
   - **Skip PCR deduplication** (`markdup` is intentionally omitted).
2. **Configurable Read-Length-Proportional SAM Filtering**:
   - Filter alignments in the stream using `samtools view` with parameters exposed in `config.yaml`:
     - Exclude unmapped, secondary, and supplementary alignments (`-F 0x904`).
     - Minimum mapping quality: `min_mapq: 30` (default `30`; configurable).
     - Minimum aligned fraction of read length: `min_align_frac: 0.80` (`(qlen - sclen) >= 0.80 * qlen`).
     - Maximum mismatch fraction: `max_mismatch_frac: 0.01` (`[NM] <= 0.01 * (qlen - sclen)`, i.e., $\ge 99\%$ identity across the aligned portion).
   - Record each passing alignment (`read_id` with its embedded `sample_id`, `bin_id`, `cr_id`, `mapq`, `EditDistance`) in a compact TSV/Parquet file Using **Polars** for fast downstream aggregation.

---

### Stage 5: Recursive Outlier Regression, Time Series Reconstruction & Overhauled Reporting
1. **Pre-Regression Sanity Gate (Minimum Bin Occupancy)**:
   - Aggregate mapped reads per bin ($Y_1 \dots Y_B$) against bin CR length ($X_1 \dots X_B$).
   - Check that at least `min_occupied_bins_frac: 0.50` (default $\ge 50\%$ of bins) have $\ge 1$ mapped read and that the initial or fitted slope is positive ($\beta_1 > 0$). (This prevents a false-positive species with 95 empty bins and 5 artifact spikes from pruning away all its empty bins during bidirectional regression).
2. **Bidirectional Recursive Linear Regression ($R \ge 0.90$)**:
   - Fit an Ordinary Least Squares (OLS) linear regression $Y = \beta_0 + \beta_1 X$ and compute Pearson correlation $R$ (plus Spearman $\rho$ for reporting).
   - If $R \ge \text{target\_r}$ (default `0.90`) and $\beta_1 > 0$: **STOP** $\rightarrow$ Status is `VALIDATED_CLEAN` (if 0 bins removed) or `VALIDATED_PRUNED` (if $\ge 1$ bins removed).
   - While $R < 0.90$ and retained CR length $\ge \text{min\_retained\_cr\_frac} \times L_{\text{CR}}$ (default `0.50`, i.e., $\ge 50\%$ of total CR sequence still retained):
     - Identify the single bin $i^*$ with the **largest absolute vertical residual from the regression line** ($\vert{}Y_i - \hat{Y}_i\vert{}$)—allowing removal of both high-copy repeat/organellar spikes *above* the line and contaminated reference assembly scaffolds *below* the line.
     - Remove bin $i^*$, log its removal step (`iteration`, `bin_id`, `cr_bp_removed`, `reads_removed`, `residual`, `R_before`, `R_after`), and refit the OLS line on the remaining bins.
   - If retained CR length drops below $50\%$ (or slope $\beta_1 \le 0$) before reaching $R \ge 0.90$, mark the species as `FAILED_LOW_CORRELATION` (flagged as unvalidated, **never** silently deleted from the output).
3. **Reconstructing the Cleaned 516-Week Time Series**:
   - Using the `sample_id` tags on the mapped reads, output for every species:
     - `raw_kraken_reads`: Weekly read counts from raw Kraken 2 output.
     - `all_cr_mapped_reads`: Weekly read counts mapping to all 100 CR bins.
     - `validated_clean_bin_reads`: Weekly read counts mapping ** strictly to the retained (non-pruned) CR bins**, together with the species' `validation_status` flag.
4. **Complete Overhaul of Final Visual & Tabular Reports**:
   - Design a modern, publication-quality reporting suite (both machine-readable TSV/Parquet tables via **Polars** and a polished, self-contained **HTML visual dashboard** + high-DPI publication figures via Matplotlib/Seaborn or Plotly):
     1. **Master Validation Summary Table (`validation_master_summary.tsv`)**: Every input taxon from `config/targets.tsv` with full lifecycle metrics: clade read counts, assembly chosen, % assembly removed by the $1,000\text{ bp}$ filter, total CR length ($2\times$ and $4\times$, % of genome), initial $R$ and $\rho$, final $R$ and $\rho$, number of bins pruned, % CR length retained, % mapped reads retained, and final `validation_status` + `status_reason`.
     2. **Per-Species Multi-Panel Diagnostic Figure**:
        - **Panel A (Regression Scatter Plot)**: Bin CR Length (bp) on X vs. Mapped eDNA Reads on Y, clearly distinguishing **Retained Bins** (`●`) from **Pruned Outlier Bins** (`▲`, labeled with pruning step order), showing both the **Initial Regression Line** (dashed grey) and the **Final Validated Regression Line** (solid color with 95% CI band).
        - **Panel B (Pruning Trajectory Curve)**: Dual-trace plot showing how Pearson $R$ increases and % Retained CR Length decreases at each outlier removal step relative to the $R = 0.90$ and $50\%$ retention thresholds.
        - **Panel C (516-Week Time Series Comparison)**: Stacked or overlaid weekly time series comparing Raw Reads vs. Validated Clean-Bin Reads across the 10-year monitoring period.
     3. **Interactive Self-Contained HTML Report (`validation_report.html`)**: Clean, modern typography, executive summary cards, sortable/filterable master status table (including excluded/failed species with clear badges explaining why they failed), and embedded per-species diagnostic plots and resource benchmark summaries.

---

### HPC Resource Management, Temporary Files & Logging
1. **Dynamic SLURM Retries (`attempt`)**:
   - Configure `retries: 5` on all per-species rules (`filter_reference`, `build_cr_bins`, `map_edna_to_cr`, `validate_and_regress`) with dynamic `resources:` lambdas (`mem_mb = lambda wc, attempt: ...`, `runtime = lambda wc, attempt: ...`) so OOM or timeout kills automatically scale up on retry.
   - **Exception**: The single-pass Kraken 2 streaming rule runs on the largest node (`cpu_largemem`, 72 threads, full node RAM) with `retries: 0`.
2. **Scratch & Temporary Files**:
   - Use `$SNIC_TMP` for all FIFOs and intermediate sorting/streaming buffers, and wrap intermediate workflow outputs in Snakemake's `temp()` directive wherever appropriate.
3. **Comprehensive Logging & Benchmarking**:
   - Every rule must define both `log:` and `benchmark:` directives. All Python scripts must use Python's `logging` module with timestamps, clear progress milestones, and explicit error tracebacks.

---

## 4. How We Should Begin
1. First, review these specifications and **ask me any clarifying questions** you need regarding:
   - My existing repository structure, SLURM profile (`profiles/slurm/`), and submit script.
   - The exact schema/columns of my Kraken 2 database manifest file (`database_manifest.tsv`) and how eDNA reads (`input/fastq`) and Kraken outputs (`input/classifications`) are stored across the 516 samples.
   - How `JolTax` is imported and initialized in our environment (available by running import joltax in python, installed via pip (https://github.com/SweBiTS/JolTax)). The path to a joltax-cache you can use is: `/proj/nobackup/snic2019-35-73/private/daniel/98-MISC/joltax_cache/master_db_2025/`
2. Propose the complete **project directory tree** and the user-facing **`config/config.yaml` + `config/targets.tsv`** design for my approval before writing the implementation code.
3. Propose an implementation sequence that prioritizes complex tasks that benefit from using Argon High before simpler tasks that can be solved by Flash 3.8 (at different levels), before you start coding.
4. We'll be working in the git branch `feature/pipeline_overhaul_v2` which track  another branch (`feature/refactor_validation_pipeline`). Goal is to merge this stack into the main branch when we are done.
