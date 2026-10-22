# ValidateClassificationsSMK (v2)

High-throughput Snakemake workflow for validating putative taxonomic detections in the 10-year Ljungbyhed airborne environmental DNA (eDNA) shotgun metagenomic timeseries (**516 weekly samples, 2013–2022**), developed by the **SweBITS** project at Umeå University.

---

## Core Concept: Why & How We Validate Detections

When a species is genuinely present in an environmental sample, shotgun sequencing reads classified to that species clade by Kraken 2 should distribute across its genome in proportion to the length of its **Classifiable Regions (CRs)**—the specific genomic intervals capable of yielding a species-clade classification given our exact database (`2.6 TB hash.k2d`) and parameters (`--confidence 0.1 --minimum-hit-groups 10`).

Conversely, false-positive detections (caused by cross-mapping to conserved rRNA/housekeeping motifs, unsequenced sister species, or reference assembly contaminants) pile into a few localized spikes while leaving most classifiable regions empty.

### Pipeline Overview (6 Automated Stages)

1. **Stage 0 — Clade Resolution (`JolTax`), Pre-Flight Gates & Best-Assembly Selection**:
   - Expands each target TaxID in `config/targets.tsv` into its full descendant clade (`JolTax`).
   - **Pre-Flight Gate 1**: Requires $\ge 1,000$ total clade read pairs across the 516 weeks (`EXCLUDED_LOW_READ_COUNT` if below threshold, ensuring sufficient statistical power across 100 bins).
   - **Pre-Flight Gate 2**: Excludes clades represented only by mitochondrial/plastid genomes (`EXCLUDED_ORGANELLAR_ONLY`).
   - Selects the single best reference assembly in the clade (`RefSeq category` $\rightarrow$ `Assembly level` $\rightarrow$ `Species node` $\rightarrow$ `N50` $\rightarrow$ `Contig count` $\rightarrow$ `Release date` $\rightarrow$ `Seed 42`) and downloads/rehydrates all selected genomes in bulk via NCBI `datasets` with parallel MD5 verification.
2. **Stage 1 — Upfront Contig Filtering ($\ge 1,000\text{ bp}$) & Quality Audit**:
   - Removes short contigs ($< 1,000\text{ bp}$) and audits sequence retention in `results/00_ledger/assembly_filtering_summary.tsv`.
   - Excludes assemblies leaving $< 100\text{ kb}$ (`EXCLUDED_INSUFFICIENT_CONTIG_LENGTH`) or losing $> 50\%$ of base pairs (`EXCLUDED_FILTERED_ASSEMBLY_PROPORTION`).
3. **Stage 2 — Zero-Disk `400 bp` / `200 bp`-Step ($2\times$) Streaming Tiling & Single-Pass Kraken 2 CR Discovery**:
   - Slices all active filtered references into `400 bp` windows stepped every `200 bp` ($2\times$ coverage) using a multi-worker producer (`stream_tiled_fasta.py`) and streams them through a POSIX FIFO into Kraken 2 (`36` OpenMP threads across 2 NUMA sockets, with `intel/2025a` loaded) once.
   - Clade-classified windows are padded by the `200 bp` step size ($\pm 200\text{ bp}$), clamped to contig ends, and merged into disjoint Classifiable Regions (`{taxid}.cr_2x.bed`) with zero base-pair overlap.
4. **Stage 3 — Shortest-to-Longest Contig Sorting & Linear Variable-Size Binning (`100 Bins`)**:
   - Sorts contigs containing CRs from shortest to longest and partitions total CR sequence ($L_{\text{CR}}$) into **100 sequential bins** with linearly increasing capacities from $0.3\bar{L}$ to $1.7\bar{L}$, outputting `{taxid}.cr_contigs.fa`.
5. **Stage 4 — Sample-Aware eDNA Read Extraction & CR-Exclusive Alignment**:
   - Extracts clade reads across all 516 samples using `seqtk subseq`, tagging each header with `|sample=<sample_id>`.
   - Aligns paired-end reads exclusively to `{taxid}.cr_contigs.fa` using `minimap2 -ax sr` (without PCR deduplication) and filters primary alignments by `MAPQ >= 30`, aligned read fraction $\ge 80\%$, and edit distance $\le 1\%$ ($\ge 99\%$ identity).
6. **Stage 5 — Bidirectional Recursive Outlier Regression & 516-Week Signal Reconstruction**:
   - Enforces a minimum bin occupancy sanity gate ($\ge 50\%$ non-empty bins and positive initial slope $\beta_1 > 0$).
   - Iteratively prunes the single bin with the largest absolute vertical residual ($|Y_i - \hat{Y}_i|$) until Pearson $R \ge 0.90$ (`VALIDATED_CLEAN` or `VALIDATED_PRUNED`) or retained CR length drops below $50\%$ (`FAILED_LOW_CORRELATION`).
   - Reconstructs the cleaned 516-week time series (`raw_kraken_reads`, `all_cr_mapped_reads`, `validated_clean_bin_reads`) and builds the self-contained HTML report.

---

## Repository Layout & Input Structure

```text
ValidateClassificationsSMK/
├── config/
│   ├── config.yaml                  # Algorithmic parameters, database paths & SLURM resources
│   ├── targets.tsv                  # User-facing target species table (t_id, t_species, override_accession)
│   └── manifests/                   # Versioned database assembly manifests (*.tsv.gz)
│       └── master_db_2025_manifest.tsv.gz
├── input/                           # External inputs (git-ignored except .gitkeep)
│   ├── classifications/             # 516 weekly per-read Kraken 2 files ({sample}.kraken2.gz)
│   ├── fastq/                       # 516 weekly paired-end FASTQ files ({sample}_R{1,2}.fastq.gz)
│   └── merged_reports/              # Merged weekly canonical counts (lj_stringent_canonical_counts.tsv.gz)
├── slurm/
│   ├── hpc2n/config.yaml            # Battle-tested HPC2N Kebnekaise SLURM profile
│   └── slurm-jobscript.sh           # SLURM job wrapper script
└── workflow/
    ├── Snakefile                    # Main Snakemake entrypoint
    ├── envs/                        # Modular Conda environment definitions (*.yaml)
    ├── rules/                       # Modular Snakemake rule stages (00_*.smk .. 06_*.smk)
    └── scripts/                     # Documented Python scripts & HTML Jinja2 template
```

---

## Quickstart

### 1. Define Target Species (`config/targets.tsv`)
Add the NCBI TaxID (`t_id`), species name (`t_species`), and optional `override_accession` for each taxon you want to validate:
```tsv
t_id	t_species	override_accession
30423	Coloeus monedula	
8934	Columba palumbus	GCA_965250375.1_bColPal1.1_genomic.fna.gz
36056	Boletus edulis	
```
- **Automatic Assembly Selection (default)**: Leave `override_accession` empty and the pipeline will automatically select and download the highest-ranked nuclear assembly for that clade from `database_manifest`.
- **Manual Accession or Local FASTA Override (`override_accession`)**: Optionally specify either an NCBI assembly accession (`GCA_...` / `GCF_...`) or a local genome FASTA filename/path (`*.fa`, `*.fna`, `*.fa.gz`, `*.fna.gz`) placed in `input/local_assemblies/`.
- **Duplicate `t_id` Handling**: If the same `t_id` is listed more than once in `config/targets.tsv`, Stage 0 keeps the first occurrence and automatically skips subsequent duplicates with a `WARNING` in the log so downstream Snakemake wildcards never collide.

### 2. Review Parameters (`config/config.yaml`)
All thresholds (`min_total_reads: 1000`, `tile_window_bp: 400`, `tile_step_bp: 200`, `n_bins: 100`, `target_r: 0.90`, `min_retained_cr_frac: 0.50`) and cluster paths (`config/manifests/master_db_2025_manifest.tsv.gz`, `joltax_cache_dir`, `kraken2_db_path`, `kraken2_toolchain_module: "intel/2025a"`) are pre-configured in [`config/config.yaml`](config/config.yaml).

### 3. Run on Kebnekaise (SLURM)
```bash
# Dry-run to inspect the planned DAG
snakemake -n

# Launch on HPC2N Kebnekaise using the SLURM profile and Conda environments
snakemake --profile slurm/hpc2n --use-conda

# Or launch a run with alternative binning or mapping parameters (reuses Stages 0–2 automatically
# and writes Stage 3–5 outputs to parameter-tagged subdirectories):
snakemake --profile slurm/hpc2n --use-conda --config n_bins=50 min_mapq=20 min_align_frac=0.80 max_mismatch_frac=0.05
```

---

## Output Directory Guide (Parameter-Namespaced Stages 3–5)

Upstream Stages 0–2 (`results/00_ledger/`, `results/01_inputs/`, and unbinned BED files in `results/02_classifiable_regions/{t_id}/`) are shared across all runs. Stage 3 is namespaced by **`{bin_tag}` = `bins{n_bins}_s{min_scale}-{max_scale}`** (default: **`bins100_s0p3-1p7`**), and Stages 4–5 (`03_mapping`, `04_validation`, `05_reports`, plus `logs/` and `benchmarks/`) are namespaced by **`{run_tag}` = `{bin_tag}_{map_tag}`** (e.g., **`bins100_s0p3-1p7_mapq30_aln80_id99`**), so running new binning or mapping thresholds never overwrites previous results:

| Path | Description |
| :--- | :--- |
| `results/05_reports/{run_tag}/validation_report_10panel_{run_tag}.html` | **10-Panel Widescreen Self-Contained HTML Dashboard** embedding the 10-panel diagnostic suite (`{t_id}.diagnostic_10panel.png`) for every active species alongside executive KPIs, filterable ledger tables, and SLURM benchmarks. |
| `results/05_reports/{run_tag}/validation_report_4panel_{run_tag}.html` | **4-Panel (2×2) Self-Contained HTML Dashboard** (also written to `validation_report_{run_tag}.html`) embedding the classic 4-panel `2×2` diagnostic plots (`{t_id}.diagnostic_4panel.png`). |
| `results/05_reports/{run_tag}/validation_master_summary_{run_tag}.tsv` | **Master Lifecycle Table** (also saved as `.parquet`) with all Stage 0–5 metrics and final `validation_status` for every target taxon. |
| `results/05_reports/{run_tag}/validated_timeseries_long_{run_tag}.parquet` | **Unified 516-Week Time Series** (`raw_kraken_reads`, `all_cr_mapped_reads`, `validated_clean_bin_reads`) across all target taxa. |
| `results/04_validation/{run_tag}/{t_id}_{Species_name}/` | **Per-Species Results Folder** containing `{t_id}.diagnostic_10panel.{png,pdf}`, `{t_id}.diagnostic_4panel.{png,pdf}`, `{t_id}.bin_regression_data.tsv`, `{t_id}.enhanced_bin_metrics.tsv`, `{t_id}.pruning_trajectory.tsv`, and `{t_id}.weekly_timeseries_516w.tsv`. |
| `results/03_mapping/{run_tag}/{t_id}/` | **Alignment Parquet & Identity Audit** containing `{t_id}.passing_alignments.parquet`, `{t_id}.mapping_stats.tsv`, and `{t_id}.identity_histogram_lt99.tsv` (exact 516-week `<90%` through `100%` 1%-bin read-pair counts). |
| `results/02_classifiable_regions/{bin_tag}/{t_id}/` | **Binned Classifiable Region Reference & Metadata** containing `{t_id}.cr_contigs.fa` and `{t_id}.bins_metadata.tsv`. |
| `results/00_ledger/assembly_filtering_summary.tsv` | **Stage 1 Contig Filtering Audit Table** (shared across runs) recording raw vs. kept contigs and base pairs ($\ge 1,000\text{ bp}$). |

---

## 10-Panel Diagnostic Dashboard Guide (`{t_id}.diagnostic_10panel.png`)

> **Full Plain-English & Mathematical Guide**: See [**`DASHBOARD_README.md`**](DASHBOARD_README.md) for a complete panel-by-panel reference containing both a plain-English guide (for biologists and manuscript reviewers) and the full mathematical specification (formulas, normalizations, statistical cutoffs, and methodological rationale).

The widescreen figure (`32.0 × 18.2 in`, 300 DPI) uses a unified `3 × 4` column-aligned grid ordered alphabetically row-by-row from left to right (`A, B, C, D` in Row 0; `E, F, G` in Row 1; `H, I` in Row 2; and `J` spanning Rows 1–2 in Col 3) with italic species titles and report-matched colored status badges:

- **Panel A — Genomic Proportionality Regression (+ Top Structural Strip & Horizontal Cosine Colorbar)**:
  - Regresses mapped eDNA read pairs ($Y$) against Classifiable Region length ($X$, sorted shortest $\rightarrow$ longest contig across 100 bins), with vertical whiskers showing the intra-bin length-weighted coverage interquartile range (`p25–p75`).
  - A slim top marginal strip (`14.5%` height, locked to the exact same X-axis limits) displays the number of distinct contigs per bin (gray bars) and median CR segment length in kb (blue line).
  - Points are colored by their **516-week Temporal Cosine Similarity** (`Viridis`, `0.0–1.0`), with a shared **horizontal colorbar placed directly between Panel A and Panel E**.
- **Panel B — Outlier Pruning Trajectory**:
  - Tracks linear correlation (Pearson $R$, blue) and Retained CR sequence % (orange, capped at `100%`) across iterative outlier bin removal steps against the validation thresholds ($R \ge 0.90$, Retained CR $\ge 50\%$).
- **Panel C — Outlier Classification Quadrant (`Dot Size = Gini Inequality`)**:
  - Plots **Excess Read Depth** ($\log_2[\text{Observed} / \text{Expected}]$, where $\text{Expected}$ is the bin's length-proportional share of clean CR reads $E_b^{\text{prop}} = \bar{d}_{\text{clean}} \cdot L_b$) on X against **Temporal Concordance (Cosine Sim)** on Y to separate **Endogenous Multi-Copy Repeats** (high depth, concordant phenology) from **Foreign Cross-Mapping Artifacts** (high depth, discordant time series).
  - Marker size encodes the **Intra-Bin Length-Weighted Gini Coefficient ($G_b \in [0, 1]$)**, measuring how unevenly reads are concentrated across CRs within each bin (`0.20 = Uniform CRs`, `0.55 = Moderate`, `≥0.85 = CR Hotspot`).
- **Panel D1, D2 & D3 — eDNA Read Attrition, Alignment Identity Spectrum & Reference Sequence Funnel**:
  - **D1 (`eDNA Read Attrition`)**: Centered symmetrical trapezoidal funnel with strictly proportional bar widths tracking read pairs across 5 filter stages (`1. Raw Kraken Clade Pairs` $\rightarrow$ `2. CR Mapped MAPQ ≥ 30` $\rightarrow$ `3. Passed ≥ 80% Aligned` $\rightarrow$ `4. Passed ≥ 99% Identity` $\rightarrow$ `5. Clean Validated Bins`).
  - **D2 (`Alignment Identity Spectrum`)**: 1%-bin sequence identity distribution (`<90%` through `100%`) for CR-aligned read pairs passing $\ge 80\%$ read length aligned (`<99%` rejected in crimson vs. `99%` and `100%` passing in emerald/amber).
  - **D3 (`Reference Sequence Funnel`, stacked under D2)**: Centered symmetrical trapezoidal funnel with strictly proportional bar widths tracking reference sequence retention (`1. Total Genome Assembly` $\rightarrow$ `2. Classifiable Regions (CR)` $\rightarrow$ `3. Clean CRs (Post-Pruning)`).
- **Panel E — Genome-Wide CR Coverage Map (Manhattan Plot)**:
  - Plots read density (`Reads / kb`) for every individual CR segment across cumulative CR coordinates (`bin_001` $\rightarrow$ `bin_100`) on a `symlog(linthresh=1.0)` y-axis (shaded slate band for the linear `0 ≤ y ≤ 1` regime; $\log_{10}$ above `y = 1`), with rose vertical bands highlighting pruned bins.
- **Panel F — Weekly Signal Retention ($\log_{10}(\text{Reads}+1)$)**:
  - Compares weekly raw Kraken 2 clade read pairs against post-validation clean-bin read pairs across all 516 weeks, reporting linear Pearson $R$, log-space $R_{\log}$, and Spearman $\rho$.
- **Panel G — Bin Temporal Occupancy (`%` of 516 Weeks Occupied)**:
  - Plots each bin's **Weekly Sample Occupancy (% of 516 weeks with $\ge 1$ read pair)** against its **Temporal Concordance (Cosine Sim)**, paired with a top marginal occupancy histogram, Clean Median Occupancy (`%` and `weeks`), and Effective Active Weeks ($N_{\text{eff}} = 1/\sum p_t^2$).
- **Panel H — 10-Year Longitudinal eDNA Time Series (2013–2022, 516 Weeks)**:
  - Stacked directly below the Manhattan plot (`Panel E`) at identical column width, comparing weekly `Raw Kraken Clade Pairs` (gray), `All-CR Mapped Pairs` (blue), and `Validated Clean-Bin Pairs` (status color).
- **Panel I — 10-Year Seasonal Phenology Profile (Weeks 1–52, 3w Circular Smooth)**:
  - Spans the combined width of `Panels B + C` (`F + G`), comparing the 10-year pre-validation vs. post-validation seasonal cycle smoothed with a `[0.25, 0.50, 0.25]` circular moving average across the Dec–Jan boundary, with dual Y-axes (`Share of Annual Signal %/week` on left; `Mean Clean Pairs / Week` on right).
- **Panel J — Bin × Week Hierarchical Clustermap**:
  - Clusters all 100 bins by their $\log_{10}(\text{Reads/kb} + 1)$ temporal profiles across active weeks using Ward's minimum-variance linkage (`method="ward"`, `metric="euclidean"`, `optimal_ordering=True`), with bottom-labeled vertical strips (`CosSim`, `Occ%`, `Status`), a left-aligned panel title matching all other panels, and 2 slim horizontal colorbars (`Occupancy % of 516 Wks` and `Read Density`).

### Validation Status Codes
- `VALIDATED_CLEAN`: Passed occupancy gate ($\ge 50\%$ bins) and achieved Pearson $R \ge 0.90$ ($\beta_1 > 0$) across all 100 bins with **0 bins pruned**.
- `VALIDATED_PRUNED`: Passed occupancy gate and achieved Pearson $R \ge 0.90$ ($\beta_1 > 0$) after bidirectional outlier pruning while retaining $\ge 50\%$ of total CR sequence.
- `FAILED_LOW_OCCUPANCY`: Fewer than $50\%$ of CR bins contained $\ge 1$ mapped read pair, or initial OLS slope was $\le 0$.
- `FAILED_LOW_CORRELATION`: Could not reach $R \ge 0.90$ before hitting the $50\%$ CR sequence retention floor.
- `EXCLUDED_LOW_READ_COUNT`: Fewer than `min_total_reads` ($1,000$) clade read pairs across the 516-week series.
- `EXCLUDED_ORGANELLAR_ONLY`: Clade represented solely by mitochondrial/chloroplast/plastid sequences in the Kraken 2 database.
- `EXCLUDED_INSUFFICIENT_CONTIG_LENGTH`: Fewer than $100\text{ kb}$ of sequence remained after removing $< 1,000\text{ bp}$ contigs.
- `EXCLUDED_FILTERED_ASSEMBLY_PROPORTION`: More than $50\%$ of assembly base pairs were discarded by the $< 1,000\text{ bp}$ filter.

---

## Technical Documentation & In-Depth Analyses

- [`PIPELINE_OVERVIEW.md`](PIPELINE_OVERVIEW.md): **Pipeline Overview & Methodology Guide (Plain-English + Technical Walkthrough)**.
  A colleague-friendly overview of what the pipeline does, why each of the 6 stages (`Stage 0`–`Stage 5`) is needed, and how the parameters work under the hood.
- [`DASHBOARD_README.md`](DASHBOARD_README.md): **10-Panel Diagnostic Figure — Reader's Guide**.
  A page-by-page guide to every panel (`A` through `J`) and validation status code, with a clickable panel index, plain-English explanations, and concise technical notes.
- [`PIPELINE_VIDEO_PROMPT.md`](PIPELINE_VIDEO_PROMPT.md): **Instructional Video Storyboard & AI Prompt Pack**.
  10-scene schematic storyboard, voiceover script, and copy-paste prompts for generating an animated explainer video of the pipeline.
- [`docs/NCBI_ASSEMBLY_SUPPRESSION_ANALYSIS.md`](docs/NCBI_ASSEMBLY_SUPPRESSION_ANALYSIS.md): **Empirical Analysis of NCBI Assembly Suppressions & Grand-Scale Workarounds**.
  An empirical audit of all 107,575 GenBank and RefSeq assemblies in `input/database_manifest.tsv`, explaining why re-annotation and submitter requests cause suppressions, evaluating underlying GenBank availability, and outlining automated workarounds (such as `GCF_` $\rightarrow$ `GCA_` fallback and handling DToL v2 chromosome updates) for future grand-scale executions.

