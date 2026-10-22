# ValidateClassificationsSMK — Pipeline Overview & Methodology Guide

*A plain-English and technical guide to how we validate species detections in the 10-year SweBITS Ljungbyhed airborne eDNA time series (516 weekly air-filter samples, 2013–2022).*

---

## 1. The Big Picture (In Plain English)

### What Problem Does This Pipeline Solve?
Across 10 years of weekly air-filter sampling at Ljungbyhed (516 samples), we sequence hundreds of millions of environmental DNA (eDNA) fragments and classify them using **Kraken 2** against a massive reference database. Starting from a list of candidate species of interest (`config/targets.tsv`), we want to know which of those detections are genuinely present in the air and which are classification artifacts.

While Kraken 2 is fast and sensitive, **not every species detection is real**:
* **What a True Positive looks like**: When an organism (a bird, insect, plant, or fungus) is genuinely present in the environment, its cells break down randomly. As a result, the DNA fragments caught on our air filters come from **all across its genome**.
* **What a False Positive looks like**: Many organisms share highly conserved DNA sequences (such as ribosomal RNA genes, mobile elements, or low-complexity repeats), and some reference genomes contain bacterial or human contamination. Reads from a common microbe in the air can match one tiny conserved spot in a rare species' genome, creating a **massive stack of reads in a single location** while the rest of that species' genome is completely empty.

### How Does the Pipeline Separate True Species from False Positives?
Instead of trusting raw read counts, this pipeline tests a simple physical principle for each candidate species on our target list:
> **If a species is truly present, the number of eDNA reads mapping to a region of its genome should be directly proportional to the size of that region.**

However, we cannot simply check whether reads cover the *entire* reference genome, because parts of every genome are conserved across many species and **cannot be classified to the species level by Kraken 2 in the first place**.

So the pipeline works in two main phases:
1. **Map the "Target Zones" (Stages 0–3)**: For each target species, it selects the reference genome that was included in our Kraken 2 database, shreds it into simulated reads, runs them through Kraken 2, and maps the exact intervals of the genome that Kraken 2 is actually capable of assigning to that species—called its **Classifiable Regions (CRs)**. It then divides those Classifiable Regions into 100 bins of gradually increasing size.
2. **Test Proportionality & Clean Out Spikes (Stages 4–5)**: It maps our real air-filter eDNA reads onto those 100 bins and checks whether larger bins collect proportionally more reads (a straight diagonal line). If a few bins have abnormally high read spikes (from cross-reactive repeats), the pipeline iteratively prunes those noisy bins one by one. If the remaining genome shows a tight linear relationship ($R \ge 0.90$), the species is **validated** and we reconstruct its clean 10-year seasonal time series.

---

## 2. High-Level Pipeline Diagram

```text
  [Target Species List]                      [516 Weekly Air-Filter Samples]
   (config/targets.tsv)                        (Kraken 2 + Paired FASTQs)
            │                                               │
            ▼                                               │
┌─────────────────────────────────────────┐                 │
│ STAGE 0: Pre-Flight & Reference Choice  │                 │
│ • Check >= 1,000 total reads            │                 │
│ • Pick best genome from Kraken DB       │                 │
└───────────────────┬─────────────────────┘                 │
                    ▼                                       │
┌─────────────────────────────────────────┐                 │
│ STAGE 1: Reference Contig Filtering     │                 │
│ • Remove tiny contigs (< 1,000 bp)      │                 │
│ • Check >= 100 kb & >= 50% genome kept  │                 │
└───────────────────┬─────────────────────┘                 │
                    ▼                                       │
┌─────────────────────────────────────────┐                 │
│ STAGE 2: Find Classifiable Regions (CR) │                 │
│ • Shred reference into 400 bp windows   │                 │
│ • Classify with Kraken 2 & merge hits   │                 │
└───────────────────┬─────────────────────┘                 │
                    ▼                                       │
┌─────────────────────────────────────────┐                 │
│ STAGE 3: Variable-Size Linear Binning   │                 │
│ • Sort contigs shortest -> longest      │                 │
│ • Cut CRs into 100 bins of growing size │                 │
└───────────────────┬─────────────────────┘                 │
                    │                                       │
                    └───────────────┐       ┌───────────────┘
                                    ▼       ▼
                    ┌───────────────────────────────────────┐
                    │ STAGE 4: Extract & Map eDNA Reads     │
                    │ • Pull species reads from 516 samples │
                    │ • Stringently align to the 100 CR bins│
                    └───────────────────┬───────────────────┘
                                        ▼
                    ┌───────────────────────────────────────┐
                    │ STAGE 5: Prune Outliers & Validate    │
                    │ • Fit Read Count (Y) vs Bin Size (X)  │
                    │ • Prune repeat spikes until R >= 0.90 │
                    │ • Output clean 10-year time series    │
                    └───────────────────────────────────────┘
```

---

## 3. Step-by-Step Guide (Plain English, Why & Technical Details)

### Stage 0 — Pre-Flight Checks & Reference Genome Selection

* **In Plain English (What happens)**:
  Starting from our list of target species (`config/targets.tsv`), the pipeline first checks each species to make sure it has enough sequencing reads across the 10-year dataset to be statistically testable, and that a nuclear genome (not just a tiny mitochondrial or chloroplast sequence) is present in our Kraken 2 database. If multiple assemblies for that species clade were included in the Kraken 2 database, it automatically selects the highest-quality one and downloads its sequence from NCBI.
* **Why We Do This**:
  1. *Matching the Kraken 2 database*: We must select our reference genome strictly from the manifest of assemblies that were built into the Kraken 2 database, because those are the exact sequences Kraken 2 used to classify our air-filter reads.
  2. *Statistical power*: If a species has fewer than 1,000 reads across all 10 years, dividing those reads into 100 genome bins leaves almost every bin empty—making linear regression unreliable.
  3. *No organellar-only traps*: Mitochondria and chloroplasts are tiny (16–150 kb) and exist in hundreds of copies per cell; they cannot represent whole-genome coverage.
* **Technical Details**:
  - **Clade expansion**: Expands each target species TaxID (`t_id`) in `config/targets.tsv` into its full descendant clade (the species plus any subspecies or strains under it).
  - **Pre-Flight Gate 1 (`EXCLUDED_LOW_READ_COUNT`)**: Queries the merged weekly count table across all 516 samples and requires $\ge 1,000$ total clade read pairs (`min_total_reads: 1000`).
  - **Pre-Flight Gate 2 (`EXCLUDED_ORGANELLAR_ONLY`)**: Checks the Kraken 2 database manifest (`master_db_2025_manifest.tsv.gz`) to verify that at least one nuclear genome assembly exists in the clade.
  - **Deterministic 4-tier assembly selection**: Among the assemblies present in the Kraken 2 database manifest for that clade, ranks candidates by (1) RefSeq category (`reference genome` > `representative genome`), (2) Assembly level (`Complete Genome` > `Chromosome` > `Scaffold` > `Contig`), (3) Exact species TaxID match over sub-taxids, and (4) Contiguity (`contig N50`, fewest contigs, newest release date).
  - **Bulk download**: Fetches the selected assemblies via NCBI `datasets` (`--dehydrated` + `rehydrate`) with parallel MD5 checksum verification.

---

### Stage 1 — Reference Contig Filtering & Assembly Quality Audit

* **In Plain English (What happens)**:
  Draft genome assemblies often contain thousands of tiny, unplaced DNA scraps alongside their main chromosomes. Stage 1 strips out all genome fragments shorter than 1,000 base pairs and checks that at least half of the genome (and at least 100,000 base pairs in total) survived the cleanup.
* **Why We Do This**:
  Short contigs ($< 1\text{ kb}$) in draft eukaryotic genomes are disproportionately made up of uncollapsed high-copy repeats, sequencing adapter artifacts, or foreign bacterial contamination accidentally sequenced with the host. Furthermore, our Stage 2 sliding window is `400 bp` long, so tiny contigs suffer severe edge-truncation bias. If throwing away $< 1\text{ kb}$ scraps removes more than half of the assembly, the reference genome is too fragmented to serve as a reliable reference.
* **Technical Details**:
  - **Length filter**: Retains only contigs $\ge 1,000\text{ bp}$ (`min_contig_length: 1000`) and writes `{t_id}.{acc}.filtered.fna` plus `{t_id}.{acc}.contig_lengths.tsv`.
  - **Stage 1 Gate 1 (`EXCLUDED_INSUFFICIENT_CONTIG_LENGTH`)**: Requires `kept_assembly_bp >= 100,000` (`min_kept_assembly_bp: 100000`), ensuring at least ~1 kb of sequence per bin across 100 bins.
  - **Stage 1 Gate 2 (`EXCLUDED_FILTERED_ASSEMBLY_PROPORTION`)**: Excludes assemblies where removing $< 1\text{ kb}$ contigs discards $> 50\%$ of the total base pairs (`max_discarded_bp_pct: 50.0`).
  - **Snakemake Checkpoint**: Emits `results/00_ledger/active_species.tsv` so only species passing both Stage 0 and Stage 1 trigger heavy downstream compute.

---

### Stage 2 — In-Silico Genome Shearing & Classifiable Region (CR) Discovery

* **In Plain English (What happens)**:
  We take the cleaned reference genome and shred it on the computer into millions of overlapping 400-base-pair pieces ("simulated reads"), stepping forward 200 base pairs at a time. We run all of those simulated pieces through Kraken 2 using the exact same database and settings we use for our air-filter samples. Overlapping pieces that Kraken 2 successfully classifies as our target species are merged into contiguous blocks, and we then pad the left and right ends of each block by 200 base pairs to form the species' **Classifiable Regions (CRs)**.
* **Why We Do This**:
  1. *Defining the fair denominator*: Large parts of any genome (conserved genes, low-complexity regions, or sequences shared with close relatives) are assigned by Kraken 2 to higher taxonomic ranks (like genus, family, or Kingdom) or left unclassified. Since our air-filter reads for a species were originally selected *because* Kraken 2 assigned them to that species clade, those reads can **only** come from the parts of the genome that Kraken 2 recognizes as species-specific.
  2. *Why we pad the ends of each CR by 200 bp (the step size)*: Because we step the 400 bp window forward by 200 bp at a time, the true boundary where sequence becomes species-specific could start up to 200 bp before the first classified window or end up to 200 bp after the last one. In addition, a real paired-end sequencing fragment in our air filter might overlap the edge of a classifiable region enough to get classified by Kraken 2 while its mate extends into the adjacent flanking sequence. Padding the outer boundaries by $\pm 200\text{ bp}$ (matching the step size) ensures those legitimate boundary-spanning read pairs can align cleanly without falling off the edge in Stage 4.
* **Technical Details**:
  - **Sliding-window tiling & classification**: Slices all active filtered genomes into `400 bp` windows stepped every `200 bp` ($2\times$ tiling coverage) and streams them through `kraken2` (`--confidence 0.1 --minimum-hit-groups 10`), loading the 2.6 TB database into RAM only once.
  - **Interval merging & flank extension**: Contiguous/overlapping windows classified to the target species clade are merged and extended on both outer flanks by the `200 bp` step size (`[start - 200, end + 200)`, clamped to `[0, contig_length]` and coalesced if adjacent padded flanks touch), producing `{t_id}.cr_2x.bed`.

---

### Stage 3 — Shortest-to-Longest Contig Sorting & Variable-Size Linear Binning

* **In Plain English (What happens)**:
  Now that we have the species' Classifiable Regions (CRs), we divide them into **100 genomic bins** so we can plot **Bin CR Size (X-axis)** against **Mapped Read Count (Y-axis)**.
  Instead of making all 100 bins the exact same size, we make them grow steadily from smallest (`bin_001`) to largest (`bin_100`):
  - Suppose we add up all of a species' Classifiable Regions into one total CR length ($L_{\text{CR}}$). If we split that CR sequence into **100 equal bins**, each bin would hold **1.0% of the total CR sequence** (that 1.0% is the "average bin size").
  - Instead, we scale the 100 bins linearly from **0.3× the average bin size (0.3% of total CR sequence)** in `bin_001` up to **1.7× the average bin size (1.7% of total CR sequence)** in `bin_100`.
  - **Concrete example**: If a species has **100 Mb** of total Classifiable Regions (so the average bin size is **1.0 Mb**):
    * `bin_001` is assigned **0.30 Mb** (300 kb, or **0.3%** of total CRs)
    * `bin_002` is assigned **0.314 Mb** (314 kb)
    * ...
    * `bin_050` is assigned **~1.00 Mb** (1,000 kb, or **1.0%** of total CRs)
    * ...
    * `bin_100` is assigned **1.70 Mb** (1,700 kb, or **1.7%** of total CRs)
    * Because $(0.3\% + 1.7\%) / 2 = 1.0\%$, all 100 bins still add up to **100% of the species' Classifiable Regions (100 Mb)** to the exact single base pair!
* **Why We Do This**:
  1. *Why equal-sized bins fail*: If all 100 bins held 1.0% of the CR sequence (e.g., 1.0 Mb each), every single point on our scatter plot would sit at the exact same X-axis value ($X = 1.0\text{ Mb}$). They would form a vertical column of dots with zero horizontal spread, making it mathematically impossible to fit a regression line.
  2. *Why whole-contig binning fails*: Using whole chromosomes/contigs as bins fails at both extremes—a chromosome-level bird genome might only have 30 chromosomes (where one giant chromosome dominates the regression), while a draft moss genome might have 50,000 tiny contigs full of zero-read noise.
  3. *Why sorting contigs shortest-to-longest + linear bin scaling works*: Making `bin_100` ($1.7\%$) nearly six times wider than `bin_001` ($0.3\%$) guarantees a **$5.67\times$ horizontal spread** on the X-axis for *every* species with zero random numbers. Sorting contigs shortest-to-longest before filling the bins naturally groups small draft contigs into the smaller initial bins while keeping large chromosomal CRs together in the larger bins.
* **Technical Details**:
  - **Target bin capacities**: Given total CR length $L_{\text{CR}}$ and $B = 100$ (`n_bins: 100`), mean bin length is $\bar{L} = L_{\text{CR}} / B$. Bin $k \in \{1 \dots B\}$ receives target capacity:
    $$W_k = \bar{L} \cdot \left(0.3 + 1.4 \times \frac{k - 1}{B - 1}\right)$$
    Cumulative integer rounding guarantees $\sum_{k=1}^B W_k = L_{\text{CR}}$ to the exact single base pair.
  - **Exact boundary splitting**: [`build_cr_bins.py`](workflow/scripts/build_cr_bins.py) walks linearly along the sorted CR intervals, splitting an interval cleanly whenever a bin capacity boundary is reached, and writes the excised CR sequences to `results/02_classifiable_regions/{BIN_TAG}/{t_id}/{t_id}.cr_contigs.fa` with bin metadata encoded in each header (`>cr_000001|bin=bin_001|contig=...|coords=s-e`).

---

### Stage 4 — Sample-Aware eDNA Read Extraction & Stringent CR Mapping

* **In Plain English (What happens)**:
  The pipeline scans all 516 weekly air-filter samples, pulls out every paired-end sequencing read that Kraken 2 originally assigned to our target species clade, and labels each read with the week it came from. It then maps those reads directly against the species' excised Classifiable Regions using a full sequence aligner (`minimap2`) and discards any read pair that aligns with low confidence, only matches partially, or has too many mismatches.
* **Why We Do This**:
  1. *Why Kraken 2 alone isn't enough*: Kraken 2 only checks short 35-bp k-mers (minimizers); a read can trigger a Kraken 2 match off a couple of short k-mer hits even if the rest of the 150-bp read diverges substantially from the reference genome. Full sequence alignment verifies the entire read pair across both mates.
  2. *Why we map ONLY to CR contigs*: Since the reads were pre-selected by Kraken 2, mapping them to a compact FASTA containing only the Classifiable Regions is vastly faster and prevents reads from scattering into unclassifiable repeat regions.
* **Technical Details**:
  - **Two-stage extraction (`02_extract_edna_reads.smk`)**: Streams all 516 weekly Kraken 2 files + FASTQs once via `seqtk subseq`, tagging headers with `|sample=Lj_YYYY_WW_...` and pooling per-species paired FASTQs (`results/01_inputs/fastq/taxID-{t_id}_R{1,2}.fq.gz`).
  - **Alignment & proportional SAM filtering (`05_map_edna_to_cr.smk`, `filter_sam_to_parquet.py`)**: Aligns with `minimap2 -ax sr` and requires proper primary pairs (`-f 0x2 -F 0x90C`, `RNEXT == "="`) where **both R1 and R2** map within the **same CR** and individually pass all three criteria:
    * **Mapping Quality**: $\text{MAPQ} \ge 30$ (`min_mapq: 30` by default).
    * **Aligned Read Fraction**: $(q_{\text{len}} - s_{\text{clip}}) \ge 0.80 \times q_{\text{len}}$ (`min_align_frac: 0.80`, i.e., $\ge 80\%$ of the read must align without soft/hard clipping).
    * **Sequence Identity**: Edit distance $\text{NM} \le \text{max\_mismatch\_frac} \times (q_{\text{len}} - s_{\text{clip}})$ (`0.01` = $\ge 99\%$ identity by default).
  - **Output**: Writes passing read pairs directly to a compressed Polars Parquet table (`{t_id}.passing_alignments.parquet`).

---

### Stage 5 — Iterative Outlier Pruning, Regression Validation & Time-Series Reconstruction

* **In Plain English (What happens)**:
  Finally, for each species, we plot the 100 bins on a graph: **Bin CR Length on the X-axis** versus **Mapped Read Pairs on the Y-axis**.
  - First, we check a basic sanity gate: at least half of the 100 bins must have reads in them, and larger bins must generally have more reads than smaller bins (positive slope).
  - If the points already form a tight straight line (correlation $R \ge 0.90$), the species immediately passes as **Validated Clean**.
  - If a few bins have huge read spikes (for example, a conserved rRNA repeat or mobile element that slipped into the CRs) that pull the line off-course, the pipeline removes the single worst outlier bin and re-checks the line. It repeats this one bin at a time until either the line becomes straight ($R \ge 0.90 \rightarrow$ **Validated Pruned**) or we have thrown away half of the species' Classifiable Regions without finding a linear trend ($\rightarrow$ **Failed**).
  - Using **only the clean, retained bins**, the pipeline rebuilds the species' 10-year (516-week) abundance time series.
* **Why We Do This**:
  Even a genuinely present species often has 3 to 10 genomic bins containing multi-copy repeats or conserved motifs that attract extra cross-mapping reads from other airborne microbes. Without outlier pruning, a single 10,000-read repeat spike in one bin would ruin the correlation ($R \approx 0.15$) of an otherwise textbook true positive, AND would inject fake background spikes into its 10-year time series. Conversely, the **50% occupancy gate** and **50% retained CR safety floor** prevent a false positive (which *only* has a few spikes and no real genome-wide coverage) from passing validation.
* **Technical Details**:
  - **Pre-regression occupancy gate (`FAILED_LOW_OCCUPANCY`)**: Requires $\ge 50\%$ of bins (`min_occupied_bins_frac: 0.50`) to contain $\ge 1$ mapped read pair, and initial OLS slope $\beta_1 > 0$.
  - **Bidirectional recursive OLS pruning (`validate_and_regress.py`)**: At each step, fits Ordinary Least Squares ($Y = \beta_0 + \beta_1 X$) across currently retained bins. If Pearson $R < 0.90$ (`target_r: 0.90`), removes the single bin $i^*$ with the largest absolute vertical residual $|Y_i - \hat{Y}_i|$, provided removing it keeps retained CR sequence $\ge 50\%$ of $L_{\text{CR}}$ (`min_retained_cr_frac: 0.50`).
  - **516-week signal reconstruction**: Sums weekly read pairs exclusively across retained clean bins (`validated_clean_bin_reads`), stripping out cross-reactive temporal spikes.

---

## 4. What the Pipeline Produces (Verdicts & Reports)

Every run is automatically saved in its own parameter-tagged folder (e.g., `bins100_s0p3-1p7_mapq30_aln80_id99`) so previous results are never overwritten when testing different binning or mapping settings.

### The 4 Validation Verdict Categories

| Verdict Badge | Color | Meaning |
| :--- | :--- | :--- |
| **`VALIDATED_CLEAN`** | **Deep Green** | **Textbook true positive**: Reached Pearson $R \ge 0.90$ across all 100 bins right out of the box with **0 bins pruned**. |
| **`VALIDATED_PRUNED`** | **Emerald Green** | **Confirmed true positive after repeat cleanup**: Reached $R \ge 0.90$ and positive slope after excising 1 or more noisy outlier bins while retaining $\ge 50\%$ of the Classifiable Region sequence. |
| **`FAILED_LOW_OCCUPANCY`**<br>**`FAILED_LOW_CORRELATION`** | **Red** | **Rejected false positive / unvalidated**: Either $>50\%$ of the genome's CR bins were completely empty (classic isolated repeat spike), or pruning all the way down to the $50\%$ CR floor still failed to uncover a linear relationship ($R < 0.90$). |
| **`EXCLUDED_*`** | **Gray** | **Pre-flight / reference gate exclusion**: Excluded before mapping due to $< 1,000$ total reads (`EXCLUDED_LOW_READ_COUNT`), organellar-only database representation, or a highly fragmented reference assembly. |

### Interactive Reports & Data Tables (`results/05_reports/<RUN_TAG>/`)
1. **`validation_report_10panel_<RUN_TAG>.html`**: Self-contained interactive HTML dashboard embedding the **10-Panel Widescreen Diagnostic Suite (Panels A–J)** for every target species—showing the regression scatter plot, pruning trajectory, repeat vs. foreign-contaminant quadrant, proportional read/genome funnels, alignment identity spectrum, genome-wide Manhattan coverage map, 516-week time series, 52-week seasonal phenology curve, and 100-bin $\times$ active-week hierarchical clustermap.
2. **`validation_report_4panel_<RUN_TAG>.html`**: Lighter interactive HTML dashboard embedding the classic $2 \times 2$ **4-Panel Diagnostic Figure** for each species.
3. **`validation_master_summary_<RUN_TAG>.{tsv,parquet}`**: Machine-readable summary table with all pre-flight, assembly, mapping, and regression metrics for every target species in `config/targets.tsv`.
4. **`validated_timeseries_long_<RUN_TAG>.parquet`**: Tidy 516-week time-series table (`raw_kraken_reads`, `all_cr_mapped_reads`, and `validated_clean_bin_reads`) ready for downstream ecological modeling.
