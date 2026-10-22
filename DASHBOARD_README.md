<a id="top"></a>
# 10-Panel Diagnostic Dashboard — Reader's Guide

*A plain-English and technical guide to interpreting the 10-panel per-species diagnostic dashboard in the HTML merged report.*

**HINT**: _press the dashboard figure to zoom and pan_.

---

## 1. Layout at a Glance & Panel Index

```text
┌─────────────────────┬────────────────┬────────────────┬─────────┬──────┐
│ A. Genomic          │ B. Outlier     │ C. Outlier     │ D1.     │ D2.  │
│    Proportionality  │    Pruning     │    Quadrant    │ Read    │ % ID │
│    Regression       │    Trajectory  │    (Repeat vs  │ Funnel  ├──────┤
│    (100 Bins)       │    Curve       │     Artifact)  │         │ D3.  │
├─────────────────────┴────────────────┴────────────────┴─────────┴──────┤
│         [ Shared Colorbar: 516-Week Temporal Cosine Similarity ]       │
├─────────────────────┬────────────────┬────────────────┬────────────────┤
│ E. Genome-Wide CR   │ F. Weekly      │ G. Bin         │                │
│    Coverage Map     │    Signal      │    Temporal    │ J. Bin × Week  │
│    (Manhattan Plot) │    Retention   │    Occupancy   │    Hierarchical│
├─────────────────────┼────────────────┴────────────────┤    Clustermap  │
│ H. 10-Year Weekly   │ I. 10-Year Seasonal Phenology   │    (Rows 2–3)  │
│    Time Series      │    Profile (Weeks 1–52)         │                │
└─────────────────────┴─────────────────────────────────┴────────────────┘
```

### Table of Contents (Click Any Panel to Jump to Its Page)

| Panel | Section Link | Core Question Answered |
| :--- | :--- | :--- |
| **A** | [**Panel A — Genomic Proportionality Regression**](#panel-a) | *Do larger Classifiable Region (CR) bins attract proportionally more reads?* |
| **B** | [**Panel B — Outlier Pruning Trajectory**](#panel-b) | *How many noisy bins had to be removed to reach $R \ge 0.90$, and how much CR DNA remained?* |
| **C** | [**Panel C — Outlier Classification Quadrant**](#panel-c) | *Are high-depth outlier bins endogenous repeats from the target species or foreign cross-mapping?* |
| **D** | [**Panels D1, D2 & D3 — Read Funnel, % Identity & Reference Funnel**](#panel-d) | *How many reads survive each filter, what is their % identity, and how much of the genome is CR?* |
| **E** | [**Panel E — Genome-Wide CR Coverage Map (Manhattan Plot)**](#panel-e) | *What does read density (`Reads / kb`) look like across every individual CR segment?* |
| **F** | [**Panel F — Weekly Signal Retention**](#panel-f) | *Does validation scale weekly counts proportionally, or scrub out false-positive weeks?* |
| **G** | [**Panel G — Bin Temporal Occupancy**](#panel-g) | *In how many of the 516 weeks does each bin appear, and how many active weeks ($N_{\text{eff}}$) drive the signal?* |
| **H** | [**Panel H — 10-Year Longitudinal eDNA Time Series (2013–2022)**](#panel-h) | *When did the species appear over the 516 consecutive weeks before and after validation?* |
| **I** | [**Panel I — 10-Year Seasonal Phenology Profile (Weeks 1–52)**](#panel-i) | *What is the species' average Jan–Dec seasonal cycle, and were off-season noise spikes removed?* |
| **J** | [**Panel J — Bin × Week Hierarchical Clustermap**](#panel-j) | *Do all 100 genomic bins light up together in lockstep across active monitoring weeks?* |

### Figure Header & Verdict Badges
At the top of each figure is the species name, its NCBI TaxID, a color-coded **Validation Status Badge**, and a two-line summary of key genome and read counts:
* **`VALIDATED_CLEAN` (Deep Green)**: Passed validation immediately across all bins ($R \ge 0.90$) with **0 bins pruned**.
* **`VALIDATED_PRUNED` (Emerald Green)**: Passed validation ($R \ge 0.90$) after removing $\ge 1$ noisy outlier bins while keeping $\ge 50\%$ of the Classifiable Regions.
* **`FAILED_LOW_OCCUPANCY` (Red)**: Rejected because $< 50\%$ of bins contained any reads (or initial slope $\le 0$).
* **`FAILED_LOW_CORRELATION` (Red)**: Rejected because pruning down to the $50\%$ CR floor still failed to reach $R \ge 0.90$.

<div id="panel-a" style="page-break-before: always; break-before: page;"></div>

## Panel A — Genomic Proportionality Regression

* **In Plain English**:
  - **What it shows**: This is the core validation test. Each dot is one of the species' **100 genomic bins**, plotted by the **size of its Classifiable Regions on the X-axis** (from the smallest `0.3×` bin on the left to the largest `1.7×` bin on the right) against the **number of mapped eDNA read pairs on the Y-axis**.
  - **How to read it**:
    * If the species is genuinely present, larger bins catch proportionally more reads, forming a clear upward diagonal line.
    * **Circles (`●`)** are **clean, retained bins**; **red-outlined triangles (`▲`, numbered `#1, #2, ...`)** are **pruned outlier bins** removed by the pipeline because they had abnormal read spikes.
    * The **dashed gray line** shows the initial trend line across all 100 bins before pruning; the **solid colored line (with shaded 95% confidence band)** shows the final trend line across the clean retained bins.
    * **Vertical whiskers on each dot**: Show how evenly reads are distributed among the individual CR segments *inside* that bin. Short whiskers mean even coverage across the bin; long whiskers mean a few segments inside the bin have much higher or lower coverage than the rest.
    * **Dot color (purple to yellow-green)**: Shows **Temporal Concordance (Cosine Similarity)**—whether the weeks in which this bin appears over the 10 years match the rest of the clean genome (`1.0` = yellow-green, identical weekly pattern; `0.0` = dark purple, completely different weeks).
    * **Slim strip at the top of Panel A**: The gray bars show how many different reference contigs went into each bin, and the blue line shows the median length (`kb`) of the CR segments inside that bin.
* **Technical Note**:
  - **Leave-one-out Cosine Similarity**: Each bin $i$'s 516-week read vector is compared via cosine similarity against the sum of all *other* retained clean bins ($\mathbf{C}_{-i} = \sum_{j \in \mathcal{R} \setminus \{i\}} \mathbf{v}_j$) so a high-read bin cannot inflate its own similarity score.
  - **Whiskers (`p25–p75`)**: Represent the length-weighted 25th and 75th percentiles of per-kb read density across the CR segments inside each bin, scaled to the bin's width so they can be plotted directly on the read-count axis.

[↑ Back to Panel Index](#top)

<div id="panel-b" style="page-break-before: always; break-before: page;"></div>

## Panel B — Outlier Pruning Trajectory

* **In Plain English**:
  - **What it shows**: A step-by-step history of the outlier pruning process (`Step 0, 1, 2, ...` on the X-axis).
  - **How to read it**:
    * The **blue curve (left Y-axis)** tracks the linear correlation (**Pearson $R$**) as noisy bins are removed one by one, climbing toward the **green dashed validation threshold ($R = 0.90$)**.
    * The **orange curve (right Y-axis)** tracks the **percentage of Classifiable Region (CR) DNA remaining**, which is not allowed to drop below the **red dotted safety floor ($50\%$)**.
    * For a true positive with a few repeat spikes, the blue curve jumps sharply above $0.90$ after removing just a handful of bins while the orange curve stays high ($>85\%$). For a false positive, the blue curve stays low until the orange curve hits the $50\%$ floor and fails.
* **Technical Note**:
  - At each step, Ordinary Least Squares (OLS) regression is fit to the currently retained bins, and the single bin with the largest absolute vertical residual ($|Y_i - \hat{Y}_i|$) is removed—catching both extreme upward repeat spikes and extreme downward coverage dropouts.

[↑ Back to Panel Index](#top)

<div id="panel-c" style="page-break-before: always; break-before: page;"></div>

## Panel C — Outlier Classification Quadrant

* **In Plain English**:
  - **What it shows**: Explains *why* certain bins had too many reads by separating **repeats that belong to the target species itself** from **foreign reads cross-mapping from a different organism**.
  - **How to read it**:
    * **X-axis (`Excess Read Depth`)**: How many more (or fewer) reads a bin caught compared to what we would expect from its size alone (`0` = exact expected read depth; `+1` = $2\times$ expected; `+2` = $4\times$; `+3` = $8\times$).
    * **Y-axis (`Temporal Concordance / Cosine Similarity`)**: Whether the bin's reads appear in the same weeks of the 10-year time series as the rest of the species' clean genome (`1.0` = same weeks; `0.0` = different weeks).
    * This splits high-coverage outlier bins (`▲` on the right side) into two meaningful categories:
      1. **Upper-Right (`Endogenous Multi-Copy Repeat`)**: The bin has way too many reads, **but they show up in the exact same weeks as the target species**. This points to a multi-copy repeat (like ribosomal DNA or a transposable element) that is genuinely part of the target organism's genome.
      2. **Lower-Right (`Foreign Cross-Mapping Artifact`)**: The bin has way too many reads, **and they show up in completely different weeks from the target species**. This points to foreign DNA from another airborne microbe or assembly contamination cross-mapping to a conserved spot.
    * **Dot size (`Gini Inequality`)**: Indicates whether the bin's reads are spread evenly across the bin (**small dot**) or piled into one tiny hotspot inside the bin (**large dot**).
* **Technical Note**:
  - **Expected Depth**: Computed from the length-proportional clean read density ($\bar{d}_{\text{clean}} \times \text{bin\_bp}$) rather than the unconstrained OLS line so negative intercepts cannot distort fold-changes in small bins.
  - **Gini Coefficient ($G_b \in [0, 1]$)**: Length-weighted Lorenz inequality of read density across all CR segments within bin $b$ ($0$ = perfectly uniform density across segments; $1$ = all reads concentrated in a single base pair).

[↑ Back to Panel Index](#top)

<div id="panel-d" style="page-break-before: always; break-before: page;"></div>

## Panels D1, D2 & D3 — Read Attrition Funnel, Alignment Identity Spectrum & Reference Funnel

* **In Plain English**:
  - **Panel D1 (`eDNA Read Attrition` — left funnel)**:
    * Shows what fraction of the species' initial Kraken 2 reads survive each quality filter, from top to bottom:
      1. **`1. Raw Kraken Clade Pairs`**: All read pairs assigned to the species clade by Kraken 2.
      2. **`2. CR Mapped (MAPQ ≥ 30)`**: Pairs that align confidently (`minimap2`) to the species' Classifiable Regions.
      3. **`3. Passed ≥ 80% Aligned`**: Pairs where at least 80% of the read length aligns without being clipped off at the ends.
      4. **`4. Passed ≥ 99% Identity`**: Pairs matching the reference genome with at least 99% sequence identity (or the active run's identity cutoff).
      5. **`5. Clean Validated Bins`**: Passing read pairs that land inside the final **retained clean bins**.
    * The width of each bar is drawn strictly to scale with the percentage of reads remaining.
  - **Panel D2 (`Alignment Identity Spectrum` — top-right histogram)**:
    * Shows the sequence match percentage (`<90%` up to `100%` identity) of all read pairs that aligned to the Classifiable Regions (Stage 3 of the funnel).
    * **Red bars (left of the dashed cutoff line)**: Reads rejected because they had too many mismatches. A big hump around `93%–97%` typically indicates reads from a closely related sister species that isn't in the reference database.
    * **Green / Orange stacked bars (right of the cutoff line)**: Reads passing the identity cutoff, split into **clean retained bins (green)** vs. **pruned outlier bins (orange)**.
  - **Panel D3 (`Reference Sequence Funnel` — bottom-right funnel)**:
    * Shows how much of the species' reference genome is actually usable for classification (drawn strictly to scale):
      1. **`1. Total Genome Assembly`**: Full size of the reference genome (`100%`).
      2. **`2. Classifiable Regions (CR)`**: The portion of the genome unique enough for Kraken 2 to assign to this species clade.
      3. **`3. Clean CRs (Post-Pruning)`**: The portion of the genome remaining after pruning outlier bins.
* **Technical Note**:
  - In **Panel D2**, each aligned read pair is placed into a single integer identity bin (`<90%`, `90%`, ..., `100%`) using the **lower identity of its two mates** ($\min(\text{ID}_{R1}, \text{ID}_{R2})$), since both mates must pass the identity cutoff for a read pair to be accepted in Stage 4. This makes the total height of every bar in D2 invariant across runs that only change the identity threshold.

[↑ Back to Panel Index](#top)

<div id="panel-e" style="page-break-before: always; break-before: page;"></div>

## Panel E — Genome-Wide CR Coverage Map (Manhattan Plot)

* **In Plain English**:
  - **What it shows**: Zooms inside the 100 bins to show the read density (`Reads / kb`) of **every individual Classifiable Region (CR) segment** across the entire genome from `bin_001` (left) to `bin_100` (right).
  - **How to read it**:
    * Each dot is one CR segment, colored by its 516-week **Temporal Concordance** (using the shared colorbar directly above Panel E).
    * **Vertical pink/rose shaded bands** highlight the regions belonging to **pruned outlier bins**—letting you see at a glance whether a pruned bin contained one towering spike or a broad cluster of noisy segments.
    * **Shaded gray band at the bottom (`0` to `1 Read / kb`)**: Uses a standard linear scale so zero-read CRs are visible, while everything above `y = 1` uses a logarithmic scale (`10`, `100`, `1,000 Reads / kb`) so huge repeat spikes fit on the chart.
* **Technical Note**:
  - **Cosine Similarity**: Each individual CR segment's 516-week profile is compared directly against the total 516-week signal of all retained clean bins (after pruning).
  - **Y-Axis Scale**: Uses a symmetrical logarithmic (`symlog`, `linthresh = 1.0`) Y-axis transformation—linear on $[0, 1]$ and $\log_{10}$ above $1\text{ Read/kb}$.

[↑ Back to Panel Index](#top)

<div id="panel-f" style="page-break-before: always; break-before: page;"></div>

## Panel F — Weekly Signal Retention

* **In Plain English**:
  - **What it shows**: Compares the weekly read counts **before validation (X-axis)** against **after validation (Y-axis)** across all 516 weekly samples.
  - **How to read it**:
    * Each dot is one week (`516` dots total).
    * **Dashed diagonal line (`Y = X`)**: Where weeks would sit if 100% of raw Kraken 2 reads survived validation.
    * **Solid dark curve (`Y = c · X`)**: Where weeks sit if validation retains a steady proportion ($c$) of reads every week.
    * Dots sitting tightly along the solid curve mean the weekly abundance pattern was preserved; dots sitting flat along the bottom (`Y = 0` despite high `X`) are **false-positive weeks that were completely scrubbed out** by validation.
* **Technical Note**:
  - Plotted on $\log_{10}(\text{Reads} + 1)$ axes and reports three correlations in the legend: untransformed linear Pearson $R$ (`Lin R`), log-space Pearson $R_{\log}$ (`R_log`), and non-parametric Spearman rank correlation ($\rho$).

[↑ Back to Panel Index](#top)

<div id="panel-g" style="page-break-before: always; break-before: page;"></div>

## Panel G — Bin Temporal Occupancy

* **In Plain English**:
  - **What it shows**: Plots how frequently each of the 100 bins appears across the 10-year dataset (**X-axis: `% of the 516 weeks` where the bin caught at least 1 read pair**) against its **Temporal Concordance (Y-axis: Cosine Similarity)**.
  - **How to read it**:
    * **Why this matters**: Tells you whether the species appears across many weeks of the 10-year series or only in 1 or 2 sudden bursts.
    * Check the **`Effective active weeks (N_eff)`** number in the subtitle:
      - If a species (like a rare migrating insect) has 95% of all its 10-year reads in just **2 summer weeks** (`N_eff ≈ 2w`), almost every bin will show a Cosine Similarity near `1.0` simply because all reads came from those same 2 weeks. Panel G makes that immediately visible by showing all dots clustered at very low weekly occupancy (`< 2%`) on the far left.
      - Conversely, for a common species present across hundreds of weeks (`N_eff > 30w`), clean bins sit high and to the right, while any low-concordance outlier bins stand out clearly near the bottom.
* **Technical Note**:
  - **Effective Active Weeks ($N_{\text{eff}}$)**: Calculated as the inverse Simpson diversity index of the clean 516-week time series, $N_{\text{eff}} = 1 / \sum_{t=1}^{516} p_t^2$, where $p_t$ is the fraction of total clean reads occurring in week $t$ ($N_{\text{eff}} = 1$ if all reads occur in a single week; $N_{\text{eff}} = 516$ if reads are equal across all 516 weeks).

[↑ Back to Panel Index](#top)

<div id="panel-h" style="page-break-before: always; break-before: page;"></div>

## Panel H — 10-Year Longitudinal eDNA Time Series (2013–2022)

* **In Plain English**:
  - **What it shows**: The complete 10-year weekly abundance timeline across all **516 consecutive weeks** from 2013 to 2022 (alternating white and gray background bands mark each calendar year).
  - **How to read it**:
    * **Light gray shaded background (`Raw Kraken Clade Pairs`)**: Unfiltered Kraken 2 read counts.
    * **Blue curve (`All-CR Mapped Pairs`)**: Reads that passed stringent `minimap2` alignment to the Classifiable Regions *before* outlier bin pruning.
    * **Green/red curve (`Validated Clean-Bin Pairs`)**: Final validated signal from the **clean retained bins only**—showing how false-positive spikes are stripped away while true seasonal peaks remain.

[↑ Back to Panel Index](#top)

<div id="panel-i" style="page-break-before: always; break-before: page;"></div>

## Panel I — 10-Year Seasonal Phenology Profile (Weeks 1–52)

* **In Plain English**:
  - **What it shows**: Folds the 10 years of data into a single **average annual calendar (January / Week 1 to December / Week 52)** to show the species' seasonal life cycle (phenology).
  - **How to read it**:
    * **Dashed gray curve (`Pre-Validation Raw`)**: Seasonal rhythm of the raw Kraken 2 reads.
    * **Solid colored curve & shaded area (`Post-Validation Clean`)**: Seasonal rhythm of the validated clean-bin reads.
    * **Pink shaded gaps between the two curves**: Highlight off-season noise (for example, winter cross-mapping spikes) that was present in the raw Kraken 2 data and removed by validation.
    * **Left vs. Right Y-axis**: The left Y-axis shows the **relative share (`%`) of the species' annual reads** occurring in each week of the year (so you can compare the seasonal shape before and after validation directly), while the right Y-axis shows the **actual average number of clean read pairs per week**.
* **Technical Note**:
  - Applies a 3-week circular smoothing kernel (`[0.25, 0.50, 0.25]`) that wraps smoothly across the Week 52 $\leftrightarrow$ Week 1 New Year boundary, and reports the Pearson correlation ($r$) between the pre- and post-validation 52-week curves.

[↑ Back to Panel Index](#top)

<div id="panel-j" style="page-break-before: always; break-before: page;"></div>

## Panel J — Bin × Week Hierarchical Clustermap

* **In Plain English**:
  - **What it shows**: A bird's-eye heatmap of **all 100 genomic bins (rows)** across the **active monitoring weeks (columns, ordered chronologically from 2013 to 2022)**.
  - **How to read it**:
    * Darker blue cells mean higher read density (`Reads / kb`) in that bin during that week; pale yellow means zero reads.
    * A clustering tree (**dendrogram on the left**) automatically groups bins that behave similarly over time next to each other.
    * **Three narrow vertical strips** between the tree and the heatmap show each bin's **Temporal Concordance (`CosSim`)**, **Weekly Occupancy (`Occ%`)**, and **Validation Status (`Retained` in green vs. `Pruned` in red)**.
    * **What to look for**:
      - **True biological presence**: Clean **vertical stripes** running top-to-bottom across the green retained bins—meaning the whole genome lights up together in the exact same weeks.
      - **Contaminants / cross-mapping spikes**: Isolated **horizontal streaks** (flagged in red in the `Status` strip) that light up in weeks where the rest of the genome is completely blank.
* **Technical Note**:
  - **Normalization & Clustering**: Row values are length-normalized log read densities, $\log_{10}(\text{Reads/kb} + 1)$, so large bins (`1.7×`) and small bins (`0.3×`) are compared on equal footing. Rows are clustered using **Euclidean distance** and **Ward's minimum-variance linkage** with optimal leaf ordering, which groups bins by both their seasonal timing and coverage intensity (up to the top 65 active weeks by read count are displayed in chronological order).

[↑ Back to Panel Index](#top)
