# ValidateClassificationsSMK — Instructional Video Storyboard & AI Prompt Pack

- **Target Runtime**: ~105 seconds (10 scenes, 8–12 seconds per scene)
- **Target Audience**: Scientific peers / conference presentation
- **Visual Style**: Clean 2D flat vector scientific schematic on a crisp publication background (`#F8FAFC`), inspired by animated *Nature Methods* / BioRender figures (zero photorealism).
- **Narrative Arc**:
  1. Physical SweBITS air-filter stations & random whole-genome DNA shedding
  2. The Kraken2 specificity dilemma (uniform coverage vs. cross-reactive repeat spikes & unclassifiable genome gaps)
  3. In-silico `400 bp` sliding-window reference shearing (`200 bp` step, $2\times$ coverage), Kraken2 classification, and `200 bp` flanking extension to define **Classifiable Regions (CRs)**
  4. Equal-sized vs. variable-size linear binning (`100` bins scaled `0.3×` to `1.7×` mean CR length)
  5. Stringent read mapping (`minimap2`) and iterative leave-one-out OLS outlier pruning ($R \ge 0.90$, $\ge 50\%$ retained CR safety floor)
  6. Pre-validation vs. post-validation 516-week SweBITS time series and final True/False Positive verdict

---

## 0. Global Style Lock & Color Palette

Use this consistent four-color schematic palette across all tools:
- **Background**: Crisp Off-White (`#F8FAFC`) with a subtle light-gray geometric grid (`#E2E8F0`)
- **Deep Slate Navy (`#1E293B`)**: Reference genome backbones, coordinate axes, and OLS regression line
- **Sky Blue (`#0EA5E9`)**: Airborne eDNA fragments and SweBITS air-filter streamlines
- **Emerald Teal (`#10B981`)**: Classifiable Regions (CRs), retained bins, and validated seasonal time-series signals
- **Coral Red (`#EF4444`)**: Cross-reactive repeat spikes, pruned outlier bins, and raw pre-validation noise spikes

**Reusable Style Prefix (already pre-pended in Section 3 below):**
> `2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B) for genome backbones and axes, Sky Blue (#0EA5E9) for airborne eDNA reads, Emerald Teal (#10B981) for Classifiable Regions and validated data, and Coral Red (#EF4444) for cross-reactive outlier spikes.`

---

## 1. Storyboard, Text Overlays & Voiceover Script

| Scene | Time | Pipeline Stage | Visual Schematic Action | On-Screen Text Overlay (Post-Edit) | Voiceover (VO) Script |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | `0:00–0:10` | **SweBITS Air Filters & Random Shedding** | Schematic air-filter station draws in blue wind streamlines; a stylized cellular silhouette fragments randomly into dozens of short Sky Blue DNA bars collected evenly on a filter grid. | **SweBITS Air-Filter Stations**<br>True organism presence $\rightarrow$ random whole-genome DNA shedding | *"Across Sweden, SweBITS air-filter stations capture airborne DNA over decades. When a species is truly present in the environment, its cells break down randomly—shedding DNA fragments evenly across its entire genome."* |
| **2** | `0:10–0:20` | **The Kraken2 Dilemma (True vs. False Positive)** | Split horizontal genome bars: Top bar shows Sky Blue reads landing uniformly; Bottom bar shows Coral Red reads piling up into a single towering spike over a conserved repeat while the rest of the genome is empty. | **The Taxonomic Classification Trap**<br>Uniform Coverage (True) vs. Isolated Repeat Spike (False Positive) | *"Fast k-mer classifiers like Kraken2 assign millions of reads, but face a critical trap: conserved genes or repeats shared with other microbes create massive, localized read spikes that mimic species presence."* |
| **3** | `0:20–0:30` | **Why Not the Whole Genome?** | A single Navy reference genome bar with greyed-out "blind spots" (conserved/low-complexity regions) where blue reads bounce off, and bright Emerald Teal zones where species-specific reads stick. | **Not All Regions Are Classifiable**<br>We must measure uniformity *only* where Kraken2 can actually see the species | *"However, even for a true positive, reads won't map uniformly across the whole assembly—because Kraken2 can only classify regions that are unique to that species."* |
| **4** | `0:30–0:42` | **In-Silico Shearing & Kraken2 Scan** | A digital cutter slices the Navy reference genome via a sliding `400 bp` window (`200 bp` step, $2\times$ coverage) into an overlapping cascade of pseudo-reads; a scanner beam turns species-specific fragments Emerald Teal and discards grey ones. | **Stage 1–2: In-Silico Genome Shearing**<br>`400 bp` sliding window (`200 bp` step) $\rightarrow$ Kraken2 target clade filter | *"To map these classifiable zones, the pipeline first shreds the reference genome in silico—sliding a 400-base-pair window in 200-base-pair steps—and runs every pseudo-read through Kraken2."* |
| **5** | `0:42–0:52` | **Merging & 200 bp Flanking Extension (CRs)** | Overlapping Emerald Teal pseudo-reads fuse into solid blocks; arrows on both ends expand each block outward by `200 bp` so incoming paired-end eDNA reads spanning the boundary fit completely inside. | **Defining Classifiable Regions (CRs)**<br>Merge windows + `200 bp` flanking extension for real mapped reads | *"Classified windows are merged into contiguous intervals and padded by 200 base pairs on each flank—ensuring real sequencing reads that overlap the edges can align end-to-end. These form the species' Classifiable Regions."* |
| **6** | `0:52–1:04` | **Equal vs. Linear Binning (`0.3x` to `1.7x`)** | Left: Equal-sized Teal blocks drop onto an X-Y plot at the exact same X-value, forming a useless vertical column of dots. Right: Contigs sort shortest-to-longest and slice into 100 linearly growing bins (`0.3x` to `1.7x`), spreading dots horizontally across the X-axis. | **Stage 3: Variable-Size Linear Binning**<br>100 bins scaled `0.3×` $\rightarrow$ `1.7×` mean CR length ($5.7\times$ X-axis spread) | *"To test whether read count scales linearly with region size, equal-sized bins fail—they collapse onto a single X-axis value. Instead, we sort contigs shortest-to-longest and slice Classifiable Regions into 100 linearly growing bins, from 0.3x to 1.7x mean length."* |
| **7** | `1:04–1:14` | **Mapping eDNA & Initial Skewed Regression** | Sky Blue eDNA reads from the air filter rain down (`minimap2`) onto the 100 graduated Teal bins. An X-Y scatter plot (`Bin CR Length` vs. `Read Count`) forms: 95 Teal dots follow a diagonal, but 5 Coral Red outlier dots tower at the top, pulling the dashed regression line upward. | **Stage 4–5: Read Mapping & Raw OLS Fit**<br>Cross-reactive repeat bins skew linear regression ($R < 0.90$) | *"Next, real air-filter reads are stringently mapped to the 100 bins. On a plot of read count versus bin size, true coverage forms a clean diagonal—but cross-reactive repeats shoot up as extreme red outliers, skewing the regression line."* |
| **8** | `1:14–1:26` | **Iterative Leave-One-Out Pruning** | One by one, the 5 high Coral Red outlier dots pulse, turn hollow, and fade out. With each removal, the regression line pivots downward and tightens around the Teal dots until a counter hits `R ≥ 0.90`, while a side progress bar stays green (`Retained CR: 94% ≥ 50%`). | **Iterative Leave-One-Out Pruning**<br>Prune worst outlier bin until $R \ge 0.90$ (requiring $\ge 50\%$ CR retained) | *"The pipeline iteratively prunes the worst outlier bin one at a time. As those contaminated regions are removed, the regression line snaps onto the true genome-wide slope until Pearson R reaches 0.90—while safeguarding that at least half of the classifiable genome remains."* |
| **9** | `1:26–1:36` | **Pre- vs. Post-Validation Time Series** | The cleaned bins project onto a 516-week horizontal timeline. First, a jagged Coral Red curve with erratic winter/noise spikes ("Pre-Validation") is shown; as the pruned bins vanish, the red spikes melt away, revealing a smooth, repeating Emerald Teal seasonal wave ("Post-Validation"). | **10-Year SweBITS Time Series (516 Weeks)**<br>Coral Red: Raw unvalidated noise $\rightarrow$ Emerald Teal: Clean seasonal signal | *"Across the 516-week SweBITS time series, the impact is unmistakable: pruning removes erratic false-positive spikes in red, uncovering the clean, biologically coherent seasonal rhythm of the species in teal."* |
| **10** | `1:36–1:45` | **Final Verdict Split-Screen** | Split screen: Left card (**VALIDATED TRUE POSITIVE**) shows a tight Teal linear slope ($R = 0.95$) and clean seasonal wave. Right card (**REJECTED FALSE POSITIVE**) shows only a few isolated red spikes near zero and a flat/failed regression. | **Automated Validation Verdict**<br>True Positive ($R \ge 0.90$) vs. Rejected False Positive | *"Species with genuine whole-genome coverage pass as validated true positives, while classifications driven only by isolated repeat spikes are automatically rejected—turning noisy k-mer hits into high-confidence biodiversity data."* |

---

## 2. Master Unified Prompt (For End-to-End AI Video / Storyboard Platforms)

```text
Create a 105-second 2D flat-vector scientific explainer animation (clean off-white #F8FAFC background, subtle technical grid, zero photorealism, Nature Methods / BioRender infographic style) explaining the "ValidateClassificationsSMK" bioinformatics pipeline for validating environmental DNA (eDNA) detections from SweBITS air-filter stations.

Use a strict, consistent four-color schematic palette throughout:
- Deep Slate Navy (#1E293B) for reference genome backbones and plot axes
- Sky Blue (#0EA5E9) for airborne eDNA sequencing reads
- Emerald Teal (#10B981) for Classifiable Regions (CRs), retained bins, and validated seasonal signals
- Coral Red (#EF4444) for cross-reactive repeat spikes and pruned outlier bins

Follow this exact 10-step visual progression:
1. SWEBITS AIR FILTERS (0:00-0:10): A minimalist vector air-filter station pulls in blue air streamlines. A stylized biological cell fragments randomly into short Sky Blue DNA bars that distribute evenly across the filter mesh.
2. THE KRAKEN2 DILEMMA (0:10-0:20): Two horizontal Navy genome bars compare a True Positive (Sky Blue reads landing evenly across the bar) against a False Positive (Coral Red reads piling into a single tall vertical spike over a conserved repeat region).
3. CLASSIFIABLE VS. UNCLASSIFIABLE GENOME (0:20-0:30): On a single Navy genome bar, greyed-out conserved sections deflect reads, while unique Emerald Teal segments capture them—showing why uniformity must be measured only within Kraken2-classifiable regions.
4. IN-SILICO GENOME SHEARING (0:30-0:42): A geometric caliper slides along the Navy reference genome in a 400 bp window (200 bp step), slicing it into a cascade of overlapping tiles. A Kraken2 scanner beam highlights target-species tiles in Emerald Teal and fades out non-specific grey tiles.
5. MERGING & 200 BP EXTENSION (0:42-0:52): The Emerald Teal tiles merge into solid intervals, and outward arrows pad both ends by 200 bp so full paired-end eDNA reads spanning the edges align cleanly. These form the Classifiable Regions (CRs).
6. EQUAL VS. LINEAR BINNING (0:52-1:04): First, equal-sized Teal bins drop onto an X-Y chart at the same X coordinate, forming a vertical stack of dots with zero horizontal spread. Then, the CRs sort shortest-to-longest and slice into 100 progressively wider bins (scaled linearly from 0.3x to 1.7x mean length), spreading the dots horizontally across the X-axis.
7. READ MAPPING & SKEWED REGRESSION (1:04-1:14): Sky Blue eDNA reads map onto the 100 graduated bins. On a scatter plot of Bin CR Length (X) vs. Read Count (Y), 95 Teal dots follow a linear diagonal, but 5 Coral Red outlier dots sit high above the trend, pulling the regression line off-course.
8. ITERATIVE OUTLIER PRUNING (1:14-1:26): One by one, the 5 Coral Red outlier dots pulse and vanish. With each step, the regression line pivots downward and locks tightly onto the remaining Emerald Teal dots (R >= 0.90), while a vertical retention gauge stays green above the 50% minimum CR threshold.
9. PRE- VS. POST-VALIDATION TIME SERIES (1:26-1:36): The plot transitions into a 10-year (516-week) time-series chart. A noisy Coral Red curve with sharp erratic spikes ("Pre-Validation") morphs as the pruned bins disappear, leaving a clean, rhythmic Emerald Teal seasonal wave ("Post-Validation").
10. FINAL VERDICT SPLIT-SCREEN (1:36-1:45): Split-screen finale comparing a "Validated True Positive" on the left (tight Teal linear scatter plot + seasonal wave, green checkmark) against a "Rejected False Positive" on the right (scattered Coral Red spikes + collapsed regression, red cross).
```

---

## 3. Ready-to-Paste Scene-by-Scene Prompts (For Veo, Sora, Runway, Kling)

Each prompt below is self-contained (with the global style lock included) so you can paste them one by one into a generative video model.

### Scene 1 — SweBITS Air Filters & Random Shedding (`0:00–0:10`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). A minimalist 2D schematic of an environmental air-sampling station on the left with smooth Sky Blue (#0EA5E9) wind streamlines flowing left-to-right through a vertical geometric filter grid. Inside the airflow, a simple green botanical/microbial circle dissolves evenly along a long Deep Slate Navy (#1E293B) horizontal bar into dozens of uniform, tiny Sky Blue rectangular DNA fragments that land evenly across the filter grid. Smooth flat vector animation, constant speed.
```

### Scene 2 — The Kraken2 Dilemma: Uniformity vs. Spike (`0:10–0:20`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). Split horizontal view comparing two Deep Slate Navy (#1E293B) genome bars on an off-white grid. On the top bar, dozens of small Sky Blue (#0EA5E9) rectangles rain down and settle at an even, uniform height across the entire length. On the bottom bar, dozens of Coral Red (#EF4444) rectangles fall and stack exclusively onto one narrow spot, building a single tall, sharp vertical red tower while the rest of the bottom bar remains completely empty.
```

### Scene 3 — Unclassifiable Blind Spots (`0:20–0:30`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). Close-up orthographic view of a single horizontal Deep Slate Navy (#1E293B) genome bar divided into alternating segments of muted light-gray hatched zones and vibrant Emerald Teal (#10B981) solid zones. Small Sky Blue (#0EA5E9) read rectangles fall from above: the ones hitting the gray hatched zones fade away transparently, while the ones landing on the Emerald Teal zones lock smoothly into place.
```

### Scene 4 — In-Silico 400 bp Shearing & Kraken2 Classification (`0:30–0:42`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). A horizontal Deep Slate Navy (#1E293B) genome bar across the top of the frame. A bright geometric bracket slides smoothly from left to right along the bar, peeling off a dense diagonal cascade of overlapping short rectangular tiles below it like a waterfall of stepped bars. A vertical scanning line sweeps across the cascading tiles, turning groups of target tiles bright Emerald Teal (#10B981) while non-target tiles fade to faint gray and disappear.
```

### Scene 5 — Merging & 200 bp Flanking Extension (`0:42–0:52`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). Clusters of overlapping Emerald Teal (#10B981) rectangular tiles align horizontally on a Deep Slate Navy (#1E293B) track and fuse together into solid, seamless Emerald Teal blocks. Once fused, small geometric arrows appear at the left and right edges of each Emerald Teal block, smoothly stretching both ends outward by a fixed margin. Pairs of Sky Blue (#0EA5E9) read bars then drop down and fit neatly inside the expanded teal boundaries.
```

### Scene 6 — Equal Bins vs. `0.3x–1.7x` Variable Linear Binning (`0:52–1:04`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). A 2D coordinate axis in Deep Slate Navy (#1E293B). First, ten equal-width Emerald Teal (#10B981) blocks drop onto the X-axis at the exact same horizontal position, stacking their data dots into a single vertical line where a regression line spins uncertainly. Then, the Emerald Teal blocks rearrange along the top into a smooth staircase of 100 segments that grow progressively wider from left to right—narrow slivers on the left to wide bars on the right—causing the dots on the chart below to fan out horizontally across the entire X-axis.
```

### Scene 7 — Stringent Read Mapping & Skewed OLS Fit (`1:04–1:14`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). A clean 2D scatter plot on an off-white grid with Deep Slate Navy (#1E293B) axes. Along the bottom X-axis, graduated Emerald Teal (#10B981) bars increase in width from left to right. Data points populate the chart: ninety-five Emerald Teal circles form a neat upward diagonal band from bottom-left to top-right, while five bright Coral Red (#EF4444) circles shoot high up near the top of the chart. A dashed Navy regression line tries to fit all points, getting pulled steeply upward out of alignment by the five red circles.
```

### Scene 8 — Iterative Leave-One-Out Outlier Pruning (`1:14–1:26`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). Close-up of the 2D scatter plot with ninety-five Emerald Teal (#10B981) circles along a diagonal and five high Coral Red (#EF4444) outlier circles near the top. One by one, in sequential steps, each Coral Red circle pulses with a thin ring and vanishes. With each red circle removed, the straight Navy regression line pivots smoothly downward, finally snapping tightly through the center of the Emerald Teal diagonal band, which glows softly. On the right side, a vertical Emerald Teal retention bar steps down slightly but stays well above a horizontal midpoint line.
```

### Scene 9 — Pre-Validation vs. Post-Validation 516-Week Time Series (`1:26–1:36`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). A wide horizontal time-series chart spanning left to right on a clean off-white grid. Initially, a jagged Coral Red (#EF4444) line graph displays chaotic, tall vertical spikes scattered across the timeline. Smoothly, a translucent wipe sweeps across the chart as the sharp Coral Red spikes shrink and dissolve, revealing a smooth, rhythmic Emerald Teal (#10B981) wave with regular seasonal peaks and valleys underneath.
```

### Scene 10 — Final Verdict: Validated True Positive vs. Rejected False Positive (`1:36–1:45`)
```text
2D flat vector scientific motion graphics on a clean off-white background (#F8FAFC) with a subtle light-gray geometric grid. Minimalist publication-figure aesthetic (Nature Methods / BioRender style), crisp geometric lines, smooth orthographic camera motion, zero photorealism, no 3D shading, no text or gibberish letters rendered inside the animation. Strict color palette: Deep Slate Navy (#1E293B), Sky Blue (#0EA5E9), Emerald Teal (#10B981), and Coral Red (#EF4444). A balanced 2D split-screen composition on an off-white background divided by a crisp vertical line. On the left side, an Emerald Teal (#10B981) framed panel displays a tight diagonal scatter plot with a solid regression line above a smooth seasonal wave, accented by a clean geometric green checkmark. On the right side, a Coral Red (#EF4444) framed panel displays an empty scatter plot with only three isolated red dots near the origin and a flat dashed line, accented by a geometric red cross.
```
