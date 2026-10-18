# NCBI Assembly Lifecycle & Suppression Analysis

## Overview
During the execution of **ValidateClassificationsSMK (v2)**, an assembly selection failure occurred for *Columba palumbus* (TaxID `8934`, accession `GCA_965250375.1`) because the assembly was suppressed on NCBI/EBI after the SweBITS Kraken 2 database (`master_db_2025`) was constructed in mid-2025.

To understand the scope and impact of assembly deprecation across the entire pipeline, we performed an empirical audit across all **107,575 GenBank and RefSeq assemblies** present in `input/database_manifest.tsv` by streaming them through the live NCBI `datasets summary genome accession` CLI.

---

## 1. Global Audit Results (107,575 Assemblies)

| Status | GenBank (`GCA_`) | RefSeq (`GCF_`) | Total Assemblies | Percentage | Downloadable via `datasets`? |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`current`** | 92,363 | 14,655 | **107,018** | **99.48%** | **Yes** (Primary target) |
| **`previous`** | 376 | 5 | **381** | **0.35%** | **Yes** (Downloadable by exact accession) |
| **`suppressed`** | 105 | 69 | **174** | **0.16%** | **No** (FASTA omitted by NCBI `datasets`) |
| **`not_returned`** | 2 | 0 | **2** | **0.00%** | **No** (Accession unrecognized) |
| **Total** | **92,846** | **14,729** | **107,575** | **100.00%** | — |

**Key Takeaway**: **99.83%** of all database assemblies (current + previous) remain immediately downloadable from NCBI today. Only **0.16% (174 assemblies)** are suppressed.

---

## 2. In-Depth Analysis of Suppressed Assemblies

Of the 174 suppressed assemblies, two primary categories account for almost 50% of all suppressions:

### Category A: Standard Eukaryotic Genome Re-Annotation (66 Assemblies)
- **Manifest Source**: **100% (66 / 66) are RefSeq (`GCF_`) accessions**.
- **NCBI Reason**: *"This record was removed as a result of standard genome annotation processing. Please see www.ncbi.nlm.nih.gov/genome/annotation_euk/process/ for more information."*
- **Why this happens**: RefSeq is a secondary curated database. When NCBI computes a new annotation release for an organism (updating gene models with newer RNA-seq data or selecting a newer representative assembly), the older `GCF_` record is retired/suppressed in RefSeq.
- **Underlying Primary Sequence (GenBank `GCA_`) Status**:
  When we queried the corresponding GenBank accessions (`GCA_...` instead of `GCF_...`) for these 66 organisms:
  - **49 (74.2%)** of the underlying GenBank assemblies are **still `current` and fully downloadable**!
  - **10 (15.2%)** are `previous` (still downloadable by exact accession).
  - Only **2** are suppressed in GenBank.
- **Recommended Grand-Scale Workaround**:
  Because RefSeq eukaryotic assemblies are built directly from their GenBank counterparts (identical nucleotide coordinates and contigs, differing only in NCBI-curated feature annotation), an automated fallback rule can convert a suppressed `GCF_XXXXXXXXX.N` accession to `GCA_XXXXXXXXX.N`. This immediately recovers **89.4%** of these suppressed assemblies with zero sequence discrepancy.

---

### Category B: Submitter Takedown & Replacement Requests (20 Assemblies)
- **Manifest Source**: **100% (20 / 20) are GenBank (`GCA_`) accessions**.
- **NCBI Reason**: *"This record was removed at the submitter's request. Please contact update@ebi.ac.uk for further details."*
- **Submitter Pattern**:
  - **16 of the 20 (80%)** were submitted to **EBI** (`update@ebi.ac.uk`) by the **Darwin Tree of Life (DToL) Project** at the **Wellcome Sanger Institute**.
  - **Taxonomic Focus**: Almost exclusively British and European bird species (*Columba palumbus*, *Anas acuta*, *Anser fabalis*, *Branta canadensis*, *Cygnus columbianus*, *Grus grus*, *Larus argentatus*, *Numenius arquata*, *Phalacrocorax carbo*, *Bucephala clangula*, *Clangula hyemalis*, *Caprimulgus europaeus*).
- **Why Sanger Requested Removal**:
  DToL initially released PacBio + Arima Hi-C draft assemblies in 2023–2024 (e.g. `GCA_965250375.1` for *Columba palumbus*). In late 2025/2026, Sanger re-sequenced these specimens with Oxford Nanopore ultra-long reads, performed manual TreeVal chromosome curation, and published **phased dual-haplotype assemblies** (`hap1.2` and `hap2.2`, e.g. `GCA_988290405.1` and `GCA_988291125.1`). Upon releasing the new haplotype assemblies, Sanger requested EBI to suppress the older v1 pseudohaplotype records.
- **Replacement Assessment for Grand-Scale Runs**:
  - For each of these organisms, live replacement genomes exist on NCBI.
  - *Feasibility of using replacement genomes*: Because the Kraken 2 database was indexed against the v1 assemblies in 2025, using a v2 assembly introduces minor coordinate shifts (e.g., in *C. palumbus*, chromosome 30 was split into chromosomes 30 and 35 in v2). However, for species-level Classifiable Region (CR) discovery and eDNA read mapping, synteny and k-mer content between v1 and v2 are typically >98% identical.
  - *Hierarchy of resolution*:
    1. If other current non-suppressed assemblies exist within the JolTax species clade, select the best live assembly (already implemented in Stage 0).
    2. If all clade assemblies are suppressed, check if an official replacement accession is listed in `assembly_data_report.jsonl` (`linkedAssemblies` or `currentAccession`).
    3. If strict historical reproducibility is required, maintain a verified local archive of pre-downloaded database fastas.
    4. Otherwise, flag the species gracefully as `EXCLUDED_NCBI_UNAVAILABLE` in the status ledger.

---

## 3. Implemented Protections in ValidateClassificationsSMK (v2)

To ensure the pipeline never halts or fails unexpectedly due to live NCBI metadata drift:
1. **Live Pre-Flight Query in Stage 0 (`resolve_clades_and_assemblies.py`)**:
   Before running the 4-tier assembly selection cascade, Stage 0 queries NCBI `datasets summary genome accession` in a single bulk request for all candidate accessions in target clades, filtering out suppressed accessions in real time.
2. **Clade Fallback**:
   If the top-ranked assembly in `database_manifest.tsv` has been suppressed, the cascade automatically falls back to the next best live assembly in that clade.
3. **Graceful Status Ledger Tracking**:
   If all assemblies in a clade are suppressed, the taxon is marked `EXCLUDED_NCBI_UNAVAILABLE` with the exact submitter suppression reason, allowing all other species to continue uninterrupted.
