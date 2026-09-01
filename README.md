# ValidateClassificationsSMK

Map taxon-specific reads (e.g., previously classified by Kraken 2) to reference genomes using BBMap and validate classifications through simulated read benchmarking and an interactive joint HTML report.

This pipeline is designed for high-performance computing clusters using the SLURM Workload Manager and Conda package manager. Profiles are provided for PDC Dardel (`slurm/pdc`) and HPC2N Kebnekaise (`slurm/hpc2n`).

---

## Workflow Overview

The pipeline consists of two main branches:

### 1. Mapping Branch (`rules/map.smk`)
For each TaxID and reference genome defined in `supporting_files/mapping_specification.tsv`:
1. **Index Reference Genome**: Builds BBMap index files for the reference FASTA genome.
2. **Map Reads**: Maps real paired-end reads (`input/fastq/taxID-{tax_id}_R{1,2}.fq.gz`) to the reference genome using `bbmap.sh` with strict pairing filters (`killbadpairs=t`, `pairedonly=t`).
3. **Sort & Mark Duplicates**: Coordinate-sorts BAMs with `sambamba sort` and flags duplicates with `sambamba markdup`.
4. **Coverage & Stats Calculation**: Calculates mapping coverage across genomes and computes insert size / read length distributions.

### 2. Simulation & Validation Branch (`rules/simulation_validation.smk`)
*(Enabled when `RUN_VALIDATION: True` in `config.yaml`)*
1. **Read Simulation (pIRS)**: Generates synthetic error-free reads matching empirical insert sizes and read lengths.
2. **Kraken 2 Classification**: Concatenates and classifies simulated reads against the reference Kraken 2 database.
3. **Classifiable Regions BED**: Identifies regions with sufficient reference coverage to be classified correctly.
4. **Overlap & Correlation Analysis**: Intersects real read alignments with classifiable regions and measures per-contig correlations.
5. **Joint Validation Report**: Compiles workflow metadata, mapping metrics, and interactive Plotly scatterplots into a consolidated report (`output/validation/joint_validation_report.html`).

---

## Usage

### Prerequisites & Setup

#### 1. Snakemake Environment
Ensure you have Snakemake (v9+) installed with the SLURM executor plugin:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake>=9.1.3 snakemake-executor-plugin-slurm pyyaml
conda activate snakemake
```

#### 2. Reference Genomes
Place or symlink uncompressed reference genomes (`.fasta`, `.fna`, or `.fa`) into:
```
supporting_files/genomes/
```

#### 3. Mapping Specification File (`supporting_files/mapping_specification.tsv`)
A tab-delimited TSV file linking taxonomic IDs to reference genomes and species names:
* Path configurable via `MAPPING_SPECIFICATION` in `config.yaml`.
* Template provided in `supporting_files/mapping_specification_TEMPLATE.tsv`.

**Format (tab-delimited):**
```tsv
tax_id	name	fasta_file
30423	Coloeus monedula	GCA_013407035.1_ASM1340703v1_genomic.fna
37598	Fringilla coelebs	GCA_963513975.1_bFriCoe1.1_genomic.fna
9187	Turdus merula	GCA_046127255.1_ASM4612725v1_genomic.fna
```
* **`tax_id`**: Integer taxonomic identifier.
* **`name`**: Descriptive taxon name (alphanumerics, spaces, underscores, and hyphens allowed).
* **`fasta_file`**: Filename of the uncompressed reference FASTA in `supporting_files/genomes/`.

#### 4. Input FASTQ Files
Place or symlink paired-end FASTQ files in `input/fastq/`.
* File naming pattern must match `FASTQ_FILE_PATTERN` in `config.yaml` (default: `fastq/taxID-{tax_id}_R{read_pair}.fq.gz`, resolved relative to `input/`).
* Example: `input/fastq/taxID-30423_R1.fq.gz` and `input/fastq/taxID-30423_R2.fq.gz`.

#### 5. Kraken 2 Reference Database
Set `KRAKEN2_DB_PATH` in `config.yaml` to point to a Kraken 2 database directory containing `taxonomy/names.dmp` and `taxonomy/nodes.dmp`.

---

### Running the Pipeline

1. **Activate Environment**:
   ```bash
   conda activate snakemake
   ```

2. **Dry Run (Verify Workflow Plan)**:
   ```bash
   snakemake -n --use-conda
   ```

3. **Execute via SLURM Profile**:
   * For PDC Dardel:
     ```bash
     snakemake --profile slurm/pdc --use-conda
     ```
   * For HPC2N Kebnekaise:
     ```bash
     snakemake --profile slurm/hpc2n --use-conda
     ```

4. **Execute Locally / Interactive Node**:
   ```bash
   snakemake --use-conda --cores <NUM_CORES>
   ```

---

### Output

Key output files produced in `output/`:
* `output/mapping/`: Indexed, deduplicated BAM files and coverage/stats reports.
* `output/validation/`: Simulated reads, classifiable BED regions, and correlation TSVs.
* **`output/validation/joint_validation_report.html`**: Interactive HTML report summarizing all runs, metrics, and correlation plots across all taxa.
