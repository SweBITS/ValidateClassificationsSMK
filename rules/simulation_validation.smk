from pathlib import Path
import json


# ========================== #
# --- FILE PATH PATTERNS --- #
# ========================== #

# --- output patterns --- #
VALIDATION_OUT_DIR = OUTPUT_DIR / "validation"
SLICED_SIM_READS_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_sliced_reads.fa.gz")
AGG_SIM_PATTERN = str(VALIDATION_OUT_DIR / "aggregated_simulated.fa.gz")
ALL_KRAKEN_OUT_PATTERN = str(VALIDATION_OUT_DIR / "aggregated_simulated.kraken.out")
ALL_KRAKEN_REPORT_PATTERN = str(VALIDATION_OUT_DIR / "aggregated_simulated.kraken.report")
KRAKEN_OUT_SIM_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated.kraken.out")           # NB: Changing this will affect the split_kraken_output rule
KRAKEN_REPORT_SIM_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated.kraken.report")     # NB: Changing this will affect the split_kraken_output rule
CORRECT_READ_IDS_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated_correctly_classified_read_ids.txt")
FILTERED_SIM_READS_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated_correctly_classified.fa.gz")
MAPPED_SIM_STATS_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapped_simulated_correctly_classified_stats.txt")
SORTED_MAPPED_SIM_BAM_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapped_simulated_correctly_classified.sorted.bam")
CLASSIFIABLE_REGIONS_BED_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_classifiable_regions.bed")
OVERLAPPING_READS_BAM_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_real_reads_in_classifiable_regions.bam")
PER_SAMPLE_SUMMARY_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_summary.tsv")
PLOT_SCATTER_HTML_PATTERN  = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_scatterplot_reads_vs_classifiable_length.html")
VALIDATION_BRANCH_REPORT_PATTERN = str(VALIDATION_OUT_DIR / "joint_validation_report.html")

# --- log patterns --- #
VALIDATION_LOG_DIR = LOG_DIR / "validation"
LOG_SLICE_READS_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_slice_reference_genome.log")
LOG_KRAKEN2_BATCH_PATTERN = str(VALIDATION_LOG_DIR / "classify_all_simulated.log")
LOG_SPLIT_KRAKEN_PATTERN = str(VALIDATION_LOG_DIR / "split_kraken_output.log")
LOG_EXTRACT_IDS_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_extract_correct_read_ids.log")
LOG_FILTER_FASTQ_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_filter_fastq_with_seqtk.log")
LOG_MAP_SIM_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_minimap2_map_simulated.log")
LOG_DEF_CR_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_define_classifiable_regions.log")
LOG_INTERSECT_READS_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_intersect_real_reads_with_classifiable_regions.log")
LOG_SUMMARIZE_VALIDATION_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_summarize_validation_coverage.log")
LOG_PLOT_SCATTER_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_plot_validation_scatter.log")
LOG_VALIDATION_BRANCH_REPORT_PATTERN = str(VALIDATION_LOG_DIR / "joint_validation_report.log")

# ================================== #
# --- HELPER FUNCS AND VARIABLES --- #
# ================================== #

# Linearized values of all tax IDs, sanitized names, and genome basenames
# to be used in the expand functions of rule validation_branch_report
ALL_TAX_IDS = MAPPING_SPEC_DF['tax_id'].tolist()
ALL_SANITIZED_NAMES = MAPPING_SPEC_DF['sanitized_name'].tolist()
ALL_GENOME_BASENAMES = MAPPING_SPEC_DF['genome_basename'].tolist()

# --- Target Expansion Helper Functions --- #
def get_sliced_simulated_reads(wildcards):
    """
    Returns the list of sliced FASTA files for all target species.
    """
    return expand(
        SLICED_SIM_READS_PATTERN,
        zip,
        tax_id=ALL_TAX_IDS,
        sanitized_name=ALL_SANITIZED_NAMES,
        genome_basename=ALL_GENOME_BASENAMES)

def get_all_kraken_outs(wildcards):
    """Gets list of all expected per-sample kraken output files."""
    return get_target_outputs(MAPPING_SPEC_DF, KRAKEN_OUT_SIM_PATTERN)

def get_all_html_files(wildcards):
    """Gets list of all per-sample HTML files."""
    return expand(
        PLOT_SCATTER_HTML_PATTERN,
        zip,
        tax_id=ALL_TAX_IDS,
        sanitized_name=ALL_SANITIZED_NAMES,
        genome_basename=ALL_GENOME_BASENAMES)


# =============== #
# --- OUTPUTS --- #
# =============== #

# --- Generate the list of expected Kraken output files (after splitting them) ---
ALL_EXPECTED_KRAKEN_OUTS = get_all_kraken_outs(None)
if not ALL_EXPECTED_KRAKEN_OUTS:
    logger.warning("No expected Kraken output files generated based on mapping spec. Check patterns/spec file.")


# ============= #
# --- RULES --- #
# ============= #

localrules:
    WF_collect_validation_branch


# --- Slice Reference Genome into overlapping pseudo-reads using seqkit sliding --- #
rule slice_reference_genome:
    input:
        fasta = get_fasta_path
    output:
        reads = SLICED_SIM_READS_PATTERN
    params:
        window = config.get("SLICE_WINDOW", 300),
        step = config.get("SLICE_STEP", 50)
    log:
        path = LOG_SLICE_READS_PATTERN
    conda:
        ".." / ENVS_DIR / "simulate_reads.yaml"
    threads: config.get("SLICE_THREADS", 4)
    resources:
        runtime = config.get("SLICE_RUNTIME", "30m"),
        mem_mb = config.get("SLICE_MEM_MB", 4000),
        cpus_per_task = config.get("SLICE_THREADS", 4)
    benchmark: LOG_SLICE_READS_PATTERN.replace("log", "benchmark")
    shell:
        r"""
        (seqkit sliding -W {params.window} -s {params.step} -g -j {threads} {input.fasta} \
            | awk -v tid={wildcards.tax_id} -v gbase={wildcards.genome_basename} '/^>/ {{print $0 "_taxid=" tid "_genome=" gbase; next}} {{print}}' \
            | pigz -p {threads} -c > {output.reads}) 2> {log.path}
        """

# --- Aggregate Sliced Simulated FASTAs --- #
rule aggregate_simulated_fastqs:
    input:
        fq_files = get_sliced_simulated_reads
    output:
        agg_fq = AGG_SIM_PATTERN
    threads: config.get("AGGREGATE_FASTQ_THREADS", 2)
    resources:
        runtime = config.get("AGGREGATE_FASTQ_RUNTIME", "1h"),
        mem_mb = config.get("AGGREGATE_FASTQ_MEM_MB", 2000),
        cpus_per_task = config.get("AGGREGATE_FASTQ_THREADS", 2)
    shell:
        "cat {input.fq_files} > {output.agg_fq}"

# --- Classify ALL Concatenated Simulated Reads --- #
rule classify_all_concatenated:
    input:
        reads = AGG_SIM_PATTERN
    output:
        kraken_out = ALL_KRAKEN_OUT_PATTERN,
        report = ALL_KRAKEN_REPORT_PATTERN
    params:
        db_path = config["KRAKEN2_DB_PATH"],
        # Determine executable path: use specified path or default 'kraken2' command
        # Path already validated in main Snakefile, if no path is given we assume conda execution
        executable = str(Path(config["KRAKEN2_BIN_DIR"]) / "kraken2") if config.get("KRAKEN2_BIN_DIR") else "kraken2",
        confidence = config.get("KRAKEN2_CONFIDENCE", 0.0),
        mhg = config.get("KRAKEN2_MINIMUM_NUM_HIT_GROUPS", 2)
    log:
        path = LOG_KRAKEN2_BATCH_PATTERN
    conda:
        # Conditionally activate based on executable path
        ".." / ENVS_DIR / "classify.yaml" if not config.get("KRAKEN2_BIN_DIR") else None
    threads: config.get("KRAKEN2_THREADS", 32)
    resources:
        runtime = config.get("KRAKEN2_RUNTIME", "6h"),
        mem_mb = config.get("KRAKEN2_MEM_MB", 500000),
        cpus_per_task = config.get("KRAKEN2_THREADS", 32)
    benchmark: LOG_KRAKEN2_BATCH_PATTERN.replace("log", "benchmark")
    shell:
        "{params.executable} "
            "--db {params.db_path} "
            "--threads {threads} "
            "--output {output.kraken_out} "
            "--report {output.report} "
            "--gzip-compressed "
            "--confidence {params.confidence} "
            "--minimum-hit-groups {params.mhg} "
            "{input.reads} "
            "2> {log.path}"

# --- Split Combined Kraken Output by Sample --- #
rule split_kraken_output:
    input:
        # Input is the single combined kraken output file
        combined_kraken_out = ALL_KRAKEN_OUT_PATTERN
    output:
        # All the kraken read-by-read classification output files,
        # split from the combined output file in classify_all_concatenated
        ALL_EXPECTED_KRAKEN_OUTS
    params:
        # Regex to extract info embedded in read ID by rewrite_fastq_headers rule
        # Assumes header format like @BaseID_taxid=XXX_genome=YYY
        # Capture group 1: TaxID digits, Capture group 2: GenomeBaseName chars
        header_regex = r"_taxid=(\d+)_genome=([a-zA-Z0-9_.\-]+)",

        # Regex to extract wildcards from expected output paths
        # Assumes pattern ".../{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated.kraken.out"
        output_filename_regex = r".*/([a-zA-Z0-9_-]+)_taxID-(\d+)_genome-([a-zA-Z0-9_.\-]+)_simulated\.kraken\.out$"
    log:
        path = LOG_SPLIT_KRAKEN_PATTERN
    conda:
        ".." / ENVS_DIR / "analysis.yaml"
    threads: config.get("SPLIT_KRAKEN_THREADS", 1)
    resources:
        runtime = config.get("SPLIT_KRAKEN_RUNTIME", "60m"),
        mem_mb = config.get("SPLIT_KRAKEN_MEM_MB", 2000),
        cpus_per_task = config.get("SPLIT_KRAKEN_THREADS", 1)
    benchmark: LOG_SPLIT_KRAKEN_PATTERN.replace("log", "benchmark")
    script:
        "../scripts/split_kraken.py"

# --- Extract Correctly Classified Read IDs using selmeout.py --- #
rule extract_correct_read_ids:
    input:
        # Kraken output from simulated reads classification
        kraken_out = KRAKEN_OUT_SIM_PATTERN
    output:
        # Output is a table of read IDs (classified_taxid, root_taxid, read_id)
        read_ids = CORRECT_READ_IDS_PATTERN
    params:
        selmeout_path = "scripts/selmeout.py",
        mode = "clade",
        names = str(Path(config["KRAKEN2_DB_PATH"]) / "taxonomy/names.dmp"),
        nodes = str(Path(config["KRAKEN2_DB_PATH"]) / "taxonomy/nodes.dmp")
    log:
        path = LOG_EXTRACT_IDS_PATTERN
    conda:
        ".." / ENVS_DIR / "extract_reads.yaml"
    threads: config.get("EXTRACT_IDS_THREADS", 1)
    resources:
        runtime = config.get("EXTRACT_IDS_RUNTIME", "60m"),
        mem_mb = config.get("EXTRACT_IDS_MEM_MB", 2000),
        cpus_per_task = config.get("EXTRACT_IDS_THREADS", 1)
    benchmark: LOG_EXTRACT_IDS_PATTERN.replace("log", "benchmark")
    shell:
        "python {params.selmeout_path} "
            "--input {input.kraken_out} "
            "--output {output.read_ids} "
            "--tax_id {wildcards.tax_id} "
            "--mode {params.mode} "
            "--nodes {params.nodes} "
            "--names {params.names} "
            "> {log.path} 2>&1"

# --- Filter correctly classified simulated reads from FASTQ files using seqtk --- #
rule filter_fastq_with_seqtk:
    input:
        reads = SLICED_SIM_READS_PATTERN,
        read_ids = CORRECT_READ_IDS_PATTERN
    output:
        reads = FILTERED_SIM_READS_PATTERN
    log:
        path = LOG_FILTER_FASTQ_PATTERN
    conda:
        ".." / ENVS_DIR / "extract_reads.yaml"
    threads: config.get("SEQTK_THREADS", 1)
    resources:
        runtime = config.get("SEQTK_RUNTIME", "60m"),
        mem_mb = config.get("SEQTK_MEM_MB", 8000),
        cpus_per_task = config.get("SEQTK_THREADS", 1)
    benchmark: LOG_FILTER_FASTQ_PATTERN.replace("log", "benchmark")
    shell:
        "seqtk subseq {input.reads} <(cut -f 3 {input.read_ids}) | pigz -c > {output.reads} 2> {log.path}"

# --- Map Correctly Classified Simulated Reads using Minimap2 --- #
rule map_correct_simulated:
    input:
        idx = MINIMAP2_INDEX_PATTERN,
        reads = FILTERED_SIM_READS_PATTERN
    output:
        bam = SORTED_MAPPED_SIM_BAM_PATTERN,
        stats = MAPPED_SIM_STATS_PATTERN
    log:
        path = LOG_MAP_SIM_PATTERN
    conda:
        ".." / ENVS_DIR / "map.yaml"
    threads: config.get("MINIMAP2_MAP_SIM_THREADS", 8)
    resources:
        mem_mb = config.get("MINIMAP2_MAP_SIM_MEM_MB", 16000),
        runtime = config.get("MINIMAP2_MAP_SIM_RUNTIME", "1h"),
        cpus_per_task = config.get("MINIMAP2_MAP_SIM_THREADS", 8)
    benchmark: LOG_MAP_SIM_PATTERN.replace("log", "benchmark")
    shell:
        r"""
        (minimap2 -ax sr -t {threads} {input.idx} {input.reads} \
            | samtools view -F 0x900 -b - \
            | samtools sort -@ {threads} -o {output.bam} -) 2> {log.path}
        samtools flagstat {output.bam} > {output.stats} 2>> {log.path}
        """


# --- Get Classifiable Regions--- #
rule define_classifiable_regions:
    input:
        # The sorted BAM file containing mapped simulated reads
        bam = SORTED_MAPPED_SIM_BAM_PATTERN
    output:
        # BED file of regions covered >= threshold (marked temp)
        covered_bed = CLASSIFIABLE_REGIONS_BED_PATTERN
    params:
        threshold = 1  # Coverage threshold (site is kept if coverage >= threshold)
    log:
        path = LOG_DEF_CR_PATTERN
    conda:
        ".." / ENVS_DIR / "analysis.yaml"
    threads: config.get("DEF_CR_THREADS", 1)
    resources:
        cpus_per_task=config.get("DEF_CR_THREADS", 1),
        mem_mb=config.get("DEF_CR_MEM_MB", 2000),
        runtime=config.get("DEF_CR_RUNTIME", "30m")
    benchmark: LOG_DEF_CR_PATTERN.replace("log", "benchmark")
    shell:
        # Use bedtools genomecov -bga, coverage filter with awk, merge resulting bed-regions
        "bedtools genomecov -bga -ibam {input.bam} "
            "| awk -v threshold={params.threshold} '$NF >= threshold' "
            "| bedtools merge -i stdin "
            "> {output.covered_bed} "
            "2> {log.path}"

# --- Intersect Real Reads with Classifiable Regions --- #
rule intersect_real_reads_with_classifiable_regions:
    input:
        # Use the deduplicated BAM of REAL reads from map.smk, and the classifiable regions BED
        dedup_bam = DEDUP_BAM_OUT_PATTERN,
        classifiable_regions = CLASSIFIABLE_REGIONS_BED_PATTERN
    output:
        # Output BAM containing unique reads overlapping classifiable regions
        overlapping_bam = OVERLAPPING_READS_BAM_PATTERN
    log:
        path = LOG_INTERSECT_READS_PATTERN
    conda:
        ".." / ENVS_DIR / "analysis.yaml"
    threads: config.get("INTERSECT_READS_THREADS", 1)
    resources:
        runtime=config.get("INTERSECT_READS_RUNTIME", "30m"),
        mem_mb=config.get("INTERSECT_READS_MEM_MB", 2000),
        cpus_per_task=config.get("INTERSECT_READS_THREADS", 1)
    benchmark: LOG_INTERSECT_READS_PATTERN.replace("log", "benchmark")
    shell:
        # Use bedtools intersect -u; output is BAM
        "bedtools intersect "
            "-a {input.dedup_bam} "             # Input BAM (real reads)
            "-b {input.classifiable_regions} "  # Classifiable regions BED
            "-u "                               # Report each alignment record in -a once if *any* overlap exists
            "-wa "                              # Write the original BAM records
            "> {output.overlapping_bam} "
            "2> {log.path}"

# --- Summarize Validation Coverage Per Contig --- #
rule summarize_validation_coverage:
    input:
        # BAM with reads in classifiable regions, and the regions BED to calculate length
        overlapping_bam = OVERLAPPING_READS_BAM_PATTERN,
        classifiable_bed = CLASSIFIABLE_REGIONS_BED_PATTERN
    output:
        summary_tsv = PER_SAMPLE_SUMMARY_PATTERN
    log:
        path = LOG_SUMMARIZE_VALIDATION_PATTERN
    conda:
        ".." / ENVS_DIR / "analysis.yaml"
    threads: config.get("SUMMARIZE_VALIDATION_THREADS", 1)
    resources:
        runtime=config.get("SUMMARIZE_VALIDATION_RUNTIME", "20m"),
        mem_mb=config.get("SUMMARIZE_VALIDATION_MEM_MB", 2000),
        cpus_per_task=config.get("SUMMARIZE_VALIDATION_THREADS", 1)
    benchmark: LOG_SUMMARIZE_VALIDATION_PATTERN.replace("log", "benchmark")
    script:
        "../scripts/summarize_validation.py"

# --- Create Validation Scatter Plot --- #
rule plot_validation_scatter:
    input:
        summary_tsv = PER_SAMPLE_SUMMARY_PATTERN
    output:
        plot_html = PLOT_SCATTER_HTML_PATTERN
    log:
        path = LOG_PLOT_SCATTER_PATTERN
    conda:
        ".." / ENVS_DIR / "analysis.yaml"
    threads: config.get("PLOT_THREADS", 1)
    resources:
        runtime=config.get("PLOT_RUNTIME", "15m"),
        mem_mb=config.get("PLOT_MEM_MB", 1000),
        cpus_per_task=config.get("PLOT_THREADS", 1)
    benchmark: LOG_PLOT_SCATTER_PATTERN.replace("log", "benchmark")
    script:
        "../scripts/plot_scatter_plotly.py"

# --- Final Report Generation --- #
rule validation_branch_report:
    input:
        # General Workflow Info Files
        ts_start = WF_TIMESTAMP_START_TARGET,
        snk_version = WF_SNAKEMAKE_VERSION_TARGET,

        # Per-Validation Files
        real_stats = expand(STATS_OUT_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
        markdup_logs = expand(LOG_MARKDUP_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
        dedup_bams = expand(DEDUP_BAM_OUT_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
        coverage_files = expand(COVERAGE_OUTPUT_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
        slice_logs = expand(LOG_SLICE_READS_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
        correct_ids_logs = expand(LOG_EXTRACT_IDS_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
        sim_stats = expand(MAPPED_SIM_STATS_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
        classifiable_beds = expand(CLASSIFIABLE_REGIONS_BED_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
        summary_tsvs = expand(PER_SAMPLE_SUMMARY_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES)
    output:
        report_html = VALIDATION_BRANCH_REPORT_PATTERN
    log:
        path = LOG_VALIDATION_BRANCH_REPORT_PATTERN
    params:
        # Pass the whole mapping spec DF and relevant config values the script 
        # needs for the general section
        mapping_spec_df = MAPPING_SPEC_DF,
        report_config = lambda w, input: {
             "kraken_db": config.get("KRAKEN2_DB_PATH"),
             "kraken_conf": config.get("KRAKEN2_CONFIDENCE"),
             "kraken_mhg": config.get("KRAKEN2_MINIMUM_NUM_HIT_GROUPS"),
             "cov_cutoff": config.get("MIN_MEAN_DEPTH_COVERAGE_FOR_VALIDATION"),
             "slice_window": config.get("SLICE_WINDOW", 300),
             "slice_step": config.get("SLICE_STEP", 50),
             "fastq_pattern": config.get("FASTQ_FILE_PATTERN"),
             "mapping_spec_path": str(MAPPING_SPEC)},
        # Patterns are used to find the input files that are related to each
        # combination of tax_id and genome_basename (using MAPPING_SPEC_DF)
        path_patterns = {
             "real_stats": STATS_OUT_PATTERN,
             "markdup_logs": LOG_MARKDUP_PATTERN,
             "dedup_bams": DEDUP_BAM_OUT_PATTERN,
             "coverage_files": COVERAGE_OUTPUT_PATTERN,
             "slice_logs": LOG_SLICE_READS_PATTERN,
             "correct_ids_logs": LOG_EXTRACT_IDS_PATTERN,
             "sim_stats": MAPPED_SIM_STATS_PATTERN,
             "classifiable_beds": CLASSIFIABLE_REGIONS_BED_PATTERN,
             "summary_tsvs": PER_SAMPLE_SUMMARY_PATTERN}
    conda:
        ".." / ENVS_DIR / "analysis.yaml"
    threads: config.get("VALIDATION_REPORT_THREADS", 1)
    resources:
        runtime=config.get("VALIDATION_REPORT_RUNTIME", "15m"),
        mem_mb=config.get("VALIDATION_REPORT_MEM_MB", 2000),
        cpus_per_task=config.get("VALIDATION_REPORT_THREADS", 1)
    benchmark: LOG_VALIDATION_BRANCH_REPORT_PATTERN.replace("log", "benchmark")
    script:
        "../scripts/generate_validation_branch_report.py"

# --- Collect Simulation Validation Results --- #
rule WF_collect_validation_branch:
    input:
        html_plots = get_all_html_files,
        joint_report = VALIDATION_BRANCH_REPORT_PATTERN,
        start_ts = WF_TIMESTAMP_START_TARGET
    output:
        WF_VALIDATION_BRANCH_TARGET
    shell:
        "touch {output}"