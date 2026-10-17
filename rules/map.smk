from pathlib import Path


# ========================== #
# --- FILE PATH PATTERNS --- #
# ========================== #

# --- output patterns --- #
MAPPING_OUT_DIR = OUTPUT_DIR / "mapping"
MINIMAP2_INDEX_DIR = OUTPUT_DIR / "minimap2_indices"
MINIMAP2_INDEX_PATTERN = str(MINIMAP2_INDEX_DIR / "{genome_basename}.mmi")
RAW_BAM_OUT_PATTERN = str(MAPPING_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapping.bam")
SORTED_BAM_OUT_PATTERN = str(MAPPING_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapping.sorted.bam")
DEDUP_BAM_OUT_PATTERN = str(MAPPING_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapping_dedup.bam")
DEDUP_BAM_INDEX_OUT_PATTERN = str(MAPPING_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapping_dedup.bam.bai")
STATS_OUT_PATTERN = str(MAPPING_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapping_stats.txt")
COVERAGE_OUTPUT_PATTERN = str(MAPPING_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_coverage.tsv")

# --- log patterns --- #
MAPPING_LOG_DIR = LOG_DIR / "mapping"
LOG_MINIMAP2_INDEX_PATTERN = str(MAPPING_LOG_DIR / "{genome_basename}__minimap2_index.log")
LOG_MAP_PATTERN = str(MAPPING_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_minimap2_map.log")
LOG_SORT_PATTERN = str(MAPPING_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_sambamba_sort.log")
LOG_MARKDUP_PATTERN = str(MAPPING_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_sambamba_markdup.log")
LOG_CALC_COV_PATTERN = str(MAPPING_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_coverage.log")

# ==================== #
# --- HELPER FUNCS --- #
# ==================== #

def get_fasta_path(wildcards):
    """Get the path to the reference file, uses the mapping spec DF"""
    basename = wildcards.genome_basename
    matching_rows = MAPPING_SPEC_DF.query('genome_basename == @basename')
    fasta_file_name = matching_rows['fasta_file'].iloc[0]    # Fasta filename
    fasta_file_path = Path(GENOMES_DIR) / fasta_file_name    # Full path
    return str(fasta_file_path)

def get_all_coverage_files(wildcards):
    """Gets list of all per-sample coverage files. For collection rule at the end + checkpoint"""
    return get_target_outputs(MAPPING_SPEC_DF, COVERAGE_OUTPUT_PATTERN)

# ============= #
# --- RULES --- #
# ============= #

localrules:
    WF_collect_mapping_branch

# --- Index Reference Genome using Minimap2 --- #
rule minimap2_index:
    input:
        fasta = get_fasta_path
    output:
        idx = MINIMAP2_INDEX_PATTERN
    log:
        path = LOG_MINIMAP2_INDEX_PATTERN
    conda:
        ".." / ENVS_DIR / "map.yaml"
    threads: config.get("MINIMAP2_INDEX_THREADS", 4)
    resources:
        mem_mb=config.get("MINIMAP2_INDEX_MEM_MB", 8000),
        runtime=config.get("MINIMAP2_INDEX_RUNTIME", "30m"),
        cpus_per_task=config.get("MINIMAP2_INDEX_THREADS", 4)
    benchmark: LOG_MINIMAP2_INDEX_PATTERN.replace("log", "benchmark")
    shell:
        "minimap2 -t {threads} -d {output.idx} {input.fasta} > {log.path} 2>&1"

# --- Map Reads to Reference using Minimap2 --- #
rule minimap2_map_reads:
    input:
        idx = MINIMAP2_INDEX_PATTERN,
        r1 = STARTING_FASTQ_PATTERN.format(tax_id="{tax_id}", read_pair="1"),
        r2 = STARTING_FASTQ_PATTERN.format(tax_id="{tax_id}", read_pair="2")
    output:
        bam = RAW_BAM_OUT_PATTERN,
        stats = STATS_OUT_PATTERN
    log:
        path = LOG_MAP_PATTERN
    conda:
        ".." / ENVS_DIR / "map.yaml"
    threads: config.get("MINIMAP2_MAP_THREADS", 8)
    resources:
        mem_mb=config.get("MINIMAP2_MAP_MEM_MB", 16000),
        runtime=config.get("MINIMAP2_MAP_RUNTIME", "2h"),
        cpus_per_task=config.get("MINIMAP2_MAP_THREADS", 8)
    retries: 3
    benchmark: LOG_MAP_PATTERN.replace("log", "benchmark")
    shell:
        r"""
        (minimap2 -ax sr -t {threads} {input.idx} {input.r1} {input.r2} \
            | samtools view -q 20 -F 0x900 -b - > {output.bam}) 2> {log.path}
        samtools flagstat {output.bam} > {output.stats} 2>> {log.path}
        """

# --- Sort BAM file using Sambamba --- #
rule sambamba_sort:
    input:
        bam = RAW_BAM_OUT_PATTERN
    output:
        bam = SORTED_BAM_OUT_PATTERN
    log:
        path = LOG_SORT_PATTERN
    conda:
        ".." / ENVS_DIR / "map.yaml"
    params:
        mem_limit_gb = config.get("SAMBAMBA_SORT_MEM_LIMIT_GB", "8G")
    threads: config.get("SAMBAMBA_SORT_THREADS", 4)
    resources:
        mem_mb=config.get("SAMBAMBA_SORT_MEM_MB", 10000),
        runtime=config.get("SAMBAMBA_SORT_RUNTIME", "1h"),
        cpus_per_task=config.get("SAMBAMBA_SORT_THREADS", 4)
    benchmark: LOG_SORT_PATTERN.replace("log", "benchmark")
    shell:
        "sambamba sort "
            "-t {threads} "
            "-m {params.mem_limit_gb} "
            "-o {output.bam} "
            "{input.bam} 2> {log.path} "

# --- Mark/Remove Duplicates using Sambamba --- #
rule mark_duplicates_sambamba:
    input:
        # Input is the coordinate-sorted BAM from sambamba_sort
        sorted_bam = SORTED_BAM_OUT_PATTERN
    output:
        # Output is the BAM file with duplicates removed
        dedup_bam = DEDUP_BAM_OUT_PATTERN,
        index = DEDUP_BAM_INDEX_OUT_PATTERN
    params:
        # Define a temporary directory relative to the output
        tmpdir = lambda wildcards, output: Path(output.dedup_bam).parent / "tmp_sambamba_markdup"
    log:
        path = LOG_MARKDUP_PATTERN
    conda:
         ".." / ENVS_DIR / "map.yaml"
    threads: config.get("SAMBAMBA_MARKDUP_THREADS", 4)
    resources:
        mem_mb=config.get("SAMBAMBA_MARKDUP_MEM_MB", 8000),
        runtime=config.get("SAMBAMBA_MARKDUP_RUNTIME", "2h"),
        cpus_per_task=config.get("SAMBAMBA_MARKDUP_THREADS", 4)
    benchmark: LOG_MARKDUP_PATTERN.replace("log", "benchmark")
    shell:
        "sambamba markdup -t {threads} "
            "--tmpdir={params.tmpdir} "
            "--remove-duplicates "
            "{input.sorted_bam} "
            "{output.dedup_bam} "
            "> {log.path} 2>&1 && "
        "rm -rf {params.tmpdir} "

# --- Calculate Real Read Coverage using Samtools --- #
rule calculate_coverage:
    input:
        # Use the final deduplicated BAM
        bam = DEDUP_BAM_OUT_PATTERN
    output:
        # File containing coverage
        cov_file = COVERAGE_OUTPUT_PATTERN
    log:
        path = LOG_CALC_COV_PATTERN
    conda:
         ".." / ENVS_DIR / "map.yaml"
    threads: config.get("CALC_REAL_COV_THREADS", 1)
    resources:
        runtime=config.get("CALC_REAL_COV_RUNTIME", "60m"),
        mem_mb=config.get("CALC_REAL_COV_MEM_MB", 1000),
        cpus_per_task=config.get("CALC_REAL_COV_THREADS", 1)
    benchmark: LOG_CALC_COV_PATTERN.replace("log", "benchmark")
    shell:
        # Use a samtools/awk pipe to calculate coverage
        r"""
        coverage=$(samtools depth -a {input.bam} | awk '{{sum+=$3; count++}} END {{if (count > 0) print sum/count; else print 0}}') 2> {log.path}
        echo -e {wildcards.tax_id}\\t{wildcards.sanitized_name}\\t{wildcards.genome_basename}\\t${{coverage}} > {output.cov_file} 2>> {log.path}
        """

# --- Collect Coverage Results --- #
rule WF_collect_mapping_branch:
    input:
        cov_files = get_all_coverage_files,
        start_ts = WF_TIMESTAMP_START_TARGET
    output:
        WF_MAPPING_BRANCH_TARGET
    shell:
        "touch {output}"

