from pathlib import Path
import json


# ========================== #
# --- FILE PATH PATTERNS --- #
# ========================== #

# --- output patterns --- #
VALIDATION_OUT_DIR = OUTPUT_DIR / "validation"
CHK_AGGREGATE_COVERAGE_OUT = str(VALIDATION_OUT_DIR / "aggregated_coverage.tsv")
SAMTOOLS_STATS_OUT_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_samtools_stats.txt")
SIM_PARAMS_OUT_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulation_params.json")
SIM_READS_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated_{read_pair}.fq.gz")  # Must end in "_{read_pair}.fq.gz"
CLEANED_SIM_READS_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated_headerfixed_R{read_pair}.fq.gz")
AGG_SIM_PATTERN = str(VALIDATION_OUT_DIR / "aggregated_simulated_R{read_pair}.fq.gz")
ALL_KRAKEN_OUT_PATTERN = str(VALIDATION_OUT_DIR / "aggregated_simulated.kraken.out")
ALL_KRAKEN_REPORT_PATTERN = str(VALIDATION_OUT_DIR / "aggregated_simulated.kraken.report")
KRAKEN_OUT_SIM_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated.kraken.out")           # NB: Changing this will affect the split_kraken_output rule
KRAKEN_REPORT_SIM_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated.kraken.report")     # NB: Changing this will affect the split_kraken_output rule
CORRECT_READ_IDS_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated_correctly_classified_read_ids.txt")
FILTERED_SIM_READS_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulated_correctly_classified_R{read_pair}.fq.gz")
MAPPED_SIM_BAM_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapped_simulated_correctly_classified.bam")
MAPPED_SIM_STATS_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapped_simulated_correctly_classified_stats.txt")
SORTED_MAPPED_SIM_BAM_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_mapped_simulated_correctly_classified.sorted.bam")
CLASSIFIABLE_REGIONS_BED_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_classifiable_regions.bed")
OVERLAPPING_READS_BAM_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_real_reads_in_classifiable_regions.bam")
PER_SAMPLE_SUMMARY_PATTERN = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_summary.tsv")
PLOT_SCATTER_HTML_PATTERN  = str(VALIDATION_OUT_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_scatterplot_reads_vs_classifiable_length.html")
VALIDATION_BRANCH_REPORT_PATTERN = str(VALIDATION_OUT_DIR / "joint_validation_report.html")

# --- log patterns --- #
VALIDATION_LOG_DIR = LOG_DIR / "validation"
LOG_CHK_AGGREGATE_COVERAGE = str(VALIDATION_LOG_DIR / "aggregate_coverage.log")
LOG_ESTIMATE_PARAMS_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_estimate_simulation_params.log")
LOG_SIMULATE_READS_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_simulate_reads_pirs.log")
LOG_REWRITE_HEADERS_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_R{read_pair}__rewrite_fastq_headers.log")
LOG_KRAKEN2_BATCH_PATTERN = str(VALIDATION_LOG_DIR / "classify_all_simulated.log")
LOG_SPLIT_KRAKEN_PATTERN = str(VALIDATION_LOG_DIR / "split_kraken_output.log")
LOG_EXTRACT_IDS_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_extract_correct_read_ids.log")
LOG_FILTER_FASTQ_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_R{read_pair}__filter_fastq_with_seqtk.log")
LOG_MAP_SIM_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_map_correct_simulated.log")
LOG_SORT_SIM_PATTERN = str(VALIDATION_LOG_DIR / "{sanitized_name}/{sanitized_name}_taxID-{tax_id}_genome-{genome_basename}_sort_mapped_simulated.log")
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

def get_simulation_validation_targets(wildcards, pattern_template, **fixed_wildcards):
    """
    Input function generalized for simulation_validation.smk. Waits for aggregate_coverage_results
    checkpoint, reads the list of passing samples, and expands the final target patterns only for those samples.
    Allows fixing specific wildcards (like read_pair).
    """
    final_validation_targets = []
    
    # Step 1: Wait for the coverage aggregation checkpoint to finish
    passed_samples_tsv = checkpoints.aggregate_coverage_results.get().output[0]

    # Step 2: Read the TSV file into a pandas DataFrame
    passed_df = pd.read_csv(passed_samples_tsv, sep='\t', dtype=str)
    
    # Step 3: Generate final targets using expand() only for passing samples
    if not passed_df.empty:
        
        # Extract lists needed for expand's zip mode
        tax_id_list = passed_df["TaxID"].tolist()
        sanitized_name_list = passed_df["SanitizedName"].tolist()
        genome_basename_list = passed_df["GenomeBasename"].tolist()

        # Prepare expansion keyword arguments
        expand_kwargs = {
            'tax_id': tax_id_list,
            'sanitized_name': sanitized_name_list,
            'genome_basename': genome_basename_list}

        # Add any fixed wildcards provided
        for wc_name, wc_value in fixed_wildcards.items():
            # Ensure the value is a list of the same length for zip mode
            expand_kwargs[wc_name] = [wc_value] * len(tax_id_list)
        
        # Generate targets using expand() with zip
        target_files = expand(pattern_template, zip, **expand_kwargs)
        final_validation_targets.extend(target_files)

    else:
        logger.warning(f"Passed samples file '{passed_samples_tsv}' is empty. No validation targets generated.")

    # Return the list of validation targets (may be empty)
    return final_validation_targets

# --- Input function for the aggregate_simulated_fastqs rule --- #
def get_rewritten_simulated_fastqs(wildcards):
    """
    Returns the list of R1 or R2 rewritten FASTQ files based on the
    wildcards.read_pair value ('1' or '2').
    """
    # Call the core helper function directly, passing:
    return get_simulation_validation_targets(
        wildcards,
        CLEANED_SIM_READS_PATTERN,
        read_pair=wildcards.read_pair)

def get_all_kraken_outs(wildcards):
    """Gets list of all expected per-sample kraken output files."""
    return get_target_outputs(MAPPING_SPEC_DF, KRAKEN_OUT_SIM_PATTERN)

def get_all_html_files(wildcards):
    """Gets list of all per-sample HTML files."""
    return get_simulation_validation_targets(
        wildcards,
        PLOT_SCATTER_HTML_PATTERN)


# =============== #
# --- OUTPUTS --- #
# =============== #

# TODO: The split kraken output rule should only split into the files that correspond 
# to the passed tax_id/genome_basename combinations. Currently attempts to split all
# starting combinations, but it apparently doen't cause a crash - just some empty output folders.

# --- Generate the list of expected Kraken output files (after splitting them) ---
# This is needed because we can't use a function in the output section of the rule
# But since we know what the output files will be called, we can generate a list of them
# here and use it in the rule.
ALL_EXPECTED_KRAKEN_OUTS = get_all_kraken_outs(None)
if not ALL_EXPECTED_KRAKEN_OUTS:
    logger.warning("No expected Kraken output files generated based on mapping spec. Check patterns/spec file.")


# ============= #
# --- RULES --- #
# ============= #

localrules:
    aggregate_coverage_results,
    WF_collect_validation_branch

# --- Checkpoint: Aggregate Coverage Results and Filter Samples --- #
checkpoint aggregate_coverage_results:
    input:
        # Use helper function to get list of all per-sample coverage files
        coverage_files = get_all_coverage_files
    output:
        # Single TSV listing samples that passed the threshold
        passed_samples_tsv = CHK_AGGREGATE_COVERAGE_OUT
    params:
        cov_threshold = config.get("MIN_MEAN_DEPTH_COVERAGE_FOR_VALIDATION", 0.1)
    log:
        path = LOG_CHK_AGGREGATE_COVERAGE
    conda:
        ".." / ENVS_DIR / "analysis.yaml"
    benchmark: LOG_CHK_AGGREGATE_COVERAGE.replace("log", "benchmark")
    script:
        "../scripts/aggregate_and_filter_coverage.py"

# --- Estimate Simulation Parameters using samtools stats --- #
rule estimate_simulation_params:
    input:
        # The deduplicated BAM file from the mapping step
        bam = DEDUP_BAM_OUT_PATTERN
    output:
        # Raw output from samtools stats and JSON file with simulation parameters
        stats_raw = SAMTOOLS_STATS_OUT_PATTERN,
        params_json = SIM_PARAMS_OUT_PATTERN
    log:
        path = LOG_ESTIMATE_PARAMS_PATTERN
    conda:
        ".." / ENVS_DIR / "map.yaml"
    threads: config.get("SAMTOOLS_STATS_THREADS", 1)
    resources:
        mem_mb=config.get("SAMTOOLS_STATS_MEM_MB", 1000),
        runtime=config.get("SAMTOOLS_STATS_RUNTIME", "15m"),
        cpus_per_task=config.get("SAMTOOLS_STATS_THREADS", 1)
    benchmark: LOG_ESTIMATE_PARAMS_PATTERN.replace("log", "benchmark")
    script:
        "../scripts/run_samtools_stats_parse.py"

# --- Simulate error-free reads from genome using pIRS --- #
rule simulate_reads_pirs:
    input:
        # Input is the reference genome and the simulation params
        ref_fasta = get_fasta_path,
        params_json = SIM_PARAMS_OUT_PATTERN
    output:
        # Output simulated R1 and R2 reads
        r1 = SIM_READS_PATTERN.format(
            sanitized_name="{sanitized_name}", tax_id="{tax_id}", genome_basename="{genome_basename}", read_pair="1"),
        r2 = SIM_READS_PATTERN.format(
            sanitized_name="{sanitized_name}", tax_id="{tax_id}", genome_basename="{genome_basename}", read_pair="2")
    log:
        path = LOG_SIMULATE_READS_PATTERN
    conda:
        ".." / ENVS_DIR / "simulate_reads.yaml"
    threads: config.get("PIRS_SIM_THREADS", 4)
    resources:
        runtime=config.get("PIRS_SIM_RUNTIME", "2h"),
        mem_mb=config.get("PIRS_SIM_MEM_MB", 4000),
        cpus_per_task=config.get("PIRS_SIM_THREADS", 4)
    retries: 8
    benchmark: LOG_SIMULATE_READS_PATTERN.replace("log", "benchmark")
    script:
        "../scripts/run_pirs_simulate_reads.py"

# --- Rewrite FASTQ Headers to Embed Source Info (using awk) --- #
rule rewrite_fastq_headers:
    input:
        # Input are the raw simulated reads from pIRS
        raw_fq = SIM_READS_PATTERN
    output:
        # Output are cleaned FASTQ files with tax_id and genome_basename embedded in the headers
        cleaned_fq = CLEANED_SIM_READS_PATTERN
    log:
        path = LOG_REWRITE_HEADERS_PATTERN
    params:
        pigz_threads = 2
    conda:
        ".." / ENVS_DIR / "extract_reads.yaml"
    threads: config.get("REWRITE_HEADERS_THREADS", 5)
    resources:
        runtime=config.get("REWRITE_HEADERS_RUNTIME", "2h"),
        mem_mb=config.get("REWRITE_HEADERS_MEM_MB", 2000),
        cpus_per_task=config.get("REWRITE_HEADERS_THREADS", 5)
    shell:
        r"""
        AWK_CMD='NR % 4 == 1 {{ sub(/\/[12]$/, ""); print $0 "_taxid=" tid "_genome=" gbase; next }} {{ print }}' && \
        pigz -p {params.pigz_threads} -dc {input.raw_fq} \
            | awk -v tid={wildcards.tax_id} -v gbase={wildcards.genome_basename} "$AWK_CMD" \
            | pigz -p {params.pigz_threads} -c > {output.cleaned_fq} 2> {log.path}
        """

# --- Aggregate Rewritten Simulated FASTQs --- #
rule aggregate_simulated_fastqs:
    input:
        # The rewritten simulated FASTQ files
        fq_files = get_rewritten_simulated_fastqs
    output:
        # Output is a single aggreated FASTQ file
        agg_fq = AGG_SIM_PATTERN
    threads: config.get("AGGREGATE_FASTQ_THREADS", 2)
    resources:
        runtime=config.get("AGGREGATE_FASTQ_RUNTIME", "1h"),
        mem_mb=config.get("AGGREGATE_FASTQ_MEM_MB", 2000),
        cpus_per_task=config.get("AGGREGATE_FASTQ_THREADS", 2)
    shell:
        "cat {input.fq_files} > {output.agg_fq}"

# --- Classify ALL Concatenated Simulated Reads --- #
rule classify_all_concatenated:
    input:
        # Input are the two aggregated FASTQ files
        r1 = AGG_SIM_PATTERN.format(read_pair="1"),
        r2 = AGG_SIM_PATTERN.format(read_pair="2")
    output:
        # Outputs are aggregated kraken output and report files
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
        runtime=config.get("KRAKEN2_RUNTIME", "6h"),
        mem_mb=config.get("KRAKEN2_MEM_MB", 500000),
        cpus_per_task=config.get("KRAKEN2_THREADS", 32)
    benchmark: LOG_KRAKEN2_BATCH_PATTERN.replace("log", "benchmark")
    shell:
        "export OMP_NUM_THREADS={threads} && "
        "module purge && "
        "module load intel/2025a && " 
        "{params.executable} "
            "--db {params.db_path} "
            "--threads {threads} "
            "--paired "
            "--output {output.kraken_out} "
            "--report {output.report} "
            "--gzip-compressed "
            "--confidence {params.confidence} "
            "--minimum-hit-groups {params.mhg} "
            "{input.r1} {input.r2} "
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
        runtime=config.get("SPLIT_KRAKEN_RUNTIME", "60m"),
        mem_mb=config.get("SPLIT_KRAKEN_MEM_MB", 2000),
        cpus_per_task=config.get("SPLIT_KRAKEN_THREADS", 1)
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
        runtime=config.get("EXTRACT_IDS_RUNTIME", "60m"),
        mem_mb=config.get("EXTRACT_IDS_MEM_MB", 2000),
        cpus_per_task=config.get("EXTRACT_IDS_THREADS", 1)
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
        # The rewritten simulated FASTQ files and the read IDs to extract
        fastq = CLEANED_SIM_READS_PATTERN,
        read_ids = CORRECT_READ_IDS_PATTERN
    output:
        # Final filtered FASTQ files (contains simulated reads that were correctly classified)
        fastq = FILTERED_SIM_READS_PATTERN
    log:
        path = LOG_FILTER_FASTQ_PATTERN
    conda:
        ".." / ENVS_DIR / "extract_reads.yaml"
    threads: config.get("SEQTK_THREADS", 1)
    resources:
        runtime=config.get("SEQTK_RUNTIME", "60m"),
        mem_mb=config.get("SEQTK_MEM_MB", 2000),
        cpus_per_task=config.get("SEQTK_THREADS", 1)
    benchmark: LOG_FILTER_FASTQ_PATTERN.replace("log", "benchmark")
    shell:
        "seqtk subseq {input.fastq} <(cut -f 3 {input.read_ids}) | pigz -c > {output.fastq} "
        "2> {log.path}"

# --- Map Correctly Classified Simulated Reads using BBMap --- #
# Run BBmap with the same settings as in the map.smk rule, but with different input/output patterns
# and treating ambiguous reads differently - randomly assign them to one of the best positions
use rule bbmap_map_reads as map_correct_simulated with:
    input:
        # Use the filtered FASTQ files from seqtk step
        idx_dir = BBMAP_INDEX_DIR_PATTERN,
        r1 = FILTERED_SIM_READS_PATTERN.format(
            sanitized_name="{sanitized_name}", tax_id="{tax_id}", genome_basename="{genome_basename}", read_pair="1"),
        r2 = FILTERED_SIM_READS_PATTERN.format(
            sanitized_name="{sanitized_name}", tax_id="{tax_id}", genome_basename="{genome_basename}", read_pair="2")
    output:
        # Output BAM and mapping stats
        bam = MAPPED_SIM_BAM_PATTERN,
        stats = MAPPED_SIM_STATS_PATTERN
    params:
        ambig = 'random'
    log:
        path = LOG_MAP_SIM_PATTERN
    benchmark: LOG_MAP_SIM_PATTERN.replace("log", "benchmark")

# --- Sort Mapped Simulated Reads BAM file by coordinate --- #
# Reuse the sambamba sort rule from map.smk, but with different input/output patterns
use rule sambamba_sort as sort_mapped_simulated with:
    input:
        bam = MAPPED_SIM_BAM_PATTERN
    output:
        bam = SORTED_MAPPED_SIM_BAM_PATTERN
    log:
        path = LOG_SORT_SIM_PATTERN
    benchmark: LOG_SORT_SIM_PATTERN.replace("log", "benchmark")

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
        sim_params_jsons = expand(SIM_PARAMS_OUT_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
        pirs_logs = expand(LOG_SIMULATE_READS_PATTERN, zip, sanitized_name=ALL_SANITIZED_NAMES, tax_id=ALL_TAX_IDS, genome_basename=ALL_GENOME_BASENAMES),
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
             "sim_cov_target": config.get("COVERAGE"),
             "bbmap_minid": config.get("BBMAP_MINID"),
             "fastq_pattern": config.get("FASTQ_FILE_PATTERN"),
             "mapping_spec_path": str(MAPPING_SPEC)},
        # Patterns are used to find the input files that are related to each
        # combination of tax_id and genome_basename (using MAPPING_SPEC_DF)
        path_patterns = {
             "real_stats": STATS_OUT_PATTERN,
             "markdup_logs": LOG_MARKDUP_PATTERN,
             "dedup_bams": DEDUP_BAM_OUT_PATTERN,
             "coverage_files": COVERAGE_OUTPUT_PATTERN,
             "sim_params_jsons": SIM_PARAMS_OUT_PATTERN,
             "pirs_logs": LOG_SIMULATE_READS_PATTERN,
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