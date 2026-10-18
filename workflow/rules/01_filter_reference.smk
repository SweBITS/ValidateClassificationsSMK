# ============================================================================ #
# Stage 1: Upfront Reference Contig Filtering (>= 1,000 bp) & Checkpoint Audit
# ============================================================================ #

checkpoint filter_and_audit_assemblies:
    """
    WHY THIS RULE IS A SNAKEMAKE CHECKPOINT:
    Filters all rehydrated reference assemblies to remove contigs < 1,000 bp (unresolved
    repeats, assembly graph tips, or trace contamination), writes the mandatory audit table
    `assembly_filtering_summary.tsv`, and evaluates two Stage 1 quality gates:
      - `kept_assembly_bp < min_kept_assembly_bp` -> `EXCLUDED_INSUFFICIENT_CONTIG_LENGTH`
      - `discarded_bp_pct > max_discarded_bp_pct` -> `EXCLUDED_FILTERED_ASSEMBLY_PROPORTION`

    Emitting `active_species.tsv` as a Snakemake checkpoint ensures that only species passing
    BOTH Stage 0 and Stage 1 trigger downstream read extraction, Kraken 2 streaming, CR binning,
    and Minimap2 alignment, while excluded species proceed directly to the final master report.
    """
    input:
        stage0_ledger = "results/00_ledger/stage0_preflight_and_assemblies.tsv",
        acc_manifest = "results/01_inputs/ncbi_datasets/accession_to_fasta.tsv",
        rehydrate_flag = "results/01_inputs/ncbi_datasets/rehydrate.done"
    output:
        audit_tsv = "results/00_ledger/assembly_filtering_summary.tsv",
        active_species = "results/00_ledger/active_species.tsv"
    params:
        script = str(SCRIPTS_DIR / "filter_and_audit_assembly.py"),
        refs_dir = "results/01_inputs/references",
        min_contig_len = config.get("min_contig_length", 1000),
        min_kept_bp = config.get("min_kept_assembly_bp", 100000),
        max_disc_bp_pct = config.get("max_discarded_bp_pct", 50.0),
        max_disc_contig_pct = config.get("max_discarded_contig_pct", 100.0)
    log:
        "logs/01_filter_reference/filter_and_audit_assemblies.log"
    benchmark:
        "benchmarks/01_filter_reference/filter_and_audit_assemblies.tsv"
    conda:
        str(ENVS_DIR / "core_polars_joltax.yaml")
    threads: get_threads("filter_reference", 4)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("filter_reference", 4000),
        runtime = scale_runtime("filter_reference", 30),
        cpus_per_task = get_threads("filter_reference", 4)
    shell:
        r"""
        python {params.script} \
            --stage0-ledger {input.stage0_ledger} \
            --accession-manifest {input.acc_manifest} \
            --out-references-dir {params.refs_dir} \
            --out-audit-tsv {output.audit_tsv} \
            --out-active-tsv {output.active_species} \
            --min-contig-length {params.min_contig_len} \
            --min-kept-bp {params.min_kept_bp} \
            --max-discarded-bp-pct {params.max_disc_bp_pct} \
            --max-discarded-contig-pct {params.max_disc_contig_pct} \
            --threads {threads} \
            --log {log}
        """
