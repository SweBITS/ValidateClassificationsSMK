# ============================================================================ #
# Stage 5: Recursive Outlier Regression, Time Series Reconstruction & Reports
# ============================================================================ #

rule validate_and_regress_species:
    """
    WHY:
    Evaluates the minimum bin occupancy sanity gate (>= 50% of bins with >= 1 mapped read
    pair and positive initial slope), runs bidirectional recursive OLS outlier pruning
    (removing the bin with largest |Y_i - Y_hat_i| until R >= 0.90 or retained CR drops
    below 50%), and reconstructs the 516-week time series (`raw_kraken_reads`,
    `all_cr_mapped_reads`, `validated_clean_bin_reads`).
    """
    input:
        bins_tsv = "results/02_classifiable_regions/{t_id}/{t_id}.bins_metadata.tsv",
        parquet = "results/03_mapping/{t_id}/{t_id}.passing_alignments.parquet",
        map_stats = "results/03_mapping/{t_id}/{t_id}.mapping_stats.tsv",
        weekly_long = "results/01_inputs/raw_weekly_counts/raw_weekly_counts_long.tsv"
    output:
        bin_reg_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.bin_regression_data.tsv",
        traj_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.pruning_trajectory.tsv",
        weekly_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.weekly_timeseries_516w.tsv",
        metrics_json = "results/04_validation/{t_id}_{species_slug}/{t_id}.validation_metrics.json"
    params:
        script = str(SCRIPTS_DIR / "validate_and_regress.py"),
        species_name = lambda wc: get_active_species_map(wc)[str(wc.t_id)]["scientific_name"],
        raw_weekly_tsv = lambda wc: f"results/01_inputs/raw_weekly_counts/{wc.t_id}.raw_weekly.tsv",
        min_occ_frac = config.get("min_occupied_bins_frac", 0.50),
        target_r = config.get("target_r", 0.90),
        min_ret_cr_frac = config.get("min_retained_cr_frac", 0.50)
    log:
        "logs/06_regress_and_report/validate_and_regress_{t_id}_{species_slug}.log"
    benchmark:
        "benchmarks/06_regress_and_report/validate_and_regress_{t_id}_{species_slug}.tsv"
    conda:
        str(ENVS_DIR / "core_polars_joltax.yaml")
    threads: get_threads("validate_and_regress", 2)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("validate_and_regress", 8000),
        runtime = scale_runtime("validate_and_regress", 30),
        cpus_per_task = get_threads("validate_and_regress", 2)
    shell:
        r"""
        python {params.script} \
            --taxid {wildcards.t_id} \
            --species-name "{params.species_name}" \
            --species-slug "{wildcards.species_slug}" \
            --bins-tsv {input.bins_tsv} \
            --alignments-parquet {input.parquet} \
            --mapping-stats-tsv {input.map_stats} \
            --raw-weekly-tsv {params.raw_weekly_tsv} \
            --out-bin-regression-tsv {output.bin_reg_tsv} \
            --out-trajectory-tsv {output.traj_tsv} \
            --out-weekly-tsv {output.weekly_tsv} \
            --out-metrics-json {output.metrics_json} \
            --min-occupied-frac {params.min_occ_frac} \
            --target-r {params.target_r} \
            --min-retained-cr-frac {params.min_ret_cr_frac} \
            --log {log}
        """


rule generate_species_diagnostic_figure:
    """
    WHY:
    Produces the publication-grade 3-panel diagnostic figure for each active species:
      - Panel A: Regression Scatter Plot (Retained circles vs. Pruned triangles labeled by step,
                 Initial dashed regression line vs. Final solid line with 95% CI band).
      - Panel B: Dual-axis Pruning Trajectory Curve (Pearson R vs. % Retained CR Length).
      - Panel C: 516-Week Time Series Comparison (Raw Kraken Reads vs. Validated Clean-Bin Reads).
    """
    input:
        bin_reg_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.bin_regression_data.tsv",
        traj_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.pruning_trajectory.tsv",
        weekly_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.weekly_timeseries_516w.tsv",
        metrics_json = "results/04_validation/{t_id}_{species_slug}/{t_id}.validation_metrics.json"
    output:
        png = "results/04_validation/{t_id}_{species_slug}/{t_id}.diagnostic_3panel.png",
        pdf = "results/04_validation/{t_id}_{species_slug}/{t_id}.diagnostic_3panel.pdf"
    params:
        script = str(SCRIPTS_DIR / "generate_species_plots.py"),
        target_r = config.get("target_r", 0.90),
        min_ret_cr_frac = config.get("min_retained_cr_frac", 0.50)
    log:
        "logs/06_regress_and_report/plot_species_{t_id}_{species_slug}.log"
    benchmark:
        "benchmarks/06_regress_and_report/plot_species_{t_id}_{species_slug}.tsv"
    conda:
        str(ENVS_DIR / "core_polars_joltax.yaml")
    threads: 1
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("validate_and_regress", 4000),
        runtime = scale_runtime("validate_and_regress", 20),
        cpus_per_task = 1
    shell:
        r"""
        python {params.script} \
            --taxid {wildcards.t_id} \
            --bin-regression-tsv {input.bin_reg_tsv} \
            --trajectory-tsv {input.traj_tsv} \
            --weekly-tsv {input.weekly_tsv} \
            --metrics-json {input.metrics_json} \
            --target-r {params.target_r} \
            --min-retained-cr-frac {params.min_ret_cr_frac} \
            --out-png {output.png} \
            --out-pdf {output.pdf} \
            --log {log}
        """


def get_all_active_metrics_jsons(wildcards):
    records = get_active_species_records(wildcards)
    return [
        f"results/04_validation/{r['t_id']}_{r['species_slug']}/{r['t_id']}.validation_metrics.json"
        for r in records
    ]


def get_all_active_diagnostic_pngs(wildcards):
    records = get_active_species_records(wildcards)
    return [
        f"results/04_validation/{r['t_id']}_{r['species_slug']}/{r['t_id']}.diagnostic_3panel.png"
        for r in records
    ]


def get_all_active_diagnostic_pdfs(wildcards):
    records = get_active_species_records(wildcards)
    return [
        f"results/04_validation/{r['t_id']}_{r['species_slug']}/{r['t_id']}.diagnostic_3panel.pdf"
        for r in records
    ]


rule generate_master_validation_report:
    """
    WHY:
    Aggregates lifecycle metrics across EVERY target species in `config/targets.tsv`—including
    taxa excluded at Stage 0 (`EXCLUDED_LOW_READ_COUNT`, `EXCLUDED_ORGANELLAR_ONLY`), Stage 1
    (`EXCLUDED_INSUFFICIENT_CONTIG_LENGTH`, `EXCLUDED_FILTERED_ASSEMBLY_PROPORTION`), or
    Stage 5 (`FAILED_LOW_OCCUPANCY`, `FAILED_LOW_CORRELATION`, `VALIDATED_CLEAN`,
    `VALIDATED_PRUNED`)—into machine-readable Polars TSV/Parquet tables and a self-contained,
    interactive HTML dashboard (`validation_report.html`).
    """
    input:
        wf_info = "results/00_ledger/workflow_run_info.json",
        stage0_ledger = "results/00_ledger/stage0_preflight_and_assemblies.tsv",
        stage1_audit = "results/00_ledger/assembly_filtering_summary.tsv",
        active_species = "results/00_ledger/active_species.tsv",
        cr_summary = "results/02_classifiable_regions/cr_streaming_summary.tsv",
        raw_weekly_long = "results/01_inputs/raw_weekly_counts/raw_weekly_counts_long.tsv",
        metrics_jsons = get_all_active_metrics_jsons,
        diagnostic_pngs = get_all_active_diagnostic_pngs,
        diagnostic_pdfs = get_all_active_diagnostic_pdfs
    output:
        master_tsv = "results/05_reports/validation_master_summary.tsv",
        master_parquet = "results/05_reports/validation_master_summary.parquet",
        timeseries_parquet = "results/05_reports/validated_timeseries_long.parquet",
        html_report = "results/05_reports/validation_report.html"
    params:
        script = str(SCRIPTS_DIR / "generate_master_report.py"),
        template = str(SCRIPTS_DIR / "templates" / "validation_report.html.j2"),
        validation_dir = "results/04_validation",
        benchmarks_dir = "benchmarks",
        config_path = "config/config.yaml"
    log:
        "logs/06_regress_and_report/generate_master_validation_report.log"
    benchmark:
        "benchmarks/06_regress_and_report/generate_master_validation_report.tsv"
    conda:
        str(ENVS_DIR / "core_polars_joltax.yaml")
    threads: get_threads("master_report", 4)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("master_report", 16000),
        runtime = scale_runtime("master_report", 30),
        cpus_per_task = get_threads("master_report", 4)
    shell:
        r"""
        python {params.script} \
            --wf-info {input.wf_info} \
            --config-yaml {params.config_path} \
            --stage0-ledger {input.stage0_ledger} \
            --stage1-audit {input.stage1_audit} \
            --cr-summary {input.cr_summary} \
            --raw-weekly-long {input.raw_weekly_long} \
            --validation-dir {params.validation_dir} \
            --benchmarks-dir {params.benchmarks_dir} \
            --template {params.template} \
            --out-master-tsv {output.master_tsv} \
            --out-master-parquet {output.master_parquet} \
            --out-timeseries-parquet {output.timeseries_parquet} \
            --out-html {output.html_report} \
            --log {log}
        """
