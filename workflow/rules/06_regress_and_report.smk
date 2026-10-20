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
        active_species = "results/00_ledger/active_species.tsv",
        stage1_audit = "results/00_ledger/assembly_filtering_summary.tsv",
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
        mem_mb = scale_mem_mb("validate_and_regress", 4000, input_mb_mult=30.0),
        runtime = scale_runtime("validate_and_regress", 15),
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
    Produces the classic 2x2 4-panel diagnostic figure for each active species:
      - Panel A: Regression Scatter Plot (Retained circles vs. Pruned triangles labeled by step,
                 Initial dashed regression line vs. Final solid line with 95% CI band).
      - Panel B: Dual-axis Pruning Trajectory Curve (Pearson R vs. % Retained CR Length).
      - Panel C: 516-Week Time Series Comparison (Raw Kraken Reads vs. Validated Clean-Bin Reads).
      - Panel D: 1-to-1 Weekly Scatter Plot (Raw Kraken Clade Reads vs. Validated Clean-Bin Reads).
    """
    input:
        bin_reg_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.bin_regression_data.tsv",
        traj_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.pruning_trajectory.tsv",
        weekly_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.weekly_timeseries_516w.tsv",
        metrics_json = "results/04_validation/{t_id}_{species_slug}/{t_id}.validation_metrics.json"
    output:
        png = "results/04_validation/{t_id}_{species_slug}/{t_id}.diagnostic_3panel.png",
        pdf = "results/04_validation/{t_id}_{species_slug}/{t_id}.diagnostic_3panel.pdf",
        png_2x2 = "results/04_validation/{t_id}_{species_slug}/{t_id}.diagnostic_2x2.png",
        pdf_2x2 = "results/04_validation/{t_id}_{species_slug}/{t_id}.diagnostic_2x2.pdf"
    params:
        script = str(SCRIPTS_DIR / "generate_species_plots.py"),
        target_r = config.get("target_r", 0.90),
        min_ret_cr_frac = config.get("min_retained_cr_frac", 0.50)
    log:
        "logs/06_regress_and_report/plot_species_2x2_{t_id}_{species_slug}.log"
    benchmark:
        "benchmarks/06_regress_and_report/plot_species_2x2_{t_id}_{species_slug}.tsv"
    conda:
        str(ENVS_DIR / "core_polars_joltax.yaml")
    threads: get_threads("plot_species", 1)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("plot_species", 2000),
        runtime = scale_runtime("plot_species", 10),
        cpus_per_task = get_threads("plot_species", 1)
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
        cp {output.png} {output.png_2x2}
        cp {output.pdf} {output.pdf_2x2}
        """


rule generate_species_widescreen_figure:
    """
    WHY:
    Produces the extended Widescreen Diagnostic Suite (`32.0 x 18.2 in`, 300 DPI PNG + PDF)
    for each active species on a unified 3x4 column-aligned grid:
      - Panel A: Genomic Proportionality Regression (+ top structural strip & horizontal Cosine colorbar)
      - Panel B: Outlier Pruning Trajectory (Pearson R vs. Retained CR %)
      - Panel C: Genome-Wide CR Coverage Map / Manhattan Plot (symlog y-axis with shaded 0-1 linear band)
      - Panel D: Outlier Classification Quadrant (Excess Read Depth vs. Temporal Concordance, sized by Gini)
      - Panel E1, E2, E3: eDNA Read Attrition Funnel, Alignment Identity Spectrum, & Reference Sequence Funnel
      - Panel F: Weekly Signal Retention (Log10(Reads+1) 1:1 Correlation across 516 weeks)
      - Panel G: 10-Year Longitudinal eDNA Time Series (2013-2022, 516 weeks)
      - Panel H: 10-Year Seasonal Phenology Profile (Weeks 1-52, 3-week circular moving average)
      - Panel I: Bin Temporal Occupancy (% of 516 weeks occupied vs. Cosine Similarity + N_eff)
      - Panel J: Bin x Week Hierarchical Clustermap (100 CR bins x active weeks, Ward linkage + 3 strips)
    """
    input:
        cr_fasta = "results/02_classifiable_regions/{t_id}/{t_id}.cr_contigs.fa",
        parquet = "results/03_mapping/{t_id}/{t_id}.passing_alignments.parquet",
        map_stats = "results/03_mapping/{t_id}/{t_id}.mapping_stats.tsv",
        bin_reg_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.bin_regression_data.tsv",
        traj_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.pruning_trajectory.tsv",
        weekly_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.weekly_timeseries_516w.tsv",
        metrics_json = "results/04_validation/{t_id}_{species_slug}/{t_id}.validation_metrics.json"
    output:
        png_wide = "results/04_validation/{t_id}_{species_slug}/{t_id}.diagnostic_widescreen.png",
        pdf_wide = "results/04_validation/{t_id}_{species_slug}/{t_id}.diagnostic_widescreen.pdf",
        enhanced_tsv = "results/04_validation/{t_id}_{species_slug}/{t_id}.enhanced_bin_metrics.tsv"
    params:
        script = str(SCRIPTS_DIR / "generate_enhanced_species_dashboard.py"),
        target_r = config.get("target_r", 0.90),
        min_ret_cr_frac = config.get("min_retained_cr_frac", 0.50)
    log:
        "logs/06_regress_and_report/plot_species_widescreen_{t_id}_{species_slug}.log"
    benchmark:
        "benchmarks/06_regress_and_report/plot_species_widescreen_{t_id}_{species_slug}.tsv"
    conda:
        str(ENVS_DIR / "core_polars_joltax.yaml")
    threads: get_threads("plot_species_widescreen", 2)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("plot_species_widescreen", 6000, input_mb_mult=15.0),
        runtime = scale_runtime("plot_species_widescreen", 20),
        cpus_per_task = get_threads("plot_species_widescreen", 2)
    shell:
        r"""
        python {params.script} \
            --taxid {wildcards.t_id} \
            --species-slug "{wildcards.species_slug}" \
            --target-r {params.target_r} \
            --min-retained-cr-frac {params.min_ret_cr_frac} \
            --out-png {output.png_wide} \
            --out-pdf {output.pdf_wide} \
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
        f"results/04_validation/{r['t_id']}_{r['species_slug']}/{r['t_id']}.diagnostic_2x2.png"
        for r in records
    ]


def get_all_active_diagnostic_pdfs(wildcards):
    records = get_active_species_records(wildcards)
    return [
        f"results/04_validation/{r['t_id']}_{r['species_slug']}/{r['t_id']}.diagnostic_2x2.pdf"
        for r in records
    ]


def get_all_active_widescreen_pngs(wildcards):
    records = get_active_species_records(wildcards)
    return [
        f"results/04_validation/{r['t_id']}_{r['species_slug']}/{r['t_id']}.diagnostic_widescreen.png"
        for r in records
    ]


def get_all_active_widescreen_pdfs(wildcards):
    records = get_active_species_records(wildcards)
    return [
        f"results/04_validation/{r['t_id']}_{r['species_slug']}/{r['t_id']}.diagnostic_widescreen.pdf"
        for r in records
    ]


rule generate_master_validation_report:
    """
    WHY:
    Aggregates lifecycle metrics across EVERY target species in `config/targets.tsv` into
    machine-readable Polars TSV/Parquet tables and BOTH self-contained interactive HTML dashboards:
      1. `validation_report_2x2.html` (and `validation_report.html`): embeds the compact 2x2 figures.
      2. `validation_report_widescreen.html`: embeds the widescreen 10-panel diagnostic figures.
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
        diagnostic_pdfs = get_all_active_diagnostic_pdfs,
        widescreen_pngs = get_all_active_widescreen_pngs,
        widescreen_pdfs = get_all_active_widescreen_pdfs
    output:
        master_tsv = "results/05_reports/validation_master_summary.tsv",
        master_parquet = "results/05_reports/validation_master_summary.parquet",
        timeseries_parquet = "results/05_reports/validated_timeseries_long.parquet",
        html_report = "results/05_reports/validation_report.html",
        html_report_2x2 = "results/05_reports/validation_report_2x2.html",
        html_report_wide = "results/05_reports/validation_report_widescreen.html"
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
            --out-html-2x2 {output.html_report_2x2} \
            --out-html-widescreen {output.html_report_wide} \
            --log {log}
        """
