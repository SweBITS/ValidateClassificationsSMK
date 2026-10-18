# ============================================================================ #
# Stage 0: Clade Resolution (JolTax), Pre-Flight Gates & Best-Assembly Selection
# ============================================================================ #

rule workflow_start_info:
    output:
        info_json = "results/00_ledger/workflow_run_info.json"
    log:
        "logs/00_preflight/workflow_start_info.log"
    shell:
        r"""
        python -c '
import json, subprocess, datetime
from pathlib import Path
snk_ver = subprocess.check_output(["snakemake", "--version"], text=True).strip()
git_rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip() or "unknown"
git_branch = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True, text=True).stdout.strip() or "unknown"
Path("{output.info_json}").parent.mkdir(parents=True, exist_ok=True)
Path("{output.info_json}").write_text(json.dumps({{
    "start_timestamp": datetime.datetime.now().astimezone().isoformat(),
    "snakemake_version": snk_ver,
    "git_branch": git_branch,
    "git_commit": git_rev
}}, indent=2))
' > {log} 2>&1
        """


rule resolve_clades_and_assemblies:
    """
    WHY:
    Expands each target species TaxID in `config/targets.tsv` into its full descendant clade
    using JolTax, queries `lj_stringent_canonical_counts.tsv.gz` for Pre-Flight Gate 1
    (>= 1,000 total clade reads across 516 weeks) and baseline weekly read counts, checks
    `database_manifest.tsv` for Pre-Flight Gate 2 (excluding organellar-only clades), and
    applies the 4-tier deterministic cascade to select the single best reference assembly.
    """
    input:
        targets = config["targets_tsv"],
        manifest = config["database_manifest"],
        canonical_counts = config["merged_canonical_counts_tsv"]
    output:
        ledger = "results/00_ledger/stage0_preflight_and_assemblies.tsv",
        accessions = "results/01_inputs/ncbi_datasets/selected_accessions.txt",
        clade_map = "results/00_ledger/stage0_clade_map.tsv",
        weekly_long = "results/01_inputs/raw_weekly_counts/raw_weekly_counts_long.tsv"
    params:
        script = str(SCRIPTS_DIR / "resolve_clades_and_assemblies.py"),
        joltax_cache = config["joltax_cache_dir"],
        min_reads = config.get("min_total_reads", 1000),
        force_taxids = " ".join(str(x) for x in config.get("force_include_taxids", [])),
        seed = config.get("random_seed", 42),
        weekly_dir = "results/01_inputs/raw_weekly_counts"
    log:
        "logs/00_preflight/resolve_clades_and_assemblies.log"
    benchmark:
        "benchmarks/00_preflight/resolve_clades_and_assemblies.tsv"
    conda:
        str(ENVS_DIR / "core_polars_joltax.yaml")
    threads: get_threads("resolve_clades", 2)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("resolve_clades", 4000),
        runtime = scale_runtime("resolve_clades", 15),
        cpus_per_task = get_threads("resolve_clades", 2)
    shell:
        r"""
        python {params.script} \
            --targets {input.targets} \
            --manifest {input.manifest} \
            --joltax-cache {params.joltax_cache} \
            --canonical-counts {input.canonical_counts} \
            --min-total-reads {params.min_reads} \
            --force-include-taxids {params.force_taxids} \
            --random-seed {params.seed} \
            --out-ledger {output.ledger} \
            --out-accessions {output.accessions} \
            --out-clade-map {output.clade_map} \
            --out-weekly-dir {params.weekly_dir} \
            --log {log}
        """


rule download_and_rehydrate_assemblies:
    """
    WHY:
    Fetches all Stage 0 selected reference genomes in a single bulk NCBI `datasets`
    dehydrated download (`--dehydrated`), followed by multi-threaded rehydration and
    parallel MD5 checksum verification with automatic retry on corrupt/truncated files.
    """
    input:
        accessions = "results/01_inputs/ncbi_datasets/selected_accessions.txt"
    output:
        manifest = "results/01_inputs/ncbi_datasets/accession_to_fasta.tsv",
        flag = "results/01_inputs/ncbi_datasets/rehydrate.done"
    params:
        script = str(SCRIPTS_DIR / "rehydrate_assemblies.py"),
        work_dir = "results/01_inputs/ncbi_datasets",
        max_attempts = config.get("rehydrate_max_attempts", 5),
        api_key = config.get("ncbi_api_key", "")
    log:
        "logs/00_preflight/download_and_rehydrate_assemblies.log"
    benchmark:
        "benchmarks/00_preflight/download_and_rehydrate_assemblies.tsv"
    conda:
        str(ENVS_DIR / "ncbi_datasets.yaml")
    threads: get_threads("download_rehydrate", 16)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("download_rehydrate", 8000),
        runtime = scale_runtime("download_rehydrate", 180),
        cpus_per_task = get_threads("download_rehydrate", 16)
    shell:
        r"""
        python {params.script} \
            --accessions-file {input.accessions} \
            --work-dir {params.work_dir} \
            --out-manifest {output.manifest} \
            --out-flag {output.flag} \
            --max-attempts {params.max_attempts} \
            --threads {threads} \
            --api-key "{params.api_key}" \
            --log {log}
        """
