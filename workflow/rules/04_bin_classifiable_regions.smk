# ============================================================================ #
# Stage 3: Shortest-to-Longest Contig Sorting & Linear Variable-Size CR Binning
# ============================================================================ #

def get_species_filtered_fasta(wildcards):
    rec = get_active_species_map(wildcards)[str(wildcards.t_id)]
    return rec["filtered_fasta"]


def get_species_contig_lengths(wildcards):
    rec = get_active_species_map(wildcards)[str(wildcards.t_id)]
    return rec["contig_lengths_tsv"]


def scale_cr_bins_mem_mb(wildcards, attempt):
    base = int(RES_CFG.get("build_cr_bins", {}).get("mem_mb", 4000))
    est_mb = base
    try:
        rec = get_active_species_map(wildcards)[str(wildcards.t_id)]
        kept_mb = int(rec.get("kept_assembly_bp", 0)) / 1_000_000.0
        est_mb = max(base, int(base / 2 + kept_mb * 2.0))
    except Exception:
        est_mb = base
    return int(round(est_mb * (MEM_RETRY_BASE ** (attempt - 1))))


rule build_cr_bins:
    """
    WHY:
    Takes the merged Classifiable Regions (CRs) for a species from the active coverage tier
    (`2x` by default, or `4x`), sorts contigs containing CRs from shortest total contig length
    to longest, partitions the total CR sequence into `n_bins` of linearly increasing target
    capacity ($0.3\\bar{L}$ to $1.7\\bar{L}$), and extracts ONLY those CR sequences into a compact
    reference FASTA (`results/02_classifiable_regions/{BIN_TAG}/{t_id}/{t_id}.cr_contigs.fa`)
    for Stage 4 alignment.

    NOTE: `--cleanup-filtered-fasta` is intentionally omitted so that `{t_id}.{acc}.filtered.fna`
    remains on disk across multiple runs with different `n_bins` or `bin_min_scale`/`bin_max_scale`.
    """
    input:
        active_species = "results/00_ledger/active_species.tsv",
        cr_sentinel = "results/02_classifiable_regions/cr_stream.done",
        contig_lengths = get_species_contig_lengths
    output:
        cr_fasta = f"results/02_classifiable_regions/{BIN_TAG}/{{t_id}}/{{t_id}}.cr_contigs.fa",
        bins_tsv = f"results/02_classifiable_regions/{BIN_TAG}/{{t_id}}/{{t_id}}.bins_metadata.tsv"
    params:
        script = str(SCRIPTS_DIR / "build_cr_bins.py"),
        filtered_fasta = get_species_filtered_fasta,
        cr_bed = lambda wc: f"results/02_classifiable_regions/{wc.t_id}/{wc.t_id}.cr_2x.bed",
        n_bins = config.get("n_bins", 100),
        min_scale = config.get("bin_min_scale", 0.3),
        max_scale = config.get("bin_max_scale", 1.7)
    log:
        f"logs/04_bin_cr/{BIN_TAG}/build_cr_bins_taxID-{{t_id}}.log"
    benchmark:
        f"benchmarks/04_bin_cr/{BIN_TAG}/build_cr_bins_taxID-{{t_id}}.tsv"
    conda:
        str(ENVS_DIR / "core_polars_joltax.yaml")
    threads: get_threads("build_cr_bins", 1)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_cr_bins_mem_mb,
        runtime = scale_runtime("build_cr_bins", 15),
        cpus_per_task = get_threads("build_cr_bins", 1)
    shell:
        r"""
        python {params.script} \
            --taxid {wildcards.t_id} \
            --filtered-fasta {params.filtered_fasta} \
            --contig-lengths {input.contig_lengths} \
            --cr-bed {params.cr_bed} \
            --out-cr-fasta {output.cr_fasta} \
            --out-bins-tsv {output.bins_tsv} \
            --n-bins {params.n_bins} \
            --min-scale {params.min_scale} \
            --max-scale {params.max_scale} \
            --log {log}
        """
