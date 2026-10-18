# ============================================================================ #
# Stage 3: Shortest-to-Longest Contig Sorting & Linear Variable-Size CR Binning
# ============================================================================ #

def get_species_filtered_fasta(wildcards):
    rec = get_active_species_map(wildcards)[str(wildcards.t_id)]
    return rec["filtered_fasta"]


def get_species_contig_lengths(wildcards):
    rec = get_active_species_map(wildcards)[str(wildcards.t_id)]
    return rec["contig_lengths_tsv"]


rule build_cr_bins:
    """
    WHY:
    Takes the merged Classifiable Regions (CRs) for a species from the active coverage tier
    (`2x` by default, or `4x`), sorts contigs containing CRs from shortest total contig length
    to longest, partitions the total CR sequence into 100 bins of linearly increasing target
    capacity ($0.3\\bar{L}$ to $1.7\\bar{L}$), and extracts ONLY those CR sequences into a compact
    reference FASTA (`{t_id}.cr_contigs.fa`) for Stage 4 alignment.
    """
    input:
        cr_sentinel = "results/02_classifiable_regions/cr_stream.done",
        filtered_fasta = get_species_filtered_fasta,
        contig_lengths = get_species_contig_lengths
    output:
        cr_fasta = "results/02_classifiable_regions/{t_id}/{t_id}.cr_contigs.fa",
        bins_tsv = "results/02_classifiable_regions/{t_id}/{t_id}.bins_metadata.tsv"
    params:
        script = str(SCRIPTS_DIR / "build_cr_bins.py"),
        tier = config.get("coverage_tier", "2x"),
        cr_bed = lambda wc: f"results/02_classifiable_regions/{wc.t_id}/{wc.t_id}.cr_{config.get('coverage_tier', '2x')}.bed",
        n_bins = config.get("n_bins", 100),
        min_scale = config.get("bin_min_scale", 0.3),
        max_scale = config.get("bin_max_scale", 1.7)
    log:
        "logs/04_bin_cr/build_cr_bins_taxID-{t_id}.log"
    benchmark:
        "benchmarks/04_bin_cr/build_cr_bins_taxID-{t_id}.tsv"
    conda:
        str(ENVS_DIR / "core_polars_joltax.yaml")
    threads: get_threads("build_cr_bins", 2)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("build_cr_bins", 8000),
        runtime = scale_runtime("build_cr_bins", 30),
        cpus_per_task = get_threads("build_cr_bins", 2)
    shell:
        r"""
        python {params.script} \
            --taxid {wildcards.t_id} \
            --filtered-fasta {input.filtered_fasta} \
            --contig-lengths {input.contig_lengths} \
            --cr-bed {params.cr_bed} \
            --out-cr-fasta {output.cr_fasta} \
            --out-bins-tsv {output.bins_tsv} \
            --n-bins {params.n_bins} \
            --min-scale {params.min_scale} \
            --max-scale {params.max_scale} \
            --log {log}
        """
