# ============================================================================ #
# Stage 4: Mapping eDNA Reads Exclusively to CR Contigs (No PCR Deduplication)
# ============================================================================ #

rule map_edna_to_cr:
    """
    WHY:
    Aligns a species' pooled paired-end eDNA reads (`taxID-{t_id}_R{1,2}.fq.gz`) exclusively
    to its excised Classifiable Region reference (`{t_id}.cr_contigs.fa`) using `minimap2 -ax sr`
    without PCR deduplication (`markdup` is intentionally omitted so high-abundance biological
    peaks do not saturate compact CRs).

    Alignments are stream-filtered on the fly (`-f 0x2 -F 0x90C`, `RNEXT == "="`, and BOTH
    mates passing `MAPQ >= 30`, `(qlen - sclen) >= 0.80 * qlen`, `NM <= 0.01 * (qlen - sclen)`
    within the exact same CR) directly into a compact Polars Parquet file without writing
    large intermediate SAM/BAM files to disk.
    """
    input:
        cr_fasta = "results/02_classifiable_regions/{t_id}/{t_id}.cr_contigs.fa",
        r1 = "results/01_inputs/fastq/taxID-{t_id}_R1.fq.gz",
        r2 = "results/01_inputs/fastq/taxID-{t_id}_R2.fq.gz"
    output:
        parquet = "results/03_mapping/{t_id}/{t_id}.passing_alignments.parquet",
        stats_tsv = "results/03_mapping/{t_id}/{t_id}.mapping_stats.tsv",
        id_hist_tsv = "results/03_mapping/{t_id}/{t_id}.identity_histogram_lt99.tsv"
    params:
        script = str(SCRIPTS_DIR / "filter_sam_to_parquet.py"),
        min_mapq = config.get("min_mapq", 30),
        min_align_frac = config.get("min_align_frac", 0.80),
        max_mismatch_frac = config.get("max_mismatch_frac", 0.01)
    log:
        minimap2 = "logs/05_map_edna/minimap2_taxID-{t_id}.log",
        filter_sam = "logs/05_map_edna/filter_sam_taxID-{t_id}.log"
    benchmark:
        "benchmarks/05_map_edna/map_edna_to_cr_taxID-{t_id}.tsv"
    conda:
        str(ENVS_DIR / "mapping.yaml")
    threads: get_threads("map_edna_to_cr", 8)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("map_edna_to_cr", 8000, input_mb_mult=2.0),
        runtime = scale_runtime("map_edna_to_cr", 45),
        cpus_per_task = get_threads("map_edna_to_cr", 8)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.parquet}) $(dirname {log.minimap2})

        # Gracefully handle empty CR reference (if 0 CRs were found in Stage 2/3)
        if [ ! -s {input.cr_fasta} ]; then
            echo "CR FASTA {input.cr_fasta} is empty; emitting empty Parquet table." > {log.minimap2}
            python {params.script} \
                --taxid {wildcards.t_id} \
                --min-mapq {params.min_mapq} \
                --min-align-frac {params.min_align_frac} \
                --max-mismatch-frac {params.max_mismatch_frac} \
                --out-parquet {output.parquet} \
                --out-stats-tsv {output.stats_tsv} \
                --log {log.filter_sam} < /dev/null
            exit 0
        fi

        minimap2 -ax sr -t {threads} {input.cr_fasta} {input.r1} {input.r2} 2> {log.minimap2} \
            | samtools view -h -f 0x2 -F 0x90C -q {params.min_mapq} - 2>> {log.minimap2} \
            | python {params.script} \
                --taxid {wildcards.t_id} \
                --min-mapq {params.min_mapq} \
                --min-align-frac {params.min_align_frac} \
                --max-mismatch-frac {params.max_mismatch_frac} \
                --minimap2-log {log.minimap2} \
                --out-parquet {output.parquet} \
                --out-stats-tsv {output.stats_tsv} \
                --log {log.filter_sam}
        """
