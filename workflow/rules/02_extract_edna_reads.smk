# ============================================================================ #
# Stage 0/2: Single-Pass Multi-Taxon eDNA Read Extraction Across 516 Samples
# ============================================================================ #

rule extract_sample_clade_reads:
    """
    WHY:
    Scans a single weekly sample's `.kraken2.gz` classification output at C speed (`pigz + awk`)
    for reads assigned to any TaxID in any active species' JolTax clade, runs `seqtk subseq`
    once on R1 and once on R2 for the union of matching read IDs, and demultiplexes into
    per-species chunks with headers strictly rewritten to `@<read_id>|sample=<sample_id>`.
    """
    input:
        active_species = "results/00_ledger/active_species.tsv",
        clade_map = "results/00_ledger/stage0_clade_map.tsv",
        kraken_out = lambda wc: config["sample_kraken_pattern"].format(sample=wc.sample),
        r1 = lambda wc: config["raw_fastq_pattern"].format(sample=wc.sample, pair="1"),
        r2 = lambda wc: config["raw_fastq_pattern"].format(sample=wc.sample, pair="2")
    output:
        done = temp("results/01_inputs/sample_chunks/{sample}/{sample}.extract.done")
    params:
        script = str(SCRIPTS_DIR / "extract_sample_clade_reads.py"),
        out_dir = "results/01_inputs/sample_chunks/{sample}"
    log:
        "logs/02_extract_reads/{sample}.extract.log"
    benchmark:
        "benchmarks/02_extract_reads/{sample}.extract.tsv"
    conda:
        str(ENVS_DIR / "extract_reads.yaml")
    threads: get_threads("extract_sample_reads", 2)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_sample_extract_mem_mb("extract_sample_reads", 4000),
        runtime = scale_runtime("extract_sample_reads", 60),
        cpus_per_task = get_threads("extract_sample_reads", 2)
    shell:
        r"""
        python {params.script} \
            --sample-id {wildcards.sample} \
            --kraken-out {input.kraken_out} \
            --r1 {input.r1} \
            --r2 {input.r2} \
            --active-species {input.active_species} \
            --clade-map {input.clade_map} \
            --out-dir {params.out_dir} \
            --out-done {output.done} \
            --threads {threads} \
            --log {log}
        """


rule merge_species_fastq:
    """
    WHY:
    Concatenates all non-empty per-sample gzipped FASTQ chunks for a target species (`{t_id}`)
    into a single pooled R1/R2 pair (`taxID-{t_id}_R{1,2}.fq.gz`) using RFC 1952 multi-member
    gzip concatenation (zero decompression/recompression CPU cost), then removes the temporary
    per-sample chunk files for that species to keep disk usage minimal.
    """
    input:
        active_species = "results/00_ledger/active_species.tsv",
        sample_dones = expand("results/01_inputs/sample_chunks/{sample}/{sample}.extract.done", sample=SAMPLES)
    output:
        r1 = "results/01_inputs/fastq/taxID-{t_id}_R1.fq.gz",
        r2 = "results/01_inputs/fastq/taxID-{t_id}_R2.fq.gz"
    params:
        chunks_root = "results/01_inputs/sample_chunks"
    log:
        "logs/02_extract_reads/merge_species_taxID-{t_id}.log"
    benchmark:
        "benchmarks/02_extract_reads/merge_species_taxID-{t_id}.tsv"
    conda:
        str(ENVS_DIR / "extract_reads.yaml")
    threads: get_threads("merge_species_fastq", 1)
    retries: DEFAULT_RETRIES
    resources:
        mem_mb = scale_mem_mb("merge_species_fastq", 1000),
        runtime = scale_runtime("merge_species_fastq", 10),
        cpus_per_task = get_threads("merge_species_fastq", 1)
    shell:
        r"""
        (
        set -euo pipefail
        mkdir -p $(dirname {output.r1})

        # Concatenate all non-empty (>28 bytes) sample R1 chunks
        r1_list=$(mktemp "${{SNIC_TMP:-/tmp}}/r1_chunks_{wildcards.t_id}_XXXXXX.txt")
        r2_list=$(mktemp "${{SNIC_TMP:-/tmp}}/r2_chunks_{wildcards.t_id}_XXXXXX.txt")
        trap 'rm -f "$r1_list" "$r2_list"' EXIT

        find {params.chunks_root} -name "*_taxID-{wildcards.t_id}_R1.fq.gz" -size +28c | sort > "$r1_list"
        find {params.chunks_root} -name "*_taxID-{wildcards.t_id}_R2.fq.gz" -size +28c | sort > "$r2_list"

        if [ -s "$r1_list" ]; then
            xargs -a "$r1_list" cat > {output.r1}
        else
            printf "" | gzip -c > {output.r1}
        fi

        if [ -s "$r2_list" ]; then
            xargs -a "$r2_list" cat > {output.r2}
        else
            printf "" | gzip -c > {output.r2}
        fi

        # Clean up individual sample chunks for this species after successful concatenation
        find {params.chunks_root} -name "*_taxID-{wildcards.t_id}_R1.fq.gz" -delete
        find {params.chunks_root} -name "*_taxID-{wildcards.t_id}_R2.fq.gz" -delete
        echo "Merged $(wc -l < "$r1_list") non-empty weekly chunks for TaxID {wildcards.t_id}"
        ) > {log} 2>&1
        """
