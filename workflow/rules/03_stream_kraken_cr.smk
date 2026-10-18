# ============================================================================ #
# Stage 2: Zero-Disk Streaming 400 bp Tiling, Single-Pass Kraken 2 & Flanked CRs
# ============================================================================ #

rule stream_kraken_cr:
    """
    WHY THIS RULE USES A ZERO-DISK POSIX FIFO AND SINGLE-PASS KRAKEN 2 STREAMING:
    1. Loading the 2.6 TB Kraken 2 database (`hash.k2d`) into RAM on a 72-core `cpu_largemem`
       node is the single heaviest I/O operation in the workflow. Streaming all active filtered
       reference genomes through a single Kraken 2 invocation ensures the 2.6 TB database is
       loaded into memory only ONCE.
    2. Writing tens of millions of 400 bp simulated reads (4x coverage) and their massive
       Kraken 2 k-mer hit strings to disk would consume hundreds of gigabytes of Lustre I/O.
       Instead, `seqkit sliding -W 400 -s 100 -g` feeds a POSIX FIFO (`$SNIC_TMP/sim_reads.fa`),
       and Kraken 2 stdout is piped directly into C-speed `awk -F'\\t' '$1 == "C" {print $2 "\\t" $3}'`
       before entering `consume_kraken_stream_to_cr.py`.
    3. We intentionally omit `--report` on Kraken 2: summarizing taxonomy trees across 2.6 TB
       of minimizers at job completion adds substantial walltime with zero benefit for coordinate
       discovery.
    4. `retries: 0` is enforced because this rule already allocates a full `largemem` node.
    """
    input:
        active_species = "results/00_ledger/active_species.tsv",
        clade_map = "results/00_ledger/stage0_clade_map.tsv"
    output:
        summary_tsv = "results/02_classifiable_regions/cr_streaming_summary.tsv",
        sentinel = "results/02_classifiable_regions/cr_stream.done"
    params:
        consumer_script = str(SCRIPTS_DIR / "consume_kraken_stream_to_cr.py"),
        kraken2_bin = (
            str(Path(config["kraken2_bin_dir"]) / "kraken2")
            if config.get("kraken2_bin_dir")
            else "kraken2"
        ),
        db_path = config["kraken2_db_path"],
        confidence = config.get("kraken2_confidence", 0.1),
        mhg = config.get("kraken2_minimum_hit_groups", 10),
        window_bp = config.get("tile_window_bp", 400),
        step_4x_bp = config.get("tile_step_4x_bp", 100),
        pad_2x = config.get("tile_step_4x_bp", 100) * 2,
        pad_4x = config.get("tile_step_4x_bp", 100),
        out_cr_dir = "results/02_classifiable_regions"
    log:
        kraken = "logs/03_stream_kraken_cr/kraken2_stream.log",
        consumer = "logs/03_stream_kraken_cr/consume_kraken_stream_to_cr.log"
    benchmark:
        "benchmarks/03_stream_kraken_cr/stream_kraken_cr.tsv"
    conda:
        str(ENVS_DIR / "seqkit_kraken.yaml")
    threads: get_threads("stream_kraken_cr", 72)
    retries: 0
    resources:
        mem_mb = int(RES_CFG.get("stream_kraken_cr", {}).get("mem_mb", 2800000)),
        runtime = f"{int(RES_CFG.get('stream_kraken_cr', {}).get('runtime_min', 1800))}m",
        cpus_per_task = get_threads("stream_kraken_cr", 72),
        constraint = RES_CFG.get("stream_kraken_cr", {}).get("constraint", "largemem")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.summary_tsv}) $(dirname {log.kraken})

        # Check if any active species exist after Stage 0 + Stage 1
        n_active=$(awk 'NR > 1 && NF > 0 {{c++}} END {{print c+0}}' {input.active_species})
        if [ "$n_active" -eq 0 ]; then
            echo "No active species passed Stage 0 + Stage 1; skipping Kraken 2 streaming." > {log.kraken}
            python {params.consumer_script} \
                --active-species {input.active_species} \
                --clade-map {input.clade_map} \
                --out-cr-dir {params.out_cr_dir} \
                --out-summary-tsv {output.summary_tsv} \
                --out-sentinel {output.sentinel} \
                --pad-2x {params.pad_2x} \
                --pad-4x {params.pad_4x} \
                --log {log.consumer} < /dev/null
            exit 0
        fi

        FIFO_DIR=$(mktemp -d "${{SNIC_TMP:-/tmp}}/kraken_fifo_XXXXXX")
        FIFO_PATH="$FIFO_DIR/sim_reads.fa"
        PRODUCER_OK="$FIFO_DIR/producer.ok"
        mkfifo "$FIFO_PATH"

        cleanup() {{
            rm -rf "$FIFO_DIR"
        }}
        trap cleanup EXIT

        # Launch background FIFO producer streaming 400bp / 100bp-step (4x) windows across all active references
        (
            set -euo pipefail
            tail -n +2 {input.active_species} | while IFS=$'\t' read -r tid t_sp sci slug acc clean_acc fasta clen_tsv rest; do
                [ -z "$tid" ] && continue
                seqkit sliding -W {params.window_bp} -s {params.step_4x_bp} -g -j 4 "$fasta" \
                    | awk -v tid="$tid" -v step={params.step_4x_bp} -v win={params.window_bp} '
                        /^>/ {{
                            hdr = substr($1, 2);
                            pos = match(hdr, /_sliding:[0-9]+-[0-9]+$/);
                            contig = substr(hdr, 1, pos - 1);
                            coords = substr(hdr, pos + 9);
                            split(coords, c, "-");
                            s0 = c[1] - 1;
                            e0 = c[2] + 0;
                            if ((e0 - s0) < win && contig == prev_contig && e0 == prev_e0) {{
                                skip = 1;
                                next
                            }}
                            skip = 0;
                            prev_contig = contig;
                            prev_e0 = e0;
                            widx = int(s0 / step);
                            print ">" tid "|" contig "|" widx "|" s0 "|" e0;
                            next
                        }}
                        !skip {{ print }}
                    '
            done > "$FIFO_PATH"
            touch "$PRODUCER_OK"
        ) 2>> {log.kraken} &
        PRODUCER_PID=$!

        # Run single-pass Kraken 2 (WITHOUT --report) -> C-speed awk -> Python dual-tier CR consumer
        {params.kraken2_bin} \
            --db {params.db_path} \
            --threads {threads} \
            --confidence {params.confidence} \
            --minimum-hit-groups {params.mhg} \
            "$FIFO_PATH" 2>> {log.kraken} \
            | awk -F'\t' '$1 == "C" {{print $2 "\t" $3}}' \
            | python {params.consumer_script} \
                --active-species {input.active_species} \
                --clade-map {input.clade_map} \
                --out-cr-dir {params.out_cr_dir} \
                --out-summary-tsv {output.summary_tsv} \
                --out-sentinel {output.sentinel} \
                --pad-2x {params.pad_2x} \
                --pad-4x {params.pad_4x} \
                --log {log.consumer}

        # Explicitly wait for background FIFO producer and verify clean completion
        wait "$PRODUCER_PID"
        if [ ! -f "$PRODUCER_OK" ]; then
            echo "ERROR: Background FIFO producer did not complete cleanly." >> {log.kraken}
            exit 1
        fi
        grep -q "^STREAM_DRAINED_AND_FSYNCED_OK$" {output.sentinel}
        """
