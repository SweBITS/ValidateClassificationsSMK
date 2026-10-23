# ============================================================================ #
# Stage 2: Zero-Disk Streaming 400 bp / 200 bp-Step (2x) Tiling & Kraken 2 CRs
# ============================================================================ #

rule stream_kraken_cr:
    """
    WHY THIS RULE USES A ZERO-DISK POSIX FIFO, PARALLEL TILER, AND 36 KRAKEN THREADS:
    1. Loading the 2.6 TB Kraken 2 database (`hash.k2d`) into RAM on a 72-core `cpu_largemem`
       node is the single heaviest I/O operation in the workflow. Streaming all active filtered
       reference genomes through a single Kraken 2 invocation ensures the 2.6 TB database is
       loaded into memory only ONCE.
    2. We tile genomes into 400 bp windows stepped every 200 bp (2x coverage) using a multi-worker
       Python producer (`stream_tiled_fasta.py`, 8 workers) feeding a 1 MB-buffered POSIX FIFO
       (`$SNIC_TMP/sim_reads.fa`). This prevents Kraken 2's `#pragma omp critical(seqread)`
       (`LoadBlock(3 MB)`) from ever blocking on single-threaded regex header parsing.
    3. On Kebnekaise's 4-socket x 18-core (72-core) Broadwell `largemem` nodes, running Kraken 2
       with 36 OpenMP threads (2 sockets) avoids 4-socket NUMA interconnect saturation and OpenMP
       critical-section lock contention (`seqread`, `output_queue`, `update_taxon_counters`), while
       leaving the remaining cores for the parallel FIFO tiler, `awk`, and CR consumer.
    4. We load `intel/2025a` (`GCCcore/14.2.0` runtime toolchain matching `kraken2_SweBiTS_MHG_intel`)
       and explicitly export `OMP_NUM_THREADS={params.kraken_threads}` to override Kebnekaise's
       default `OMP_NUM_THREADS=1`.
    """
    input:
        active_species = "results/00_ledger/active_species.tsv",
        clade_map = "results/00_ledger/stage0_clade_map.tsv"
    output:
        summary_tsv = "results/02_classifiable_regions/cr_streaming_summary.tsv",
        sentinel = "results/02_classifiable_regions/cr_stream.done"
    params:
        tiler_script = str(SCRIPTS_DIR / "stream_tiled_fasta.py"),
        consumer_script = str(SCRIPTS_DIR / "consume_kraken_stream_to_cr.py"),
        kraken2_bin = (
            str(Path(config["kraken2_bin_dir"]) / "kraken2")
            if config.get("kraken2_bin_dir")
            else "kraken2"
        ),
        toolchain_module = config.get("kraken2_toolchain_module", "intel/2025a"),
        db_path = config["kraken2_db_path"],
        confidence = config.get("kraken2_confidence", 0.1),
        mhg = config.get("kraken2_minimum_hit_groups", 10),
        window_bp = config.get("tile_window_bp", 400),
        step_bp = config.get("tile_step_bp", 200),
        pad_bp = config.get("tile_step_bp", 200),
        kraken_threads = int(RES_CFG.get("stream_kraken_cr", {}).get("kraken_threads", 36)),
        tiler_workers = int(RES_CFG.get("stream_kraken_cr", {}).get("tiler_workers", 8)),
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
        runtime = f"{int(RES_CFG.get('stream_kraken_cr', {}).get('runtime_min', 360))}m",
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
                --pad-bp {params.pad_bp} \
                --log {log.consumer} < /dev/null
            exit 0
        fi

        # Preserve Conda bin directory, then load the toolchain module (intel/2025a) used to compile Kraken 2
        CONDA_BIN_DIR="$(dirname "$(command -v python)")"
        if command -v module >/dev/null 2>&1 || [ -f /hpc2n/eb/software/lmod/lmod/init/bash ]; then
            if ! command -v module >/dev/null 2>&1; then
                set +u
                source /hpc2n/eb/software/lmod/lmod/init/bash
                set -u
            fi
            set +u
            module purge > /dev/null 2>&1 || true
            module load {params.toolchain_module} > /dev/null 2>&1 || true
            set -u
        fi
        export PATH="$CONDA_BIN_DIR:$PATH"

        # Release Kebnekaise's default OMP_NUM_THREADS=1 and set optimal Kraken 2 thread count
        export OMP_NUM_THREADS={params.kraken_threads}

        FIFO_DIR=$(mktemp -d "${{SNIC_TMP:-/tmp}}/kraken_fifo_XXXXXX")
        FIFO_PATH="$FIFO_DIR/sim_reads.fa"
        PRODUCER_OK="$FIFO_DIR/producer.ok"
        mkfifo "$FIFO_PATH"

        cleanup() {{
            rm -rf "$FIFO_DIR"
        }}
        trap cleanup EXIT

        # Launch background multi-worker FIFO producer streaming 400bp / 200bp-step (2x) windows
        (
            set -euo pipefail
            python {params.tiler_script} \
                --active-species {input.active_species} \
                --window-bp {params.window_bp} \
                --step-bp {params.step_bp} \
                --workers {params.tiler_workers} > "$FIFO_PATH"
            touch "$PRODUCER_OK"
        ) 2>> {log.kraken} &
        PRODUCER_PID=$!

        # Run single-pass Kraken 2 (WITHOUT --report) -> C-speed awk -> Python 2x CR consumer
        > {log.kraken}
        /usr/bin/time -v {params.kraken2_bin} \
            --db {params.db_path} \
            --threads {params.kraken_threads} \
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
                --pad-bp {params.pad_bp} \
                --log {log.consumer}

        # Explicitly wait for background FIFO producer and verify clean completion
        wait "$PRODUCER_PID"
        if [ ! -f "$PRODUCER_OK" ]; then
            echo "ERROR: Background FIFO producer did not complete cleanly." >> {log.kraken}
            exit 1
        fi
        grep -q "^STREAM_DRAINED_AND_FSYNCED_OK$" {output.sentinel}
        """
