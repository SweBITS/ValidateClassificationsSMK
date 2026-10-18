#!/usr/bin/env python3
"""
common_ledger.py — Shared Status Ledger, Sanitization & Logging Utilities
=========================================================================

WHY THIS MODULE EXISTS:
-----------------------
In a high-throughput validation pipeline testing dozens or hundreds of putative eDNA
detections simultaneously, individual taxa will naturally fail at different biological
or statistical gates:
  - Stage 0: Too few reads across the 516-week timeseries (< 1,000 read pairs) to power
             a 100-bin Poisson regression (`EXCLUDED_LOW_READ_COUNT`), or represented
             solely by mitochondrial/plastid sequences in the Kraken 2 database
             (`EXCLUDED_ORGANELLAR_ONLY`).
  - Stage 1: Draft assemblies where removing < 1,000 bp contigs leaves < 100 kb of
             nuclear sequence (`EXCLUDED_INSUFFICIENT_CONTIG_LENGTH`) or discards > 50%
             of the assembly's base pairs (`EXCLUDED_FILTERED_ASSEMBLY_PROPORTION`).
  - Stage 5: Insufficient bin occupancy (< 50% of bins with >= 1 mapped read or negative
             initial slope: `FAILED_LOW_OCCUPANCY`) or inability to reach R >= 0.90
             before 50% of Classifiable Region (CR) sequence is pruned
             (`FAILED_LOW_CORRELATION`).

Rather than allowing Snakemake to crash when a taxon fails a pre-check, every target
species in `config/targets.tsv` is tracked through a unified Status Ledger with
standardized status codes and human-readable explanations.
"""

import logging
import re
import sys
from pathlib import Path
from typing import Optional


# ============================================================================ #
# Standardized Validation Status Codes
# ============================================================================ #
class ValidationStatus:
    """Canonical lifecycle and terminal status codes for target taxa."""

    # Intermediate pipeline states
    PASSED_STAGE0 = "PASSED_STAGE0"
    PASSED_STAGE1 = "PASSED_STAGE1"

    # Terminal validated states (Stage 5)
    VALIDATED_CLEAN = "VALIDATED_CLEAN"       # R >= 0.90 & slope > 0 with 0 outlier bins pruned
    VALIDATED_PRUNED = "VALIDATED_PRUNED"     # R >= 0.90 & slope > 0 after pruning >= 1 outlier bins

    # Terminal failed validation states (Stage 5)
    FAILED_LOW_OCCUPANCY = "FAILED_LOW_OCCUPANCY"         # < 50% occupied bins or non-positive slope
    FAILED_LOW_CORRELATION = "FAILED_LOW_CORRELATION"     # R < 0.90 when hitting 50% CR retention floor

    # Terminal pre-flight / reference exclusion states (Stages 0-3)
    EXCLUDED_TAXID_NOT_FOUND = "EXCLUDED_TAXID_NOT_FOUND"
    EXCLUDED_LOW_READ_COUNT = "EXCLUDED_LOW_READ_COUNT"
    EXCLUDED_NO_ASSEMBLY_IN_DB = "EXCLUDED_NO_ASSEMBLY_IN_DB"
    EXCLUDED_ORGANELLAR_ONLY = "EXCLUDED_ORGANELLAR_ONLY"
    EXCLUDED_NCBI_UNAVAILABLE = "EXCLUDED_NCBI_UNAVAILABLE"
    EXCLUDED_INSUFFICIENT_CONTIG_LENGTH = "EXCLUDED_INSUFFICIENT_CONTIG_LENGTH"
    EXCLUDED_FILTERED_ASSEMBLY_PROPORTION = "EXCLUDED_FILTERED_ASSEMBLY_PROPORTION"
    EXCLUDED_INSUFFICIENT_CR_LENGTH = "EXCLUDED_INSUFFICIENT_CR_LENGTH"


def setup_logger(name: str, log_file: Optional[str] = None, level: int = logging.INFO) -> logging.Logger:
    """
    Configure a timestamped logger writing to both stderr and an optional rule log file.

    WHY:
    SLURM jobs running on Kebnekaise need structured timestamps in per-rule log files
    so resource bottlenecks, retry attempts, and biological exclusion reasons can be
    inspected post-mortem.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()

    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def sanitize_species_slug(name: str) -> str:
    """
    Convert a scientific species name into a filesystem-safe slug (e.g. 'Coloeus_monedula').

    WHY:
    Output directories named `{t_id}_{species_slug}` allow researchers browsing `results/`
    to immediately recognize their organism of interest without memorizing NCBI TaxIDs,
    while avoiding spaces or special characters that break shell pipelines.
    """
    cleaned = re.sub(r"[^a-zA-Z0-9_-]+", "_", name.strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return cleaned or "Unknown_species"


def clean_ncbi_accession(raw_accession: str) -> str:
    """
    Strip GTDB prefixes ('RS_' or 'GB_') from an assembly accession string.

    WHY:
    In `database_manifest.tsv`, prokaryotic assemblies sourced from GTDB (`a_source_db == 'GTDB'`)
    store accessions prefixed with `RS_GCF_...` or `GB_GCA_...`. Stripping `RS_`/`GB_`
    recovers the canonical NCBI Assembly accession (`GCF_...` / `GCA_...`) required by
    the NCBI `datasets download genome accession` CLI.
    """
    acc = raw_accession.strip()
    if acc.startswith("RS_") or acc.startswith("GB_"):
        return acc[3:]
    return acc


def sample_id_to_year_week(sample_id: str) -> Optional[str]:
    """
    Convert a Ljungbyhed sample identifier (e.g. 'Lj_2013_5_339') to canonical 'YYYY_WW' ('2013_05').

    WHY:
    Individual classification and FASTQ files are named by filter sample ID (`Lj_{year}_{week}_{id}`),
    whereas the master canonical count matrix (`lj_stringent_canonical_counts.tsv.gz`) indexes the
    516-week time series by zero-padded `YYYY_WW` (`2013_05`). Normalizing here ensures exact
    1-to-1 alignment between raw Kraken counts and mapped CR counts across all 516 weeks.
    """
    m = re.match(r"^Lj_(\d{4})_(\d+)_\d+$", sample_id.strip())
    if m:
        return f"{int(m.group(1)):04d}_{int(m.group(2)):02d}"
    # Also accept already formatted YYYY_WW strings
    m2 = re.match(r"^(\d{4})_(\d{1,2})$", sample_id.strip())
    if m2:
        return f"{int(m2.group(1)):04d}_{int(m2.group(2)):02d}"
    return None
