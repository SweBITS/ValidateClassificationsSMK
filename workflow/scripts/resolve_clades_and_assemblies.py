#!/usr/bin/env python3
"""
resolve_clades_and_assemblies.py — Stage 0: Clade Resolution, Pre-Flight Gates & Best-Assembly Selection
========================================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
Before launching expensive multi-terabyte Kraken 2 simulations or extracting reads across
516 metagenomic samples, we must answer three foundational questions for every target
species in `config/targets.tsv`:

1. **What is the complete taxonomic clade (`JolTax`)?**
   Kraken 2 may assign reads not only to the species TaxID (`t_id`) itself, but also to
   subordinate subspecies, strains, varietas, or formae speciales within that species.
   Using `JolTax`, we expand each target `t_id` into its full descendant clade so that
   both read counting/extraction and in silico CR classification operate at the complete
   species-clade level.

2. **Does the species pass Pre-Flight Gate 1 (Minimum Read Count >= 1,000)?**
   When a reference genome's Classifiable Regions (CRs) are partitioned into 100 bins
   linearly scaled from 0.3x to 1.7x the mean bin length, the smallest bin holds 0.3%
   of total CR sequence. Under a uniform Poisson model, a minimum of 1,000 total reads
   ensures expected count lambda >= 3 in the smallest bin (>= 95% probability of being
   non-zero), preventing zero-inflation artifacts during linear regression. We query
   `lj_stringent_canonical_counts.tsv.gz` to check total reads and simultaneously extract
   the baseline 516-week `raw_kraken_reads` time series.

3. **Does the species pass Pre-Flight Gate 2 (Non-Organellar Reference) and which single
   assembly in the Kraken 2 database represents this clade?**
   - If a clade is represented in `database_manifest.tsv` solely by mitochondrial,
     chloroplast, plastid, or apicoplast sequences (`EXCLUDED_ORGANELLAR_ONLY`), whole-genome
     uniformity validation is biologically invalid because all reads map to a tiny
     16–150 kb multicopy organelle.
   - Because our Kraken 2 database includes at most one assembly per taxonomic node, a
     species clade may contain an assembly at the species node (`t_id`) plus assemblies
     at subspecies/strain nodes. We select the single best non-organellar `current`
     reference assembly via a deterministic 4-tier cascade:
       (1) `a_refseq_category` in {'reference genome', 'representative genome'}
       (2) Highest `a_level` (Complete Genome > Chromosome > Scaffold > Contig)
       (3) Species-node preference (`assembly.t_id == target_t_id`)
       (4) Highest N50 -> Lowest sequence count -> Latest `a_release_date` ->
           Deterministic seeded random choice (`random_state=42`).
"""

import argparse
import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import joltax

# Allow importing sibling module common_ledger.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import (
    ValidationStatus,
    clean_ncbi_accession,
    extract_accession_from_override,
    is_valid_ncbi_accession,
    resolve_local_genome_file,
    sanitize_species_slug,
    setup_logger,
)

# Canonical ordering of NCBI assembly levels (lower rank integer = higher contiguity)
ASSEMBLY_LEVEL_ORDER: Dict[str, int] = {
    "Complete Genome": 0,
    "Chromosome": 1,
    "Scaffold": 2,
    "Contig": 3,
}

REFSEQ_PRIORITY_CATEGORIES: Set[str] = {"reference genome", "representative genome"}
ORGANELLE_TYPES: Set[str] = {"mitochondrion", "chloroplast", "plastid", "apicoplast"}


def is_organellar_row(row: Dict[str, str]) -> bool:
    """
    Return True if a manifest row represents an organellar sequence rather than a nuclear genome.
    """
    source_db = row.get("a_source_db", "").strip().upper()
    org_type = row.get("a_organelle_type", "").strip().lower()
    return source_db == "ORGANELLE" or org_type in ORGANELLE_TYPES or bool(org_type)


def parse_float_safe(val: str) -> Optional[float]:
    """Parse a float from a TSV string, returning None for empty or 'NA' values."""
    v = (val or "").strip()
    if not v or v.upper() in {"NA", "NAN", "NONE", "NULL"}:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def parse_date_safe(val: str) -> Optional[datetime]:
    """Parse an ISO YYYY-MM-DD release date, returning None if missing or malformed."""
    v = (val or "").strip()
    if not v or v.upper() in {"NA", "NAN", "NONE", "NULL"}:
        return None
    for fmt in ("%Y-%m-%d", "%Y/%m/%d", "%Y-%m-%dT%H:%M:%S"):
        try:
            return datetime.strptime(v[:10], fmt)
        except ValueError:
            continue
    return None


def select_best_assembly_for_clade(
    target_tid: int,
    candidates: List[Dict[str, str]],
    random_seed: int = 42,
) -> Tuple[Dict[str, str], str]:
    """
    Deterministically select a single reference assembly from all `current` non-organellar
    assemblies belonging to a target species clade.

    Cascade:
      1. `a_refseq_category` in {'reference genome', 'representative genome'}
      2. Highest `a_level` (Complete Genome > Chromosome > Scaffold > Contig)
      3. Direct assignment to the target species node (`int(row['t_id']) == target_tid`)
      4. Highest N50 (`a_contig_n50` for Contig level, else `a_scaffold_n50`)
      5. Lowest sequence count (`a_num_contigs` for Contig level, else `a_num_scaffolds`)
      6. Most recent `a_release_date`
      7. Deterministic seeded random choice (`random.Random(random_seed + target_tid)`)
    """
    if len(candidates) == 1:
        return candidates[0], "only_assembly_available"

    reasons: List[str] = []
    pool = list(candidates)

    # Step 1: RefSeq category preference
    refseq_pool = [
        r for r in pool if r.get("a_refseq_category", "").strip().lower() in REFSEQ_PRIORITY_CATEGORIES
    ]
    if refseq_pool:
        pool = refseq_pool
        reasons.append("refseq_category")
        if len(pool) == 1:
            return pool[0], "+".join(reasons)

    # Step 2: Highest assembly level
    best_level_rank = min(
        ASSEMBLY_LEVEL_ORDER.get(r.get("a_level", "").strip(), 999) for r in pool
    )
    level_pool = [
        r for r in pool if ASSEMBLY_LEVEL_ORDER.get(r.get("a_level", "").strip(), 999) == best_level_rank
    ]
    pool = level_pool
    reasons.append("assembly_level")
    if len(pool) == 1:
        return pool[0], "+".join(reasons)

    # Step 3: Species-node preference (prefer assembly directly attached to target_tid over strain nodes)
    species_node_pool = [r for r in pool if int(r["t_id"]) == target_tid]
    if species_node_pool:
        pool = species_node_pool
        reasons.append("species_node_preference")
        if len(pool) == 1:
            return pool[0], "+".join(reasons)

    # Determine N50 and sequence count columns based on current assembly level
    current_level = pool[0].get("a_level", "").strip()
    if current_level == "Contig":
        n50_col = "a_contig_n50"
        num_seq_col = "a_num_contigs"
    else:
        n50_col = "a_scaffold_n50"
        num_seq_col = "a_num_scaffolds"

    # Step 4: Highest N50
    n50_vals = [parse_float_safe(r.get(n50_col, "")) for r in pool]
    valid_n50s = [v for v in n50_vals if v is not None]
    if valid_n50s:
        max_n50 = max(valid_n50s)
        pool = [r for r, v in zip(pool, n50_vals) if v == max_n50]
        reasons.append(n50_col)
        if len(pool) == 1:
            return pool[0], "+".join(reasons)

    # Step 5: Lowest number of contigs/scaffolds
    seq_vals = [parse_float_safe(r.get(num_seq_col, "")) for r in pool]
    valid_seqs = [v for v in seq_vals if v is not None]
    if valid_seqs:
        min_seqs = min(valid_seqs)
        pool = [r for r, v in zip(pool, seq_vals) if v == min_seqs]
        reasons.append(num_seq_col)
        if len(pool) == 1:
            return pool[0], "+".join(reasons)

    # Step 6: Most recent release date
    date_vals = [parse_date_safe(r.get("a_release_date", "")) for r in pool]
    valid_dates = [d for d in date_vals if d is not None]
    if valid_dates:
        max_date = max(valid_dates)
        pool = [r for r, d in zip(pool, date_vals) if d == max_date]
        reasons.append("release_date")
        if len(pool) == 1:
            return pool[0], "+".join(reasons)

    # Step 7: Deterministic seeded tie-breaker (sort by accession first for full reproducibility)
    pool.sort(key=lambda r: r.get("a_accession", ""))
    rng = random.Random(random_seed + target_tid)
    chosen = rng.choice(pool)
    reasons.append(f"random_seed_{random_seed}")
    return chosen, "+".join(reasons)


def load_targets(targets_tsv: Path, logger) -> List[Dict[str, Any]]:
    """Load user-supplied target taxa from `config/targets.tsv`."""
    if not targets_tsv.is_file():
        raise FileNotFoundError(f"Targets TSV not found: {targets_tsv}")

    targets: List[Dict[str, Any]] = []
    seen_tids: Set[int] = set()
    with open(targets_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(
            (line for line in fh if line.strip() and not line.lstrip().startswith("#")),
            delimiter="\t",
        )
        if not reader.fieldnames or "t_id" not in reader.fieldnames:
            raise ValueError(f"Targets file {targets_tsv} must contain a 't_id' column.")

        for row in reader:
            raw_tid = (row.get("t_id") or "").strip()
            if not raw_tid:
                continue
            tid = int(raw_tid)
            if tid in seen_tids:
                logger.warning(f"Skipping duplicate target TaxID {tid} in {targets_tsv}")
                continue
            seen_tids.add(tid)
            species_label = (row.get("t_species") or "").strip()
            override_acc = (row.get("override_accession") or "").strip()
            if override_acc.upper() in {"NA", "NONE", "NULL", "-"}:
                override_acc = ""
            if override_acc:
                logger.info(f"TaxID {tid} ({species_label}): manual override_accession='{override_acc}' specified in {targets_tsv.name}")
            targets.append({
                "t_id": tid,
                "t_species_input": species_label,
                "override_accession": override_acc,
            })

    logger.info(f"Loaded {len(targets)} unique target taxa from {targets_tsv}")
    return targets


def scan_canonical_counts(
    counts_tsv_gz: Path,
    clade_map: Dict[int, Set[int]],
    logger,
) -> Tuple[List[str], Dict[int, Dict[str, int]], Dict[int, Dict[str, str]]]:
    """
    Scan `lj_stringent_canonical_counts.tsv.gz` once to extract:
      - Ordered list of the 516 weekly column names (`YYYY_WW`)
      - Per-target weekly read counts summed across any matching rows in the target's clade
        (Note: for species-level targets, `lj_stringent_canonical_counts.tsv.gz` is already
        rolled up to canonical ranks so the species row itself holds the exact clade count).
      - Taxonomic metadata (`t_domain`, `t_phylum`, `t_class`, `t_order`, `t_family`, `t_genus`,
        `t_common_name`, `db_mm_clade`, `db_mm_canonical_taxon`) from the canonical table.
    """
    if not counts_tsv_gz.is_file():
        raise FileNotFoundError(f"Canonical counts file not found: {counts_tsv_gz}")

    # Map every descendant TaxID -> list of target root TaxIDs whose clade contains it
    desc_to_roots: Dict[int, List[int]] = {}
    for root_tid, cset in clade_map.items():
        for dtid in cset:
            desc_to_roots.setdefault(dtid, []).append(root_tid)

    open_fn = gzip.open if str(counts_tsv_gz).endswith(".gz") else open
    weekly_counts: Dict[int, Dict[str, int]] = {}
    tax_meta: Dict[int, Dict[str, str]] = {}
    week_cols: List[str] = []

    logger.info(f"Scanning canonical read counts from {counts_tsv_gz} ...")
    with open_fn(counts_tsv_gz, "rt", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        col_idx = {name: i for i, name in enumerate(header)}

        week_indices: List[Tuple[str, int]] = [
            (col, idx)
            for idx, col in enumerate(header)
            if len(col) == 7 and col[:4].isdigit() and col[4] == "_" and col[5:].isdigit()
        ]
        week_cols = [w[0] for w in week_indices]
        logger.info(f"Identified {len(week_cols)} weekly timeseries columns ({week_cols[0]} .. {week_cols[-1]})")

        for root_tid in clade_map:
            weekly_counts[root_tid] = {w: 0 for w in week_cols}

        for row in reader:
            if not row or not row[0].strip():
                continue
            try:
                row_tid = int(row[0].strip())
            except ValueError:
                continue

            if row_tid not in desc_to_roots:
                continue

            # Extract weekly integers
            for root_tid in desc_to_roots[row_tid]:
                wdict = weekly_counts[root_tid]
                for wcol, widx in week_indices:
                    val_str = row[widx].strip() if widx < len(row) else ""
                    if val_str and val_str != "0" and val_str.upper() != "NA":
                        wdict[wcol] += int(float(val_str))

                # Save taxonomic lineage metadata if this row is the target species itself
                if row_tid == root_tid:
                    meta_fields = [
                        "t_macro_group", "t_domain", "t_kingdom", "t_phylum",
                        "t_class", "t_order", "t_family", "t_genus", "t_species",
                        "t_scientific_name", "t_common_name", "t_rank",
                        "db_mm_clade", "db_mm_canonical_taxon",
                    ]
                    tax_meta[root_tid] = {
                        k: (row[col_idx[k]].strip() if k in col_idx and col_idx[k] < len(row) else "")
                        for k in meta_fields
                    }

    return week_cols, weekly_counts, tax_meta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 0: Resolve JolTax clades, enforce read-count & organelle gates, and select reference assemblies."
    )
    parser.add_argument("--targets", required=True, help="Path to config/targets.tsv")
    parser.add_argument("--manifest", required=True, help="Path to database_manifest.tsv")
    parser.add_argument("--joltax-cache", required=True, help="Path to JolTax cache directory")
    parser.add_argument("--canonical-counts", required=True, help="Path to lj_stringent_canonical_counts.tsv.gz")
    parser.add_argument("--min-total-reads", type=int, default=1000, help="Minimum total clade reads across 516 weeks")
    parser.add_argument("--force-include-taxids", nargs="*", type=int, default=[], help="TaxIDs exempt from min_total_reads gate")
    parser.add_argument("--random-seed", type=int, default=42, help="Deterministic tie-breaker seed")
    parser.add_argument("--out-ledger", required=True, help="Output TSV: stage0_preflight_and_assemblies.tsv")
    parser.add_argument("--out-accessions", required=True, help="Output TXT: selected_accessions.txt for NCBI datasets")
    parser.add_argument("--out-clade-map", required=True, help="Output TSV: descendant_taxid -> target_t_id mapping for Stage 0Passed species")
    parser.add_argument("--out-weekly-dir", required=True, help="Output directory for per-species and combined 516-week raw read counts")
    parser.add_argument("--local-assemblies-dir", default="input/local_assemblies", help="Directory containing locally-provided genome FASTAs")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger("resolve_clades_and_assemblies", args.log)
    force_set: Set[int] = set(args.force_include_taxids or [])
    local_assemblies_dir = Path(args.local_assemblies_dir) if args.local_assemblies_dir else None

    # 1. Load target species
    targets = load_targets(Path(args.targets), logger)

    # 2. Load JolTax taxonomy tree and expand clades
    logger.info(f"Loading JolTax taxonomy cache from {args.joltax_cache} ...")
    tree = joltax.JolTree.load(args.joltax_cache)

    clade_map: Dict[int, Set[int]] = {}
    joltax_info: Dict[int, Dict[str, Any]] = {}
    for item in targets:
        tid = item["t_id"]
        try:
            clade_list = tree.get_clade(tid)
            sci_name = tree.get_name(tid) or item["t_species_input"] or f"TaxID_{tid}"
            com_name = tree.get_common_name(tid) or ""
            rank = tree.get_rank(tid) or ""
            clade_map[tid] = set(int(x) for x in clade_list)
            joltax_info[tid] = {
                "found": True,
                "scientific_name": sci_name,
                "common_name": com_name,
                "rank": rank,
                "clade_list": [int(x) for x in clade_list],
            }
            logger.info(f"TaxID {tid} ({sci_name}): resolved clade with {len(clade_list)} TaxID(s)")
        except Exception as exc:
            logger.warning(f"TaxID {tid} ({item['t_species_input']}) not found in JolTax: {exc}")
            clade_map[tid] = set()
            joltax_info[tid] = {
                "found": False,
                "scientific_name": item["t_species_input"] or f"TaxID_{tid}",
                "common_name": "",
                "rank": "",
                "clade_list": [],
            }

    # 3. Scan lj_stringent_canonical_counts.tsv.gz for Gate 1 & 516-week baseline time series
    week_cols, weekly_counts, tax_meta = scan_canonical_counts(
        Path(args.canonical_counts), clade_map, logger
    )

    # Write 516-week raw Kraken counts per species and as a combined long TSV
    weekly_dir = Path(args.out_weekly_dir)
    weekly_dir.mkdir(parents=True, exist_ok=True)
    combined_weekly_path = weekly_dir / "raw_weekly_counts_long.tsv"
    with open(combined_weekly_path, "w", encoding="utf-8", newline="") as comb_fh:
        comb_writer = csv.writer(comb_fh, delimiter="\t")
        comb_writer.writerow(["t_id", "t_species", "year_week", "year", "week", "raw_kraken_reads"])
        for item in targets:
            tid = item["t_id"]
            sp_name = joltax_info[tid]["scientific_name"]
            sp_weekly_path = weekly_dir / f"{tid}.raw_weekly.tsv"
            with open(sp_weekly_path, "w", encoding="utf-8", newline="") as sp_fh:
                sp_writer = csv.writer(sp_fh, delimiter="\t")
                sp_writer.writerow(["t_id", "year_week", "year", "week", "raw_kraken_reads"])
                for wcol in week_cols:
                    yr, wk = wcol.split("_")
                    cnt = weekly_counts.get(tid, {}).get(wcol, 0)
                    sp_writer.writerow([tid, wcol, int(yr), int(wk), cnt])
                    comb_writer.writerow([tid, sp_name, wcol, int(yr), int(wk), cnt])

    # 4. Scan database_manifest.tsv for all assemblies belonging to any target's clade
    all_clade_tids: Set[int] = set().union(*clade_map.values()) if clade_map else set()
    clade_assemblies: Dict[int, List[Dict[str, str]]] = {item["t_id"]: [] for item in targets}

    manifest_path = Path(args.manifest)
    if not manifest_path.is_file() and Path(str(manifest_path) + ".gz").is_file():
        manifest_path = Path(str(manifest_path) + ".gz")

    logger.info(f"Scanning database manifest {manifest_path} for {len(all_clade_tids)} clade TaxIDs ...")
    open_manifest = gzip.open if str(manifest_path).endswith(".gz") else open
    with open_manifest(manifest_path, "rt", encoding="utf-8", errors="replace") as mfh:
        mreader = csv.DictReader(mfh, delimiter="\t")
        for row in mreader:
            raw_t = (row.get("t_id") or "").strip()
            if not raw_t:
                continue
            try:
                rtid = int(raw_t)
            except ValueError:
                continue
            if rtid in all_clade_tids:
                for root_tid, cset in clade_map.items():
                    if rtid in cset:
                        clade_assemblies[root_tid].append(row)

    # 4b. Query live NCBI Datasets CLI for real-time assembly_status of all candidate nuclear assemblies
    # Why we query live NCBI rather than trusting `a_status` in `database_manifest.tsv`:
    # `database_manifest.tsv` is a static snapshot from when the Kraken 2 database was built.
    # Assemblies that were `current` at build time may have subsequently been suppressed or
    # removed at the submitter's request on ENA/NCBI (e.g. GCA_965250375.1 for Columba palumbus),
    # in which case NCBI Datasets omits their FASTA files. Filtering out suppressed candidates
    # BEFORE running the 4-tier selection cascade allows clades with multiple assemblies to
    # automatically fall back to another live assembly in the clade, or cleanly flag the clade
    # as EXCLUDED_NCBI_UNAVAILABLE if all candidates in the clade have been suppressed.
    candidate_clean_accs: Set[str] = set()
    for item in targets:
        tid = item["t_id"]
        ov_acc = item.get("override_accession", "")
        if ov_acc:
            local_fp = resolve_local_genome_file(ov_acc, local_assemblies_dir)
            if local_fp:
                ext_acc = extract_accession_from_override(ov_acc, local_fp)
                if is_valid_ncbi_accession(ext_acc):
                    candidate_clean_accs.add(ext_acc)
            else:
                c_ov = clean_ncbi_accession(ov_acc)
                if is_valid_ncbi_accession(c_ov):
                    candidate_clean_accs.add(c_ov)
                    if c_ov.startswith("GCF_"):
                        candidate_clean_accs.add("GCA_" + c_ov[4:])
        for r in clade_assemblies.get(tid, []):
            if not is_organellar_row(r) and r.get("a_status", "").strip().lower() == "current":
                raw_a = (r.get("a_accession") or "").strip()
                if raw_a:
                    c_a = clean_ncbi_accession(raw_a)
                    if is_valid_ncbi_accession(c_a):
                        candidate_clean_accs.add(c_a)
                        if c_a.startswith("GCF_"):
                            candidate_clean_accs.add("GCA_" + c_a[4:])

    live_ncbi_meta: Dict[str, Dict[str, Any]] = {}
    if candidate_clean_accs:
        datasets_bin = shutil.which("datasets")
        if not datasets_bin:
            # Fallback search in local .snakemake/conda environments if running before env rebuild
            for cand_bin in Path(".snakemake/conda").glob("*/bin/datasets"):
                if cand_bin.is_file() and os.access(cand_bin, os.X_OK):
                    datasets_bin = str(cand_bin)
                    break

        def _query_ncbi_batch(accs_to_query: Set[str]) -> None:
            if not accs_to_query or not datasets_bin:
                return
            with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as tf:
                for acc in sorted(accs_to_query):
                    tf.write(f"{acc}\n")
                tmp_acc_path = tf.name
            try:
                cmd = [
                    datasets_bin, "summary", "genome", "accession",
                    "--inputfile", tmp_acc_path,
                    "--as-json-lines",
                ]
                api_key = os.environ.get("NCBI_API_KEY", "").strip()
                if api_key:
                    cmd.extend(["--api-key", api_key])
                proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
                if proc.returncode == 0 and proc.stdout:
                    for line in proc.stdout.splitlines():
                        line = line.strip()
                        if not line:
                            continue
                        try:
                            rec = json.loads(line)
                            acc = rec.get("accession") or rec.get("current_accession") or ""
                            ainfo = rec.get("assembly_info", {})
                            astats = rec.get("assembly_stats", {})
                            org = rec.get("organism", {})
                            paired_info = ainfo.get("paired_assembly", {})
                            st = (ainfo.get("assembly_status") or "unknown").strip().lower()
                            rsn = (ainfo.get("suppression_reason") or "").strip()
                            paired_acc = (
                                paired_info.get("accession")
                                or rec.get("paired_accession")
                                or ""
                            ).strip()
                            paired_st = (paired_info.get("status") or "").strip().lower()
                            src_db_raw = (rec.get("source_database") or "").upper()
                            src_db = "REFSEQ" if "REFSEQ" in src_db_raw else ("GENBANK" if "GENBANK" in src_db_raw else "NCBI")
                            if acc:
                                live_ncbi_meta[acc] = {
                                    "status": st,
                                    "reason": rsn,
                                    "paired_accession": paired_acc,
                                    "paired_status": paired_st,
                                    "level": (ainfo.get("assembly_level") or "NA").strip(),
                                    "refseq_category": (ainfo.get("refseq_category") or "none").strip(),
                                    "seqlen": str(astats.get("total_sequence_length") or "NA"),
                                    "taxid": str(org.get("tax_id") or "NA"),
                                    "source_db": src_db,
                                }
                        except Exception:
                            pass
                else:
                    logger.warning(
                        f"`datasets summary genome accession` returned code {proc.returncode}: {proc.stderr.strip()}"
                    )
            finally:
                try:
                    os.unlink(tmp_acc_path)
                except OSError:
                    pass

        if datasets_bin:
            logger.info(
                f"Querying live NCBI assembly status via `{datasets_bin} summary genome accession` "
                f"for {len(candidate_clean_accs)} candidate/fallback accession(s)..."
            )
            _query_ncbi_batch(candidate_clean_accs)
            # If any suppressed GCF_ record listed a paired_accession with a different version suffix not yet queried, fetch it
            extra_paired: Set[str] = set()
            for acc, mdict in list(live_ncbi_meta.items()):
                p_acc = mdict.get("paired_accession", "")
                if p_acc and p_acc not in live_ncbi_meta:
                    extra_paired.add(p_acc)
            if extra_paired:
                _query_ncbi_batch(extra_paired)
            logger.info(
                f"Received live NCBI metadata for {len(live_ncbi_meta)} accession(s)."
            )
        else:
            logger.warning("`datasets` binary not found on PATH; skipping live NCBI pre-check in Stage 0.")

    def _find_live_gca_for_suppressed_gcf(gcf_acc: str) -> Optional[str]:
        """
        Given a suppressed RefSeq `GCF_` accession, check if a live (`current` or `previous`)
        GenBank `GCA_` equivalent exists on NCBI (checking both `paired_accession` and `GCA_` + suffix).
        """
        if not gcf_acc.startswith("GCF_"):
            return None
        rec = live_ncbi_meta.get(gcf_acc, {})
        candidates_gca: List[str] = []
        if rec.get("paired_accession", "").startswith("GCA_"):
            candidates_gca.append(rec["paired_accession"])
        direct_gca = "GCA_" + gcf_acc[4:]
        if direct_gca not in candidates_gca:
            candidates_gca.append(direct_gca)
        for gca_cand in candidates_gca:
            g_meta = live_ncbi_meta.get(gca_cand)
            if g_meta and g_meta["status"] not in {"suppressed", "withdrawn", "retired", "missing_on_ncbi"}:
                return gca_cand
        return None

    # 5. Evaluate Pre-Flight Gate 1, Manual Override, Gate 2, Live NCBI (+ GCF->GCA Fallback), and Selection Cascade
    Path(args.out_ledger).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_accessions).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_clade_map).parent.mkdir(parents=True, exist_ok=True)

    ledger_rows: List[Dict[str, Any]] = []
    selected_ncbi_accessions: List[str] = []
    seen_accessions: Set[str] = set()
    passing_clade_pairs: List[Tuple[int, int, str]] = []

    for item in targets:
        tid = item["t_id"]
        override_acc = item.get("override_accession", "")
        info = joltax_info[tid]
        sci_name = info["scientific_name"]
        species_slug = sanitize_species_slug(item["t_species_input"] or sci_name)
        meta = tax_meta.get(tid, {})
        common_name = info["common_name"] or meta.get("t_common_name", "")
        total_reads = sum(weekly_counts.get(tid, {}).values())
        n_weeks_present = sum(1 for v in weekly_counts.get(tid, {}).values() if v > 0)

        all_rows = clade_assemblies.get(tid, [])
        org_rows = [r for r in all_rows if is_organellar_row(r)]
        nuclear_manifest_rows = [
            r for r in all_rows
            if not is_organellar_row(r) and r.get("a_status", "").strip().lower() == "current"
        ]

        # Partition nuclear_manifest_rows into live-available (with GCF->GCA fallback) vs. suppressed on NCBI
        nuclear_current_rows: List[Dict[str, str]] = []
        suppressed_details: List[str] = []
        for r in nuclear_manifest_rows:
            c_acc = clean_ncbi_accession(r.get("a_accession", ""))
            if live_ncbi_meta:
                m_entry = live_ncbi_meta.get(c_acc)
                ncbi_st = m_entry["status"] if m_entry else "missing_on_ncbi"
                ncbi_rsn = m_entry["reason"] if m_entry else "not returned by NCBI datasets"
                if ncbi_st in {"suppressed", "withdrawn", "retired", "missing_on_ncbi"}:
                    # If this is a suppressed RefSeq (GCF_) genome, check if its GenBank (GCA_) equivalent is live
                    gca_equiv = _find_live_gca_for_suppressed_gcf(c_acc)
                    if gca_equiv:
                        gca_meta = live_ncbi_meta[gca_equiv]
                        r_copy = dict(r)
                        r_copy["a_accession"] = gca_equiv
                        r_copy["a_source_db"] = "GENBANK"
                        r_copy["_fallback_note"] = f"gcf_to_gca_fallback({c_acc}->{gca_equiv})"
                        nuclear_current_rows.append(r_copy)
                        logger.info(
                            f"TaxID {tid} ({sci_name}): suppressed RefSeq candidate {c_acc} ({ncbi_st}) "
                            f"recovered via live GenBank equivalent {gca_equiv} ({gca_meta['status']})"
                        )
                        continue
                    detail_str = f"{c_acc} [status={ncbi_st}" + (f": {ncbi_rsn}]" if ncbi_rsn else "]")
                    suppressed_details.append(detail_str)
                    logger.warning(f"TaxID {tid} ({sci_name}): excluding candidate assembly {detail_str}")
                    continue
            nuclear_current_rows.append(r)

        # Default empty selection fields
        chosen_acc = "NA"
        clean_acc = "NA"
        chosen_taxid = "NA"
        chosen_source = "NA"
        chosen_level = "NA"
        chosen_refseq = "NA"
        chosen_seqlen = "NA"
        selection_reason = "NA"

        # Evaluate gates in order:
        if not info["found"]:
            status = ValidationStatus.EXCLUDED_TAXID_NOT_FOUND
            reason = f"TaxID {tid} not present in JolTax taxonomy tree"
        elif total_reads < args.min_total_reads and tid not in force_set:
            status = ValidationStatus.EXCLUDED_LOW_READ_COUNT
            reason = (
                f"Total clade reads across 516 weeks ({total_reads:,}) is below "
                f"min_total_reads threshold ({args.min_total_reads:,})"
            )
        elif override_acc:
            # Manual user override from config/targets.tsv (either a local genome filename/path or an NCBI accession)
            local_fp = resolve_local_genome_file(override_acc, local_assemblies_dir)
            if local_fp:
                resolved_ov_acc = extract_accession_from_override(override_acc, local_fp)
                ov_sel_reason = f"local_genome_file({local_fp.name})"
                ov_meta = None
                ov_st = "local"
                ov_rsn = f"provided locally at {local_fp}"
                logger.info(
                    f"TaxID {tid} ({sci_name}): manual override specifies local genome file {local_fp} "
                    f"(accession identifier: {resolved_ov_acc})"
                )
            elif not is_valid_ncbi_accession(override_acc):
                resolved_ov_acc = ""
                ov_meta = None
                ov_st = "local_file_not_found"
                ov_rsn = f"file not found directly or under {local_assemblies_dir}"
                ov_sel_reason = "NA"
            else:
                c_ov = clean_ncbi_accession(override_acc)
                ov_meta = live_ncbi_meta.get(c_ov) if live_ncbi_meta else None
                ov_st = ov_meta["status"] if ov_meta else ("unknown" if not live_ncbi_meta else "missing_on_ncbi")
                ov_rsn = ov_meta["reason"] if ov_meta else "not returned by NCBI datasets"

                resolved_ov_acc = c_ov
                ov_sel_reason = "manual_override"

                if live_ncbi_meta and ov_st in {"suppressed", "withdrawn", "retired", "missing_on_ncbi"}:
                    gca_equiv = _find_live_gca_for_suppressed_gcf(c_ov)
                    if gca_equiv:
                        resolved_ov_acc = gca_equiv
                        ov_meta = live_ncbi_meta[gca_equiv]
                        ov_st = ov_meta["status"]
                        ov_sel_reason = f"manual_override+gcf_to_gca_fallback({c_ov}->{gca_equiv})"
                        logger.info(
                            f"TaxID {tid} ({sci_name}): manual override {c_ov} is suppressed on RefSeq; "
                            f"falling back to live GenBank equivalent {gca_equiv}"
                        )
                    else:
                        resolved_ov_acc = ""

            if not resolved_ov_acc:
                status = ValidationStatus.EXCLUDED_NCBI_UNAVAILABLE
                if ov_st == "local_file_not_found":
                    reason = (
                        f"Local override genome file '{override_acc}' was not found on disk "
                        f"({ov_rsn})"
                    )
                else:
                    reason = (
                        f"Manual override assembly {override_acc} is unavailable on live NCBI "
                        f"[status={ov_st}: {ov_rsn}] and no live GCA equivalent exists"
                    )
            else:
                # Also check if manifest has metadata for this accession as fallback
                manifest_match = next(
                    (r for r in all_rows if clean_ncbi_accession(r.get("a_accession", "")) == resolved_ov_acc),
                    None,
                )
                chosen_acc = resolved_ov_acc
                clean_acc = resolved_ov_acc
                chosen_taxid = (
                    (ov_meta["taxid"] if ov_meta and ov_meta["taxid"] != "NA" else None)
                    or (manifest_match.get("t_id", "").strip() if manifest_match else None)
                    or str(tid)
                )
                chosen_source = (
                    (ov_meta["source_db"] if ov_meta else None)
                    or (manifest_match.get("a_source_db", "").strip() if manifest_match else None)
                    or ("LOCAL" if local_fp else ("REFSEQ" if resolved_ov_acc.startswith("GCF_") else "GENBANK"))
                )
                chosen_level = (
                    (ov_meta["level"] if ov_meta and ov_meta["level"] != "NA" else None)
                    or (manifest_match.get("a_level", "").strip() if manifest_match else None)
                    or ("Local" if local_fp else "NA")
                )
                chosen_refseq = (
                    (ov_meta["refseq_category"] if ov_meta and ov_meta["refseq_category"] != "none" else None)
                    or (manifest_match.get("a_refseq_category", "").strip() if manifest_match else None)
                    or "none"
                )
                chosen_seqlen = (
                    (ov_meta["seqlen"] if ov_meta and ov_meta["seqlen"] != "NA" else None)
                    or (manifest_match.get("a_seqlen", "").strip() if manifest_match else None)
                    or (str(local_fp.stat().st_size) if local_fp else "NA")
                )
                selection_reason = ov_sel_reason
                status = ValidationStatus.PASSED_STAGE0
                forced_note = " (force_included)" if (total_reads < args.min_total_reads and tid in force_set) else ""
                reason = (
                    f"Passed Stage 0 ({total_reads:,} clade reads across {n_weeks_present}/516 weeks{forced_note}; "
                    f"selected {chosen_acc} [{chosen_level}, {chosen_refseq}] via {selection_reason})"
                )
                if clean_acc not in seen_accessions:
                    seen_accessions.add(clean_acc)
                    selected_ncbi_accessions.append(clean_acc)
                for dtid in info["clade_list"]:
                    passing_clade_pairs.append((dtid, tid, species_slug))
        elif not all_rows:
            status = ValidationStatus.EXCLUDED_NO_ASSEMBLY_IN_DB
            reason = f"No assemblies found in database_manifest.tsv for clade of TaxID {tid} ({len(info['clade_list'])} taxa)"
        elif not nuclear_manifest_rows and org_rows:
            status = ValidationStatus.EXCLUDED_ORGANELLAR_ONLY
            org_summary = ", ".join(
                f"{r.get('a_accession')} ({r.get('a_organelle_type') or 'organelle'})" for r in org_rows[:5]
            )
            reason = f"Clade contains only organellar assembly/assemblies ({len(org_rows)} total: {org_summary}) and no nuclear genome"
        elif not nuclear_manifest_rows:
            status = ValidationStatus.EXCLUDED_NO_ASSEMBLY_IN_DB
            reason = f"Clade has {len(all_rows)} manifest row(s) but 0 with a_status=='current' and non-organellar"
        elif not nuclear_current_rows and suppressed_details:
            status = ValidationStatus.EXCLUDED_NCBI_UNAVAILABLE
            reason = (
                f"All {len(nuclear_manifest_rows)} candidate nuclear assembly(ies) in clade are "
                f"suppressed/unavailable on live NCBI (and no live GCA equivalent): {'; '.join(suppressed_details)}"
            )
        else:
            # Passed Pre-Flight Gate 1, Gate 2, and live NCBI availability -> Select single best assembly!
            chosen_row, sel_reason = select_best_assembly_for_clade(
                tid, nuclear_current_rows, random_seed=args.random_seed
            )
            fb_note = chosen_row.get("_fallback_note", "")
            if fb_note:
                sel_reason = f"{sel_reason}+{fb_note}"
            chosen_acc = chosen_row.get("a_accession", "").strip()
            clean_acc = clean_ncbi_accession(chosen_acc)
            chosen_taxid = chosen_row.get("t_id", "").strip()
            chosen_source = chosen_row.get("a_source_db", "").strip()
            chosen_level = chosen_row.get("a_level", "").strip()
            chosen_refseq = chosen_row.get("a_refseq_category", "").strip() or "none"
            chosen_seqlen = chosen_row.get("a_seqlen", "").strip()
            selection_reason = sel_reason

            status = ValidationStatus.PASSED_STAGE0
            forced_note = " (force_included)" if (total_reads < args.min_total_reads and tid in force_set) else ""
            supp_note = f"; skipped {len(suppressed_details)} suppressed candidate(s)" if suppressed_details else ""
            reason = (
                f"Passed Stage 0 ({total_reads:,} clade reads across {n_weeks_present}/516 weeks{forced_note}; "
                f"selected {chosen_acc} [{chosen_level}, {chosen_refseq}] via {sel_reason}{supp_note})"
            )

            if clean_acc not in seen_accessions:
                seen_accessions.add(clean_acc)
                selected_ncbi_accessions.append(clean_acc)

            for dtid in info["clade_list"]:
                passing_clade_pairs.append((dtid, tid, species_slug))

        logger.info(f"[{status}] TaxID {tid} ({sci_name}): {reason}")

        ledger_rows.append({
            "t_id": tid,
            "t_species": item["t_species_input"] or sci_name,
            "scientific_name": sci_name,
            "common_name": common_name,
            "species_slug": species_slug,
            "t_domain": meta.get("t_domain", ""),
            "t_phylum": meta.get("t_phylum", ""),
            "t_class": meta.get("t_class", ""),
            "t_order": meta.get("t_order", ""),
            "t_family": meta.get("t_family", ""),
            "t_genus": meta.get("t_genus", ""),
            "clade_size": len(info["clade_list"]),
            "clade_taxids_csv": ",".join(str(x) for x in info["clade_list"]),
            "total_clade_reads": total_reads,
            "weeks_present": n_weeks_present,
            "db_mm_clade": meta.get("db_mm_clade", ""),
            "db_mm_canonical_taxon": meta.get("db_mm_canonical_taxon", ""),
            "candidate_nuclear_assemblies": len(nuclear_current_rows),
            "organellar_assemblies": len(org_rows),
            "chosen_accession": chosen_acc,
            "clean_ncbi_accession": clean_acc,
            "chosen_assembly_taxid": chosen_taxid,
            "chosen_source_db": chosen_source,
            "chosen_level": chosen_level,
            "chosen_refseq_category": chosen_refseq,
            "manifest_seqlen_bp": chosen_seqlen,
            "selection_reason": selection_reason,
            "stage0_status": status,
            "stage0_reason": reason,
        })

    # Write stage0_preflight_and_assemblies.tsv
    fieldnames = list(ledger_rows[0].keys()) if ledger_rows else [
        "t_id", "t_species", "scientific_name", "common_name", "species_slug",
        "t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus",
        "clade_size", "clade_taxids_csv", "total_clade_reads", "weeks_present",
        "db_mm_clade", "db_mm_canonical_taxon", "candidate_nuclear_assemblies",
        "organellar_assemblies", "chosen_accession", "clean_ncbi_accession",
        "chosen_assembly_taxid", "chosen_source_db", "chosen_level",
        "chosen_refseq_category", "manifest_seqlen_bp", "selection_reason",
        "stage0_status", "stage0_reason",
    ]
    with open(args.out_ledger, "w", encoding="utf-8", newline="") as lfh:
        writer = csv.DictWriter(lfh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(ledger_rows)

    # Write selected_accessions.txt for NCBI datasets bulk dehydrated download
    with open(args.out_accessions, "w", encoding="utf-8") as afh:
        for acc in selected_ncbi_accessions:
            afh.write(f"{acc}\n")

    # Write descendant_taxid -> target_t_id clade map for fast lookup in Stage 2
    with open(args.out_clade_map, "w", encoding="utf-8", newline="") as cfh:
        cwriter = csv.writer(cfh, delimiter="\t")
        cwriter.writerow(["descendant_taxid", "target_t_id", "species_slug"])
        for dtid, rtid, slug in passing_clade_pairs:
            cwriter.writerow([dtid, rtid, slug])

    n_pass = sum(1 for r in ledger_rows if r["stage0_status"] == ValidationStatus.PASSED_STAGE0)
    logger.info(
        f"Stage 0 complete: {n_pass}/{len(ledger_rows)} species passed Pre-Flight Gates "
        f"({len(selected_ncbi_accessions)} unique NCBI assembly accessions queued for bulk download)."
    )


if __name__ == "__main__":
    main()
