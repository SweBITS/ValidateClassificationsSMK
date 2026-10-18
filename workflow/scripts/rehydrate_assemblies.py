#!/usr/bin/env python3
"""
rehydrate_assemblies.py — Bulk NCBI Datasets Dehydrated Download, Rehydration & Catalog Parser
==============================================================================================

WHY THIS SCRIPT EXISTS:
-----------------------
Downloading dozens or hundreds of eukaryotic and prokaryotic reference genomes individually
via FTP or HTTP is brittle on HPC clusters due to transient network drops and rate limits.
Instead, we use the two-step NCBI `datasets` dehydrated workflow:
  1. `datasets download genome accession --dehydrated --inputfile selected_accessions.txt`
     fetches a lightweight metadata archive (`dehydrated.zip`) containing `fetch.txt` URLs
     and `md5sum.txt` checksums for every requested accession.
  2. `datasets rehydrate` downloads the actual genomic FASTA payloads in parallel (up to 30
     concurrent workers).
  3. After each rehydration pass, this script verifies the MD5 checksum of every downloaded
     file in parallel. Any truncated or corrupted file is unlinked and selectively
     re-fetched on the next retry iteration.
  4. Finally, `dataset_catalog.json` is parsed to map each requested NCBI accession to its
     primary `GENOMIC_NUCLEOTIDE_FASTA` file on disk (`accession_to_fasta.tsv`).
"""

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_ledger import extract_accession_from_override, resolve_local_genome_file, setup_logger


def run_cmd(cmd: List[str], logger, cwd: Optional[Path] = None) -> bool:
    """Execute a subprocess command, logging stdout/stderr and returning True on exit code 0."""
    # Redact API key in logs if present
    logged_cmd = [
        "***REDACTED***" if (i > 0 and cmd[i - 1] == "--api-key") else tok
        for i, tok in enumerate(cmd)
    ]
    logger.info(f"Running command: {' '.join(logged_cmd)}")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=cwd)
        if proc.stdout and proc.stdout.strip():
            logger.info(f"stdout:\n{proc.stdout.strip()}")
        if proc.stderr and proc.stderr.strip():
            logger.info(f"stderr:\n{proc.stderr.strip()}")
        proc.check_returncode()
        return True
    except subprocess.CalledProcessError as exc:
        logger.error(f"Command failed with exit code {exc.returncode}: {' '.join(logged_cmd)}")
        return False
    except Exception as exc:
        logger.error(f"Unexpected error executing {' '.join(logged_cmd)}: {exc}")
        return False


def compute_md5_worker(filepath_str: str) -> str:
    """Calculate the hexadecimal MD5 digest of a file, or return an error sentinel."""
    fp = Path(filepath_str)
    if not fp.is_file():
        return "file_not_found"
    md5 = hashlib.md5()
    try:
        with open(fp, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                md5.update(chunk)
        return md5.hexdigest().lower()
    except OSError:
        return "io_error"


def parse_md5sum_txt(package_dir: Path, logger) -> Optional[Dict[str, str]]:
    """Parse `md5sum.txt` inside the unzipped NCBI dataset directory."""
    md5_path = package_dir / "md5sum.txt"
    if not md5_path.is_file():
        logger.error(f"Missing md5sum.txt in {package_dir}")
        return None

    expected: Dict[str, str] = {}
    with open(md5_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                digest, rel_path = parts
                rel_clean = rel_path.lstrip("./")
                abs_fp = str((package_dir / rel_clean).resolve())
                expected[abs_fp] = digest.lower()
    logger.info(f"Loaded {len(expected)} expected file checksums from {md5_path}")
    return expected


def rehydrate_with_md5_retries(
    package_dir: Path,
    max_attempts: int,
    threads: int,
    api_key: Optional[str],
    logger,
) -> bool:
    """Run `datasets rehydrate` followed by parallel MD5 verification with automatic retry."""
    expected_map = parse_md5sum_txt(package_dir, logger)
    if expected_map is None:
        return False
    if not expected_map:
        logger.info("md5sum.txt is empty; nothing to rehydrate.")
        return True

    critical_files: Set[str] = {
        str((package_dir / "ncbi_dataset" / "data" / "assembly_data_report.jsonl").resolve()),
        str((package_dir / "ncbi_dataset" / "data" / "dataset_catalog.json").resolve()),
    }
    pending_files: Set[str] = set(expected_map.keys())
    rehydrate_workers = max(1, min(threads, 30))

    for attempt in range(1, max_attempts + 1):
        logger.info(f"=== Rehydration Attempt {attempt}/{max_attempts} ===")
        cmd = [
            "datasets",
            "rehydrate",
            "--directory",
            str(package_dir),
            "--max-workers",
            str(rehydrate_workers),
            "--no-progressbar",
        ]
        if api_key:
            cmd.extend(["--api-key", api_key])

        if not run_cmd(cmd, logger):
            logger.warning(f"datasets rehydrate returned non-zero on attempt {attempt}; checking downloaded files...")

        sorted_pending = sorted(pending_files)
        with ProcessPoolExecutor(max_workers=max(1, threads)) as pool:
            actual_digests = list(pool.map(compute_md5_worker, sorted_pending))

        failed_this_round: Set[str] = set()
        for fp_str, actual_md5 in zip(sorted_pending, actual_digests):
            exp_md5 = expected_map[fp_str]
            if actual_md5 != exp_md5:
                logger.warning(f"MD5 check failed for {fp_str}: expected={exp_md5}, got={actual_md5}")
                if fp_str in critical_files:
                    logger.error(f"Critical metadata file corrupted: {fp_str}")
                    return False
                failed_this_round.add(fp_str)

        if not failed_this_round:
            logger.info(f"All {len(pending_files)} file(s) passed MD5 verification on attempt {attempt}.")
            return True

        logger.warning(f"{len(failed_this_round)} file(s) failed MD5 verification on attempt {attempt}.")
        if attempt < max_attempts:
            for bad_fp in failed_this_round:
                p = Path(bad_fp)
                if p.exists():
                    p.unlink()
            pending_files = failed_this_round
            time.sleep(5)
        else:
            logger.error(f"Exhausted {max_attempts} rehydration attempts; {len(failed_this_round)} file(s) failed.")
            return False

    return False


def build_accession_to_fasta_manifest(
    package_dir: Path,
    out_manifest: Path,
    logger,
    local_assemblies: Optional[Dict[str, Tuple[Path, int]]] = None,
) -> None:
    """
    Parse `dataset_catalog.json` inside `package_dir/ncbi_dataset/data/` and merge with any
    locally-provided assemblies, writing a 3-column TSV:
    (`clean_ncbi_accession`, `fasta_path`, `uncompressed_bytes`).
    """
    catalog_path = package_dir / "ncbi_dataset" / "data" / "dataset_catalog.json"
    rows: List[Tuple[str, str, int]] = []
    seen_accs: Set[str] = set()

    if catalog_path.is_file():
        with open(catalog_path, "r", encoding="utf-8") as fh:
            catalog = json.load(fh)

        for item in catalog.get("assemblies", []):
            acc = item.get("accession")
            if not acc:
                continue
            genomic_files = [
                f for f in item.get("files", [])
                if f.get("fileType") == "GENOMIC_NUCLEOTIDE_FASTA"
            ]
            if not genomic_files:
                logger.warning(f"No GENOMIC_NUCLEOTIDE_FASTA listed in catalog for accession {acc}")
                continue
            genomic_files.sort(key=lambda f: int(f.get("uncompressedLengthBytes") or 0))
            largest = genomic_files[-1]
            rel_path = largest["filePath"]
            abs_fasta = (package_dir / "ncbi_dataset" / "data" / rel_path).resolve()
            size_bytes = int(largest.get("uncompressedLengthBytes") or (abs_fasta.stat().st_size if abs_fasta.is_file() else 0))
            if not abs_fasta.is_file():
                raise FileNotFoundError(f"Catalog FASTA file missing on disk for {acc}: {abs_fasta}")
            rows.append((acc, str(abs_fasta), size_bytes))
            seen_accs.add(acc)
    elif not local_assemblies:
        raise FileNotFoundError(f"dataset_catalog.json not found at {catalog_path}")

    # Merge locally provided reference assemblies
    if local_assemblies:
        for acc, (lfp, size_bytes) in local_assemblies.items():
            if acc not in seen_accs:
                rows.append((acc, str(lfp.resolve()), size_bytes))
                seen_accs.add(acc)
                logger.info(f"Added local assembly record: {acc} -> {lfp} ({size_bytes:,} bytes)")

    rows.sort(key=lambda r: r[0])
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    with open(out_manifest, "w", encoding="utf-8", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["clean_ncbi_accession", "fasta_path", "uncompressed_bytes"])
        for r in rows:
            writer.writerow(r)

    logger.info(f"Wrote {len(rows)} assembly FASTA paths to {out_manifest}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bulk download dehydrated NCBI genome packages, rehydrate with MD5 checks, and emit accession->FASTA manifest."
    )
    parser.add_argument("--accessions-file", required=True, help="Path to selected_accessions.txt")
    parser.add_argument("--targets", default=None, help="Optional path to config/targets.tsv for local file overrides")
    parser.add_argument("--work-dir", required=True, help="Directory for dehydrated.zip and download_dir")
    parser.add_argument("--out-manifest", required=True, help="Output TSV: accession_to_fasta.tsv")
    parser.add_argument("--out-flag", required=True, help="Output completion flag file")
    parser.add_argument("--local-assemblies-dir", default="input/local_assemblies", help="Directory containing locally-provided genome FASTAs")
    parser.add_argument("--max-attempts", type=int, default=5, help="Max rehydration attempts")
    parser.add_argument("--threads", type=int, default=16, help="Worker threads for rehydration & MD5 checks")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    logger = setup_logger("rehydrate_assemblies", args.log)
    import os
    api_key = (args.api_key or os.environ.get("NCBI_API_KEY", "")).strip() or None

    acc_file = Path(args.accessions_file)
    accessions = [
        line.strip() for line in acc_file.read_text(encoding="utf-8").splitlines() if line.strip()
    ] if acc_file.is_file() else []

    out_manifest = Path(args.out_manifest)
    out_flag = Path(args.out_flag)
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_flag.parent.mkdir(parents=True, exist_ok=True)

    # Graceful empty handling if 0 species passed Stage 0
    if not accessions:
        logger.warning("No accessions listed in selected_accessions.txt (0 species passed Stage 0). Creating empty manifest.")
        with open(out_manifest, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["clean_ncbi_accession", "fasta_path", "uncompressed_bytes"])
        out_flag.touch()
        return

    local_dir = Path(args.local_assemblies_dir) if args.local_assemblies_dir else None
    local_assemblies: Dict[str, Tuple[Path, int]] = {}
    accessions_set: Set[str] = set(accessions)

    # Check targets.tsv for any Stage 0 passing species explicitly configured with a local genome filename/path
    if args.targets and Path(args.targets).is_file():
        targets_file = Path(args.targets)
        with open(targets_file, "r", encoding="utf-8") as tfh:
            reader = csv.DictReader(
                (line for line in tfh if line.strip() and not line.lstrip().startswith("#")),
                delimiter="\t",
            )
            for row in reader:
                raw_ov = (row.get("override_accession") or "").strip()
                if raw_ov:
                    lfp = resolve_local_genome_file(raw_ov, local_dir)
                    if lfp:
                        clean_acc = extract_accession_from_override(raw_ov, lfp)
                        if clean_acc in accessions_set:
                            local_assemblies[clean_acc] = (lfp, lfp.stat().st_size)
                            logger.info(
                                f"Target TaxID {row.get('t_id')} ({row.get('t_species', '')}) specifies local genome file "
                                f"{lfp.name} ({lfp.stat().st_size:,} bytes). Skipping NCBI download for {clean_acc}."
                            )

    accs_to_download = [acc for acc in accessions if acc not in local_assemblies]

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    dehydrated_zip = work_dir / "dehydrated.zip"
    package_dir = work_dir / "download_dir"

    if accs_to_download:
        # Clean stale partial download dir if re-running from scratch
        if package_dir.exists():
            shutil.rmtree(package_dir)
        if dehydrated_zip.exists():
            dehydrated_zip.unlink()

        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as tf:
            for acc in accs_to_download:
                tf.write(f"{acc}\n")
            dl_acc_file = tf.name

        try:
            dl_cmd = [
                "datasets",
                "download",
                "genome",
                "accession",
                "--dehydrated",
                "--inputfile",
                dl_acc_file,
                "--no-progressbar",
                "--filename",
                str(dehydrated_zip),
            ]
            if api_key:
                dl_cmd.extend(["--api-key", api_key])

            logger.info(f"Downloading dehydrated NCBI package for {len(accs_to_download)} accession(s)...")
            if not run_cmd(dl_cmd, logger):
                sys.exit(1)

            logger.info(f"Extracting {dehydrated_zip} to {package_dir} ...")
            with zipfile.ZipFile(dehydrated_zip, "r") as zf:
                zf.extractall(package_dir)

            # Remove dehydrated.zip once extracted to save space
            if dehydrated_zip.exists():
                dehydrated_zip.unlink()

            if not rehydrate_with_md5_retries(
                package_dir=package_dir,
                max_attempts=args.max_attempts,
                threads=args.threads,
                api_key=api_key,
                logger=logger,
            ):
                sys.exit(1)
        finally:
            try:
                os.unlink(dl_acc_file)
            except OSError:
                pass
    else:
        logger.info(f"All {len(accessions)} accession(s) are locally provided. Skipping NCBI Datasets download.")
        package_dir.mkdir(parents=True, exist_ok=True)

    build_accession_to_fasta_manifest(package_dir, out_manifest, logger, local_assemblies=local_assemblies)
    out_flag.touch()
    logger.info("Bulk download and rehydration completed successfully.")


if __name__ == "__main__":
    main()
