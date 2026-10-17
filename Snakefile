from pathlib import Path
import pandas as pd
import yaml
import logging
import sys
import re

# TODO: Implement dynamic resource allocation that depend on the size of the input files.
# TODO: Implement resubmission of jobs that fail, that increases the memory and/or time resources based on the "attempts" variable
# TODO: Things for final report:
    # - For each tax_id and genome_basename:
        # - Number of reads to start with
        # - Mapping information (check the stats file from BBMap):
            # - Percentage reads mapped
            # - Percentage unambiguous
            # - Percentage ambiguous
            # - Percentage unmapped
        # - Percentage duplicates (removed)
        # - Number of generated reads and coverage
        # - Percentage correctly classified (to tax_id in mapping spec)


# ============================== #
# --- Housekeeping and Setup --- #
# ============================== #

# --- Configuration file --- #
configfile: "config.yaml"

# --- Define Core Paths and Patterns--- #
LOG_DIR = Path("logs")
OUTPUT_DIR = Path("output")
INPUT_DIR = Path("input")
GENOMES_DIR = Path("supporting_files/genomes")
ENVS_DIR = Path("envs")
MAPPING_SPEC = Path(config["MAPPING_SPECIFICATION"])
STARTING_FASTQ_PATTERN = str(INPUT_DIR / config["FASTQ_FILE_PATTERN"])

# --- Validate Kraken2 Path (conditinal on workflow mode in the main config) --- #
if config.get("RUN_VALIDATION", True) and config.get("KRAKEN2_BIN_DIR"):
    kraken2_bin_dir = Path(config["KRAKEN2_BIN_DIR"]).resolve()
    if not kraken2_bin_dir.is_dir():
        sys.exit(f"ERROR: KRAKEN2_BIN_DIR is not a valid directory: {kraken2_bin_dir}")
    kraken2_exe_path = kraken2_bin_dir / "kraken2"
    if not kraken2_exe_path.is_file() or not os.access(kraken2_exe_path, os.X_OK):
         sys.exit(f"ERROR: 'kraken2' not found or not executable in {kraken2_bin_dir}")
    print(f"WORKFLOW_INFO: Will use Kraken 2 binary in directory {kraken2_bin_dir}")


# ==================== #
# --- HELPER FUNCS --- #
# ==================== #

# --- Load Mapping Specification --- #
def load_mapping_spec(config_mapping_spec_path):
    """
    Load the mapping specification from the provided path.
    Expected columns: 'tax_id', 'name', and 'fasta_file'.
    Validates the name column format (alphanumerics and spaces only).
    Validates fasta filenames (must end in .fa/.fasta/.fna).
    Checks that each derived 'genome_basename' maps to only one 'fasta_file'.
    Removes duplicate rows (based on columns tax_id and fasta_file).
    Adds a 'genome_basename' column derived from 'fasta_file'.
    Returns a pandas DataFrame.
    """
    mapping_spec_file = Path(config_mapping_spec_path)
    required_columns = ['tax_id', 'name', 'fasta_file']
    allowed_suffixes = {'.fasta', '.fa', '.fna'}
    allowed_name_pattern = re.compile(r"^[a-zA-Z0-9 _-]+$")  # Allowed in the name col (letters, numbers, space, underscores, and hyphens)
    
    if not mapping_spec_file.is_file():
        sys.exit(f"WORKFLOW_ERROR: Mapping specification file not found at '{mapping_spec_file}'. Exiting.")

    try:
        mapping_df = pd.read_csv(
            mapping_spec_file,
            sep='\t',
            comment='#',
            dtype={'tax_id': str, 'name': str, 'fasta_file': str})

        # Validate that the mapping spec file is formatted correctly and is non-empty
        missing_cols = [col for col in required_columns if col not in mapping_df.columns]
        if missing_cols:
            sys.exit(f"WORKFLOW_ERROR: Mapping specification TSV '{mapping_spec_file}' is missing required columns: {missing_cols}. Exiting.")

        mapping_df.dropna(subset=required_columns, inplace=True)
        if mapping_df.empty:
            print(f"WORKFLOW_WARNING: Mapping specification file '{mapping_spec_file}' is empty or contains only comments/header.")
            return pd.DataFrame(columns=required_columns + ['genome_basename'])

        # Validate fasta_file suffixes
        invalid_files = []
        for fname in mapping_df['fasta_file']:
            if Path(fname).suffix.lower() not in allowed_suffixes:
                invalid_files.append(fname)
        if invalid_files:
            error_msg = (f"WORKFLOW_ERROR: Found rows in '{mapping_spec_file}' with invalid 'fasta_file' suffixes "
                         f"(must end in {', '.join(allowed_suffixes)}):\n" + "\n".join(invalid_files))
            sys.exit(error_msg)

        # Validate name column format
        invalid_names = []
        for index, row in mapping_df.iterrows():
            name = row['name']
            if not allowed_name_pattern.match(name):
                invalid_names.append(f"Row index {index}: tax_id='{row['tax_id']}', name='{name}'")
        if invalid_names:
            error_msg = (f"WORKFLOW_ERROR: Found rows in '{mapping_spec_file}' with invalid characters in the 'name' column "
                         f"(only letters, numbers, underscores, hyphens, and spaces are allowed):\n" + "\n".join(invalid_names))
            sys.exit(error_msg)

        # Check if listed fasta files actually exist
        missing_fasta_files = []
        for index, row in mapping_df.iterrows():
            fasta_filename = row['fasta_file']
            expected_path = GENOMES_DIR / fasta_filename
            if not expected_path.is_file():
                missing_fasta_files.append(f"Row index {index}: tax_id='{row['tax_id']}', file='{fasta_filename}' (expected at '{expected_path}')")
        if missing_fasta_files:
             error_msg = (f"WORKFLOW_ERROR: Found rows in '{mapping_spec_file}' referencing FASTA files that can't be found in '{GENOMES_DIR}':\n"
                          + "\n".join(missing_fasta_files))
             sys.exit(error_msg)

        # Drop duplicated rows
        initial_rows = len(mapping_df)
        mapping_df.drop_duplicates(['tax_id', 'fasta_file'], inplace=True)
        dropped_count = initial_rows - len(mapping_df)
        if dropped_count > 0:
            print(f"WORKFLOW_INFO: Removed {dropped_count} duplicate rows from mapping specification.")

        # Add the genome_basename column (remove suffix)
        mapping_df['genome_basename'] = mapping_df['fasta_file'].apply(lambda f: Path(f).stem)
        
        # Add the sanitized_name column (replace spaces with underscores)
        mapping_df['sanitized_name'] = mapping_df['name'].apply(lambda name: re.sub(r'\s+', '_', name))

        # Check for genome_basename mapping to multiple unique fasta_files
        basename_counts = mapping_df.groupby('genome_basename')['fasta_file'].nunique()
        conflicting_basenames = basename_counts[basename_counts > 1]
        if not conflicting_basenames.empty:
            error_msg = (f"WORKFLOW_ERROR: Inconsistency in mapping specification '{mapping_spec_file}'. "
                         f"Multiple files with the same basename (same name, different suffix) found.")
            sys.exit()

        print(f"WORKFLOW_INFO: Found {len(mapping_df)} unique and valid entries from '{mapping_spec_file}'.")

    except Exception as e:
        sys.exit(f"WORKFLOW_ERROR: An error occurred reading or processing mapping specification TSV '{mapping_spec_file}': {e}. Exiting.")

    return mapping_df

# --- Generate Target Output Files Function --- #
def get_target_outputs(mapping_spec_df, pattern_template):
    """
    Generates a list of output files based on the mapping specification DataFrame
    and the provided file name pattern template.
    """
    target_files = []
    for index, row in mapping_spec_df.iterrows():
        # Prepare format dictionary from the row
        format_dict = row.to_dict()

        # Generate the target path string using available keys from the row
        target_path_str = pattern_template.format(**format_dict)
        target_files.append(target_path_str)

    return target_files

# --- Workflow Target Function --- #
def get_workflow_targets():
    """
    Determines the final target file for the workflow based on the RUN_VALIDATION flag in the config.
    
    Target file is always at least the mapping branch target file.

    If validation:
        Target file becomes a flag file that triggers the coverage checkpoint to run, which determines 
        the combinations of tax_id/genome_basename values to validate.
    
    Also adds the target file for tracking start time of the workflow and the one to hold the
    snakemake version that's used for running the workflow.
    """
    final_targets = [WF_MAPPING_BRANCH_TARGET]
    
    if config.get("RUN_VALIDATION", True):
        print("RUN_VALIDATION=TRUE. Running the validation branch of the workflow.")
        final_targets.append(WF_VALIDATION_BRANCH_TARGET)
    
    final_targets.append(WF_TIMESTAMP_START_TARGET)
    final_targets.append(WF_SNAKEMAKE_VERSION_TARGET)

    return final_targets


# ================================== #
# --- WORKFLOWS AND TARGET FILES --- #
# ================================== #

wildcard_constraints:
    tax_id="[0-9]+"

# --- Final Targets --- #
WF_MAPPING_BRANCH_TARGET = str(OUTPUT_DIR / "1_mapping.done")
WF_VALIDATION_BRANCH_TARGET = str(OUTPUT_DIR / "2_validation.done")
WF_TIMESTAMP_START_TARGET = str(OUTPUT_DIR / "WF_start.txt")
WF_SNAKEMAKE_VERSION_TARGET = str(OUTPUT_DIR / "WF_snk_version.txt")

# --- Load the mapping specification --- #
MAPPING_SPEC_DF = load_mapping_spec(MAPPING_SPEC)

# --- Include Rule Modules --- #
include: "rules/map.smk"
include: "rules/simulation_validation.smk"

localrules: 
    all,
    WF_timestamp_start,
    WF_snakemake_version

# --- Define the final targets of the workflow --- #
rule all:
    input:
        get_workflow_targets()

# --- Rules creating workflow-related (WF) files for final report --- #
# Capture start time
rule WF_timestamp_start:
    output:
        ts_start = WF_TIMESTAMP_START_TARGET
    shell:
        "date --iso-8601=seconds > {output.ts_start}"

# Capture Snakemake version
rule WF_snakemake_version:
    output:
        txt_version = WF_SNAKEMAKE_VERSION_TARGET
    shell:
        "snakemake --version > {output.txt_version} 2>&1"