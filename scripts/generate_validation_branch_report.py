#!/usr/bin/env python3

from datetime import datetime
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from scipy import stats
import warnings
import logging
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import pysam
import json
import sys
import os
import re


# --- Snakemake objects ---
try:
    snk_input = snakemake.input
    snk_output = snakemake.output
    snk_params = snakemake.params
    snk_log = snakemake.log
    snk_config = snakemake.config
except AttributeError as e:
    raise AttributeError(f"Missing expected Snakemake object attribute: {e}")

# --- Setup Logging ---
log_file = snk_log.path
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=log_file,
    filemode='w',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

try:
    # Inputs
    ts_start_file = snk_input.ts_start  				# Start time of pipeline
    snk_version_file = snk_input.snk_version			# Snakemake version
    real_stats_files = snk_input.real_stats				# 
    markdup_logs = snk_input.markdup_logs				# Sambamba markdup stats
    dedup_bams = snk_input.dedup_bams					# Deduplicated bam files
    coverage_files = snk_input.coverage_files			# Coverage files
    sim_params_jsons = getattr(snk_input, 'sim_params_jsons', [])
    pirs_logs = getattr(snk_input, 'pirs_logs', [])
    slice_logs = getattr(snk_input, 'slice_logs', [])
    correct_ids_logs = snk_input.correct_ids_logs		# ??
    sim_stats_files = snk_input.sim_stats				# ??
    classifiable_beds = snk_input.classifiable_beds		# Classifiable region BED files
    summary_tsvs = snk_input.summary_tsvs				# Per contig: mapped reads and classifiable region length

    # Output
    report_html_file = snk_output.report_html			# The main output HTML file

    # Params
    mapping_df = snk_params.mapping_spec_df 			# The mapping spec DF 
    report_config_dict = snk_params.report_config		# Dict with Snakemake main config file parameter values
    path_patterns = snk_params.path_patterns 			# Dict with filename patterns

except Exception as e:
    logging.error(f"FATAL: Error accessing Snakemake objects: {e}")
    sys.exit(1)

else:
    logging.info("Successfully accessed Snakemake objects.")

# --- Helper Function for Safe File Reading ---
def read_file_content(filepath, log_msg='File'):
    """Safely reads the first line or all content of a small text file."""
    filepath = Path(filepath)
    try:
        if not filepath.exists():
            logging.error(f"{log_msg} not found: {filepath}")
            return "Error: File not found"
        content = filepath.read_text(encoding='utf-8').strip()
        if not content:
            logging.warning(f"{log_msg} is empty: {filepath}")
        return content
    except Exception as e:
        logging.error(f"Error reading {log_msg} file {filepath}: {e}")
    return f"Error reading file"

# --- Process General Workflow Info ---
logging.info("Gathering general workflow information.")
general_info = {}
try:
    general_info['start_time'] = read_file_content(ts_start_file, "Start time")
    general_info['end_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S %Z')  # Current time is considered to be the end time of the pipeline
    general_info['snakemake_version'] = read_file_content(snk_version_file, "Snakemake version").replace('\n', ' ')  # Ensure single line
    general_info['workdir'] = os.getcwd()
    general_info['config'] = report_config_dict
    general_info['mapping_spec_path'] = report_config_dict.get("mapping_spec_path", "N/A")

except Exception as e:
    logging.error(f"Error gathering general pipeline info: {e}")
    sys.exit(1)
else:
    logging.info("Done.")

# --- Parsing Functions ---
def parse_bbmap_stats(filepath):
    """
    Parses relevant fields from a BBMap stats file, including
    overall read usage, pairing data, and separate R1/R2 mapping stats.

    A BBmap stats file for reference:
    Reads Used:           	8029272	(955976917 bases)

    Mapping:          	1710.545 seconds.
    Reads/sec:       	4693.98
    kBases/sec:      	558.87


    Pairing data:   	pct pairs	num pairs 	pct bases	   num bases

    mated pairs:     	 68.5950% 	  2753838 	 69.6370% 	   665713372
    bad pairs:       	  0.0000% 	        0 	  0.0000% 	           0
    insert size avg: 	  149.20


    Read 1 data:      	pct reads	num reads 	pct bases	   num bases

    mapped:          	 88.2301% 	  3542116 	 88.1644% 	   421264344
    unambiguous:     	 68.5950% 	  2753838 	 69.6620% 	   332856686
    ambiguous:       	 19.6351% 	   788278 	 18.5024% 	    88407658
    low-Q discards:  	  0.0000% 	        0 	  0.0000% 	           0

    perfect best site:	 57.5986% 	  2312375 	 57.6869% 	   275637633
    semiperfect site:	 57.6250% 	  2313436 	 57.7187% 	   275789513
    rescued:         	  0.1153% 	     4630

    Match Rate:      	      NA 	       NA 	 99.6548% 	   332301046
    Error Rate:      	 12.0809% 	   439143 	  0.3442% 	     1147792
    Sub Rate:        	 11.4038% 	   414528 	  0.1612% 	      537365
    Del Rate:        	  0.5137% 	    18672 	  0.1786% 	      595425
    Ins Rate:        	  0.3279% 	    11920 	  0.0045% 	       15002
    N Rate:          	  0.0762% 	     2771 	  0.0010% 	        3273


    Read 2 data:      	pct reads	num reads 	pct bases	   num bases

    mapped:          	 88.0163% 	  3533535 	 87.8687% 	   420152947
    unambiguous:     	 68.5950% 	  2753838 	 69.6118% 	   332855658
    ambiguous:       	 19.4214% 	   779697 	 18.2569% 	    87297289
    low-Q discards:  	  0.0000% 	        1 	  0.0000% 	         151

    perfect best site:	 55.6268% 	  2233213 	 55.3962% 	   264882654
    semiperfect site:	 55.6510% 	  2234184 	 55.4253% 	   265021910
    rescued:         	  0.1113% 	     4467

    Match Rate:      	      NA 	       NA 	 99.6058% 	   332164629
    Error Rate:      	 14.3771% 	   516636 	  0.3925% 	     1308928
    Sub Rate:        	 13.7576% 	   494374 	  0.2013% 	      671302
    Del Rate:        	  0.5002% 	    17974 	  0.1870% 	      623503
    Ins Rate:        	  0.3128% 	    11242 	  0.0042% 	       14123
    N Rate:          	  0.1336% 	     4802 	  0.0017% 	        5604
    """
    logging.debug(f"Parsing BBMap stats: {filepath}")
    stats = {
        'input_reads_processed': 'N/A',
        'input_read_pairs': 'N/A',
        
        'mapped_reads_count': 'N/A',
        'mapped_reads_percent': 'N/A',
        'properly_paired_count': 'N/A',
        'properly_paired_percent': 'N/A',
        'singletons_count': 'N/A',
        'singletons_percent': 'N/A',

        'pairing_mated_pairs_count': 'N/A',
        'pairing_mated_pairs_percent': 'N/A',
        'pairing_bad_pairs_count': 'N/A',
        'pairing_bad_pairs_percent': 'N/A',
        'pairing_insert_size_avg': 'N/A',

        'r1_mapped_count': 'N/A', 
        'r1_mapped_percent': 'N/A',
        'r1_unambiguous_count': 'N/A', 
        'r1_unambiguous_percent': 'N/A',
        'r1_ambiguous_count': 'N/A', 
        'r1_ambiguous_percent': 'N/A',
        'r1_unmapped_count': 'N/A',
        'r1_unmapped_percent': 'N/A',
        
        'r2_mapped_count': 'N/A', 
        'r2_mapped_percent': 'N/A',
        'r2_unambiguous_count': 'N/A', 
        'r2_unambiguous_percent': 'N/A',
        'r2_ambiguous_count': 'N/A', 
        'r2_ambiguous_percent': 'N/A',
        'r2_unmapped_count': 'N/A',
        'r2_unmapped_percent': 'N/A'}

    if not filepath:
        logging.warning("Stats filepath was None. Returning N/A stats.")
        return stats

    current_section = None # Will be "pairing", "read1", or "read2"
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            all_lines = [l.strip() for l in f if l.strip()]

        # Handle samtools flagstat output (from minimap2 mapping)
        if any("QC-passed reads" in l for l in all_lines):
            for line in all_lines:
                if "in total (QC-passed" in line:
                    total = int(line.split()[0])
                    stats['input_reads_processed'] = total
                    stats['input_read_pairs'] = total // 2
                elif re.search(r"^\d+ \+ \d+ mapped \(", line):
                    mapped_cnt = int(line.split()[0])
                    try:
                        pct = float(line.split('(')[1].split('%')[0].strip())
                    except (IndexError, ValueError):
                        pct = 0.0
                    stats['mapped_reads_count'] = mapped_cnt
                    stats['mapped_reads_percent'] = pct
                    stats['r1_mapped_count'] = mapped_cnt
                    stats['r1_mapped_percent'] = pct
                    stats['r1_unambiguous_count'] = mapped_cnt
                    stats['r1_unambiguous_percent'] = pct
                    stats['r2_mapped_count'] = mapped_cnt
                    stats['r2_mapped_percent'] = pct
                    stats['r2_unambiguous_count'] = mapped_cnt
                    stats['r2_unambiguous_percent'] = pct
                elif "properly paired (" in line:
                    pp_cnt = int(line.split()[0])
                    try:
                        pct = float(line.split('(')[1].split('%')[0].strip())
                    except (IndexError, ValueError):
                        pct = 0.0
                    stats['properly_paired_count'] = pp_cnt
                    stats['properly_paired_percent'] = pct
                    stats['pairing_mated_pairs_count'] = pp_cnt // 2
                    stats['pairing_mated_pairs_percent'] = pct
                elif "singletons (" in line:
                    sing_cnt = int(line.split()[0])
                    try:
                        pct = float(line.split('(')[1].split('%')[0].strip())
                    except (IndexError, ValueError):
                        pct = 0.0
                    stats['singletons_count'] = sing_cnt
                    stats['singletons_percent'] = pct
            return stats

        for line in all_lines:

                # Global stat
                if line.startswith("Reads Used:"):
                    parts = line.split()
                    reads_used_val = int(parts[2])
                    stats['input_reads_processed'] = reads_used_val
                    stats['input_read_pairs'] = reads_used_val // 2
                
                # Identify section
                elif line.startswith("Pairing data:"):
                    current_section = "pairing"
                    continue
                elif line.startswith("Read 1 data:"):
                    current_section = "read1"
                    continue
                elif line.startswith("Read 2 data:"):
                    current_section = "read2"
                    continue
                elif line.startswith("Coverage Gaps:") or line.startswith("----") or line.startswith("Match Rate:") : # Section terminators or irrelevant sections for these primary stats
                    current_section = None 
                
                parts = line.split()
                if len(parts) < 2: # Stats lines have at least key: value
                    continue
                
                # Normalize the key by removing colon and lowercasing first word(s)
                key_part = parts[0].lower()
                if len(parts) > 1 and parts[1].lower() == "pairs:":  # Handles "mated pairs:", "bad pairs:"
                     key_part = f"{parts[0].lower()} {parts[1].lower().replace(':', '')}"
                else:
                     key_part = parts[0].lower().replace(':', '')

                if current_section == "pairing":
                    try:
                        if key_part == "mated pairs":
                            stats['pairing_mated_pairs_percent'] = float(parts[2].strip('%'))
                            stats['pairing_mated_pairs_count'] = int(parts[3])
                        elif key_part == "bad pairs":
                            stats['pairing_bad_pairs_percent'] = float(parts[2].strip('%'))
                            stats['pairing_bad_pairs_count'] = int(parts[3])
                        elif key_part == "insert" and len(parts) > 2 and parts[1].lower() == "size" and parts[2].lower() == "avg":  # "insert size avg:"
                            stats['pairing_insert_size_avg'] = float(parts[-1])  # Last element is the value
                    except (IndexError, ValueError) as e:
                        logging.warning(f"Could not parse pairing data line in {filepath}: '{line}' (key: '{key_part}') - Error: {e}")
                
                elif current_section == "read1":
                    if len(parts) < 3: continue  # Expecting: "key: pct reads    num reads   pct bases   num bases"
                    try:
                        if key_part == "mapped":
                            stats['r1_mapped_percent'] = float(parts[1].strip('%'))
                            stats['r1_mapped_count'] = int(parts[2])
                        elif key_part == "unambiguous":
                            stats['r1_unambiguous_percent'] = float(parts[1].strip('%'))
                            stats['r1_unambiguous_count'] = int(parts[2])
                        elif key_part == "ambiguous":
                            stats['r1_ambiguous_percent'] = float(parts[1].strip('%'))
                            stats['r1_ambiguous_count'] = int(parts[2])
                    except (IndexError, ValueError) as e:
                        logging.warning(f"Could not parse R1 data for key '{key_part}' in {filepath}: '{line}' - Error: {e}")

                elif current_section == "read2":
                    if len(parts) < 3: continue
                    try:
                        if key_part == "mapped":
                            stats['r2_mapped_percent'] = float(parts[1].strip('%'))
                            stats['r2_mapped_count'] = int(parts[2])
                        elif key_part == "unambiguous":
                            stats['r2_unambiguous_percent'] = float(parts[1].strip('%'))
                            stats['r2_unambiguous_count'] = int(parts[2])
                        elif key_part == "ambiguous":
                            stats['r2_ambiguous_percent'] = float(parts[1].strip('%'))
                            stats['r2_ambiguous_count'] = int(parts[2])
                    except (IndexError, ValueError) as e:
                        logging.warning(f"Could not parse R2 data for key '{key_part}' in {filepath}: '{line}' - Error: {e}")

        # Calculate unmapped stats for R1
        if isinstance(stats['input_read_pairs'], int) and isinstance(stats['r1_mapped_count'], int):
            stats['r1_unmapped_count'] = stats['input_read_pairs'] - stats['r1_mapped_count']
            if stats['input_read_pairs'] > 0:
                stats['r1_unmapped_percent'] = round((stats['r1_unmapped_count'] / stats['input_read_pairs']) * 100, 4)
            elif stats['r1_unmapped_count'] == 0: # If 0 input pairs and 0 mapped, then 0% unmapped
                stats['r1_unmapped_percent'] = 0.0

        # Calculate unmapped stats for R2
        if isinstance(stats['input_read_pairs'], int) and isinstance(stats['r2_mapped_count'], int):
            stats['r2_unmapped_count'] = stats['input_read_pairs'] - stats['r2_mapped_count']
            if stats['input_read_pairs'] > 0:
                stats['r2_unmapped_percent'] = round((stats['r2_unmapped_count'] / stats['input_read_pairs']) * 100, 4)
            elif stats['r2_unmapped_count'] == 0:
                stats['r2_unmapped_percent'] = 0.0

    except FileNotFoundError:
        logging.error(f"BBMap stats file not found: {filepath}")
    except Exception as e:
        logging.error(f"Error parsing BBMap stats file {filepath}: {e}", exc_info=True)
    
    logging.debug(f"Parsed BBMap stats for {filepath}: {stats}")
    return stats

def parse_sambamba_markdup_log(filepath):
    """
    Parses the number of duplicate reads removed from Sambamba markdup log file.

    Example sambamba markdup log file:
    sambamba 1.0.1
    by Artem Tarasov and Pjotr Prins (C) 2012-2023
        LDC 1.39.0 / DMD v2.109.1 / LLVM17.0.6 / bootstrap LDC - the LLVM D compiler (1.39.0)

    finding positions of the duplicate reads in the file...
    sorted 2480006 end pairs
        and 0 single ends (among them 0 unmatched pairs)
    collecting indices of duplicate reads...   done in 290 ms
    found 1772924 duplicates
    collected list of positions in 0 min 2 sec
    removing duplicates...
    collected list of positions in 0 min 6 sec
    """
    logging.debug(f"Parsing Sambamba markdup log: {filepath}")
    stats = {'duplicates_removed': 'N/A'}

    if not filepath:
        logging.warning("Sambamba markdup log filepath was None. Returning N/A stats.")
        return stats

    try:
        log_content = Path(filepath).read_text(encoding='utf-8')
        
        # Use regex to find the line and extract the number of duplicates
        match = re.search(r"found (\d+) duplicates", log_content, re.IGNORECASE) 
        
        if match:
            try:
                stats['duplicates_removed'] = int(match.group(1))
                logging.debug(f"Found duplicates removed: {stats['duplicates_removed']}")
            except ValueError:
                logging.warning(f"Could not convert duplicate count to integer in {filepath}: Found '{match.group(1)}'")
                stats['duplicates_removed'] = 'N/A'
        else:
            logging.warning(f"Could not find the 'found X duplicates' line in Sambamba log: {filepath}")

    except FileNotFoundError:
        logging.error(f"Sambamba markdup log file not found: {filepath}")
    except Exception as e:
        logging.error(f"Error parsing Sambamba markdup log file {filepath}: {e}", exc_info=True)
    
    return stats

def get_genome_stats_from_bam(filepath):
    """
    Gets total genome size (sum of lengths) and contig count 
    from the @SQ lines in a BAM file header using pysam.
    """
    logging.debug(f"Getting genome stats from BAM: {filepath}")
    stats = {'genome_size': 'N/A', 'contig_count': 'N/A'}

    if not filepath:
        logging.warning("BAM filepath was None for genome stats extraction.")
        return stats

    try:
        genome_size = 0
        contig_count = 0
        with pysam.AlignmentFile(filepath, "rb") as bamfile:
            if not bamfile.header or 'SQ' not in bamfile.header:
                logging.warning(f"No @SQ header lines found in BAM file: {filepath}")
                return stats

            header_sq_list = bamfile.header['SQ']
            contig_count = len(header_sq_list)
            stats['contig_count'] = contig_count

            # Sum lengths from LN tag in each SQ dictionary
            for sq_record in header_sq_list:
                if 'LN' in sq_record:
                    try:
                        genome_size += int(sq_record['LN'])
                    except (ValueError, TypeError):
                        logging.warning(f"Could not parse LN value as integer in @SQ line for {sq_record.get('SN', 'Unknown Contig')} in {filepath}: Value='{sq_record['LN']}'")
                else:
                    logging.warning(f"Missing 'LN' (length) tag in @SQ line for {sq_record.get('SN', 'Unknown Contig')} in {filepath}")
            
            stats['genome_size'] = genome_size

    except FileNotFoundError:
        logging.error(f"BAM file not found for header parsing: {filepath}")
    except Exception as e:
        logging.error(f"Unexpected error parsing BAM header for {filepath}: {e}", exc_info=True)
    
    logging.debug(f"Parsed Genome Stats for {filepath}: {stats}")
    return stats

def read_coverage_value(filepath):
    """
    Reads the mean coverage float value from the 4th column of a 
    single-line, tab-separated file.
    Expected format: tax_id\tname\tgenome_basename\tcoverage
    """
    logging.debug(f"Reading coverage value from: {filepath}")
    
    if not filepath:
        logging.warning("Coverage filepath was None. Returning 'N/A'.")
        return 'N/A'

    try:
        filepath_obj = Path(filepath)
        content = filepath_obj.read_text(encoding='utf-8').strip()

        if not content:
            logging.warning(f"Coverage file is empty: {filepath}. Returning 'N/A'.")
            return 'N/A'
             
        parts = content.split('\t')
        if len(parts) == 4:
            try:
                coverage_value = float(parts[3])
                logging.debug(f"Parsed coverage value {coverage_value} from {filepath}")
                return coverage_value
            except ValueError:
                logging.error(f"Could not convert coverage value '{parts[3]}' to float in {filepath}. Returning 'N/A'.")
                return 'N/A'
        else:
            logging.error(f"Coverage file {filepath} did not contain expected 4 tab-separated columns. Found {len(parts)}. Content: '{content}'. Returning 'N/A'.")
            return 'N/A'

    except FileNotFoundError:
        logging.error(f"Coverage file not found: {filepath}. Returning 'N/A'.")
        return 'N/A'
    except Exception as e:
        logging.error(f"Unexpected error reading coverage file {filepath}: {e}", exc_info=True)
        return 'N/A'

def parse_sim_params_json(filepath):
    """
    Parses the simulation parameters JSON file created by 
    estimate_simulation_params rule. Returns a dictionary.
    """
    logging.debug(f"Parsing simulation parameters JSON: {filepath}")
    stats = {
        'sim_params_estimated_insert_size': 'N/A',
        'sim_params_estimated_std_insert_size': 'N/A',
        'sim_params_estimated_read_length': 'N/A',
        'sim_params_used_insert_size': 'N/A',
        'sim_params_used_std_insert_size': 'N/A',
        'sim_params_used_read_length': 'N/A'}

    if not filepath:
        logging.warning("Simulation parameters JSON filepath was None.")
        return stats

    try:
        filepath_obj = Path(filepath)
        with filepath_obj.open('r', encoding='utf-8') as f:
            data = json.load(f)
            if isinstance(data, dict):
                logging.debug(f"Successfully parsed JSON from {filepath}")
                stats['sim_params_estimated_insert_size'] = data.get('estimated_mean_insert_size', 'N/A')
                stats['sim_params_estimated_std_insert_size'] = data.get('estimated_std_insert_size', 'N/A')
                stats['sim_params_estimated_read_length'] = data.get('estimated_mean_read_length', 'N/A')
                stats['sim_params_used_insert_size'] = data.get('used_mean_insert_size', 'N/A')
                stats['sim_params_used_std_insert_size'] = data.get('used_std_insert_size', 'N/A')
                stats['sim_params_used_read_length'] = data.get('used_mean_read_length', 'N/A')
            else:
                logging.warning(f"JSON file did not contain a dictionary (root object): {filepath}")
                return stats

    except FileNotFoundError:
        logging.error(f"Simulation parameters JSON file not found: {filepath}")
    except Exception as e:
        logging.error(f"Unexpected error reading/parsing simulation parameters JSON {filepath}: {e}", exc_info=True)
    
    return stats

def parse_pirs_log(filepath):
    """Parses the number of simulated read pairs from a pIRS log file."""
    logging.debug(f"Parsing pIRS log: {filepath}")
    stats = {'simulated_pairs': 'N/A'}

    if not filepath:
        logging.warning("pIRS log filepath was None. Returning N/A stats.")
        return stats

    try:
        filepath_obj = Path(filepath)
        # Regex to find the line and capture the number of pairs
        line_regex = re.compile(r"^\s*\[pIRS\]\s*Read pairs simulated:\s+(\d+)", re.IGNORECASE)
        
        found_metric = False
        with filepath_obj.open('r', encoding='utf-8') as f:
            for line in f:
                match = line_regex.search(line.strip())
                if match:
                    try:
                        stats['simulated_pairs'] = int(match.group(1))
                        found_metric = True
                        logging.debug(f"Found simulated pairs: {stats['simulated_pairs']} in {filepath}")
                        break
                    except ValueError:
                        logging.warning(f"Could not convert simulated pairs value '{match.group(1)}' to integer in {filepath}")
                        stats['simulated_pairs'] = 'N/A'
                        found_metric = True
                        break
                    except IndexError:
                        logging.warning(f"Regex matched but could not find group 1 in line '{line.strip()}' in {filepath}")
                        stats['simulated_pairs'] = 'N/A'
                        found_metric = True
                        break
        
        if not found_metric:
             logging.warning(f"Did not find 'Read pairs simulated:' line in pIRS log: {filepath}")

    except FileNotFoundError:
        logging.error(f"pIRS log file not found: {filepath}")
    except Exception as e:
        logging.error(f"Error reading/parsing pIRS log file {filepath}: {e}", exc_info=True)
        stats['simulated_pairs'] = 'N/A'
    
    return stats

def parse_extract_ids_log(filepath):
    """
    Parses the selmeout.py log (extract_correct_read_ids log) to get the number of 
    correctly classified paired-end reads and total paired-end reads processed 
    by the script. Also calculates percentage correctly classified.
    """
    logging.debug(f"Parsing extract_correct_read_ids log (selmeout.py log): {filepath}")
    stats = {
        'correctly_classified_pe_reads': 'N/A',
        'correctly_classified_pe_perc': 'N/A',
        'total_pe_reads_processed': 'N/A'}

    if not filepath:
        logging.warning("Log filepath was None. Returning N/A stats.")
        return stats

    try:
        filepath_obj = Path(filepath)
        # Regex to capture the two numbers from the target line.
        # Example: "Found 3071934 reads matching the specifications (total reads: 3124521)."
        line_regex = re.compile(r"Found (\d+) reads matching the specifications \(total reads: (\d+)\)\.")
        
        found_metrics = False
        with filepath_obj.open('r', encoding='utf-8') as f:
            for line in f:
                match = line_regex.search(line)
                if match:
                    try:
                        corr_class_pe_reads = int(match.group(1))
                        total_pe_reads_processed = int(match.group(2))
                        found_metrics = True
                        
                        # Percentage calculation
                        if total_pe_reads_processed > 0:
                            corr_class_perc = round((corr_class_pe_reads / total_pe_reads_processed) * 100, 2)
                        else:
                            corr_class_perc = 'N/A'
                        
                        # Store the values in the dict
                        stats['correctly_classified_pe_reads'] = corr_class_pe_reads
                        stats['correctly_classified_pe_perc'] = corr_class_perc
                        stats['total_pe_reads_processed'] = total_pe_reads_processed
                        stats['correctly_classified_reads'] = corr_class_pe_reads
                        stats['correctly_classified_perc'] = corr_class_perc
                        stats['total_sliced_reads'] = total_pe_reads_processed
                        
                        logging.debug(
                            f"Found these values: "
                            f"Correct PE: {corr_class_pe_reads}, "
                            f"Total PE processed: {total_pe_reads_processed}, "
                            f"Percentage correct PE: {corr_class_perc}")
                        
                        break

                    except (ValueError, IndexError) as e:
                        logging.warning(
                            f"Error parsing numbers from matched line in {filepath}: "
                            f"'{line.strip()}'. Error: {e}")
                        stats['correctly_classified_pe_reads'] = 'N/A'
                        stats['total_pe_reads_processed'] = 'N/A'
                        found_metrics = True
                        break 
        
        if not found_metrics:
            logging.warning(f"Did not find the target metrics line in the log file: {filepath}")
        
    except FileNotFoundError:
        logging.error(f"Log file not found: {filepath}")
    except Exception as e:
        logging.error(f"Error reading/parsing the log file {filepath}: {e}", exc_info=True)
        stats['correctly_classified_pe_reads'] = 'N/A'
        stats['total_pe_reads_processed'] = 'N/A'
        
    return stats

def parse_bed_file(filepath, genome_size):
    """
    Calculates the total length covered by regions in a BED file.
    Assumes a tab-separated file with at least 3 columns (chrom, start, end), 
    and uses the 2nd (start, index 1) and 3rd (end, index 2) columns.
    Also calculates the percentage of genome that is classifiable.
    """
    logging.debug(f"Parsing BED file for total length: {filepath}")
    stats = {
        'classif_regions_total_length': 'N/A',
        'classif_regions_percent_genome': 'N/A',
        'non_classif_regions_length': 'N/A',
        'non_classif_regions_percent_genome': 'N/A'}

    if not filepath:
        logging.warning("BED filepath was None. Returning N/A stats.")
        return stats

    try:
        filepath_obj = Path(filepath)
        df = pd.read_csv(
            filepath_obj, 
            sep='\t', 
            header=None, 
            usecols=[1, 2],
            comment='#') 
        
        # Col headers
        df.columns = ['start', 'end']

        # Check if any data rows were actually read
        if df.empty:
             logging.warning(f"BED file contained no valid data rows: {filepath}")
             return stats

        # Convert columns to numeric, coercing errors to NaN
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')

        # Remove rows where start or end could not be converted to numeric
        initial_rows = len(df)
        df.dropna(subset=['start', 'end'], inplace=True)
        rows_dropped = initial_rows - len(df)
        if rows_dropped > 0:
             logging.warning(f"Removed {rows_dropped} rows with non-numeric start/end values from BED file: {filepath}")
        
        # Check if any valid numeric rows remain
        if df.empty:
             logging.warning(f"BED file contained no rows with valid numeric start/end values: {filepath}")
             return stats

        # Calculate total length: sum of (end - start) for all valid rows
        total_length = (df['end'] - df['start']).sum()
        
        # Ensure result is a standard Python number
        total_length = int(total_length)
        stats['classif_regions_total_length'] = total_length

        # Percentages for classifiable regions
        if isinstance(genome_size, (int)) and genome_size > 0:
            stats['classif_regions_percent_genome'] = round((total_length / genome_size) * 100, 2)
            stats['non_classif_regions_length'] = genome_size - total_length
            stats['non_classif_regions_percent_genome'] = round(((genome_size - total_length) / genome_size) * 100, 2)
        else:
            logging.warning(f"Invalid genome_size ('{genome_size}') provided for percentage calculation for BED file: {filepath}")

        logging.debug(f"Calculated total length {stats['classif_regions_total_length']} bp from BED: {filepath}")

    except FileNotFoundError:
        logging.error(f"BED file not found: {filepath}")
    except Exception as e:
        logging.error(f"Unexpected error processing BED file {filepath}: {e}", exc_info=True)
        
    return stats

def generate_plotly_fig(summary_tsv_path, sample_info):
    """
    Generates a Plotly scatter plot figure from the per-contig summary TSV file,
    including correlations and scale toggles.
    Returns a plotly.graph_objects.Figure object or None on error.
    'sample_info' is the mapping spec df (containing keys like 'name', 'tax_id', 'genome_basename').
    """
    logging.debug(f"Attempting to generate plot from: {summary_tsv_path}")
    fig_on_error = None

    if not summary_tsv_path:
        logging.warning(f"Summary TSV path was None for sample {sample_info.get('name', 'N/A')}.")
        return fig_on_error

    try:
        filepath_obj = Path(summary_tsv_path)
        try:
            df = pd.read_csv(filepath_obj, sep='\t')
        except pd.errors.EmptyDataError:
             logging.warning(f"Summary data file for plotting is empty: {summary_tsv_path}")
             return fig_on_error

        # Check for required columns
        X_COL = 'LengthClassifiableRegion'
        Y_COL = 'NumReadsOverlapClassifiable'
        HOVER_COL = 'Contig'
        required_columns = {X_COL, Y_COL, HOVER_COL}
        if df.empty or not required_columns.issubset(df.columns):
            logging.error(f"Input TSV {summary_tsv_path} is empty or missing required columns {required_columns} for plotting.")
            return fig_on_error

        logging.info(f"Generating plot for {len(df)} contigs from {summary_tsv_path}")

        # --- Data cleaning for correlation ---
        # Convert to numeric first
        df[X_COL] = pd.to_numeric(df[X_COL], errors='coerce')
        df[Y_COL] = pd.to_numeric(df[Y_COL], errors='coerce')
        df_filt = df.dropna(subset=[X_COL, Y_COL]) # Remove rows where conversion failed
        df_filt = df_filt[df_filt[X_COL] > 0]      # Filter length > 0 for correlation

        # --- Calculate Correlations ---
        pearson_text = "Pearson r: "
        spearman_text = "Spearman ρ: "
        if len(df_filt) > 1: # Need at least 2 points
            
            # Pearson
            try:
                pearson_r, p_val_p = stats.pearsonr(df_filt[X_COL], df_filt[Y_COL])
                pearson_text += f'{pearson_r:.3f} (p={p_val_p:.2g})' if pd.notna(pearson_r) else "N/A (NaN)"
            except ValueError as e:
                logging.warning(f"Could not calculate Pearson correlation for {summary_tsv_path}: {e}")
                pearson_text += "N/A (Error)"
            
            # Spearman
            try:
                spearman_r, p_val_s = stats.spearmanr(df_filt[X_COL], df_filt[Y_COL])
                spearman_text += f'{spearman_r:.3f} (p={p_val_s:.2g})' if pd.notna(spearman_r) else "N/A (NaN)"
            except ValueError as e:
                logging.warning(f"Could not calculate Spearman correlation for {summary_tsv_path}: {e}")
                spearman_text += "N/A (Error)"
            
        else: # Handle too few points
            pearson_text += "N/A (too few points)"
            spearman_text += "N/A (too few points)"
            logging.info(f"Correlation not calculated for {summary_tsv_path} (too few points with length > 0).")

        # --- Create Plot Title ---
        plot_title = (
            f"Reads vs Classifiable Length: {sample_info.get('name', 'N/A')}<br>"
            f"TaxID: {sample_info.get('tax_id', 'N/A')} | Genome: {sample_info.get('genome_basename', 'N/A')}<br>"
            f"{pearson_text} | {spearman_text}")

        # --- Create Plot using Plotly Express ---
        fig = px.scatter(
                df,
                x=X_COL,
                y=Y_COL,
                title=plot_title,
                labels={
                    X_COL: "Length of Classifiable Region (bp)",
                    Y_COL: "Number of Overlapping Real Reads (Paired-End)"},
                hover_data=[HOVER_COL])

        # Customize layout and traces
        fig.update_layout(
            title_font_size=14,
            title_x=0.5, # Center title
            hovermode='closest',
            margin=dict(l=50, r=20, t=60, b=50))
        fig.update_traces(marker=dict(size=8, opacity=0.7), selector=dict(mode='markers'))

        # Add Buttons for Axis Scale Toggle
        fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons", direction="left",
                    buttons=list([
                        dict(args=[{"xaxis.type": "linear", "yaxis.type": "linear"}], label="Linear Scale", method="relayout"),
                        dict(args=[{"xaxis.type": "log", "yaxis.type": "linear"}], label="Log X", method="relayout"),
                        dict(args=[{"xaxis.type": "linear", "yaxis.type": "log"}], label="Log Y", method="relayout"),
                        dict(args=[{"xaxis.type": "log", "yaxis.type": "log"}], label="Log X & Y", method="relayout")
                    ]),
                    pad={"r": 10, "t": 10}, showactive=True,
                    x=0.01, xanchor="left", y=1.15, yanchor="top"
                ),
            ]
        )

        logging.debug(f"Successfully generated Plotly figure from {summary_tsv_path}")
        return fig

    # --- Error Handling ---
    except FileNotFoundError:
         logging.error(f"Summary TSV for plotting not found: {summary_tsv_path}")
         return fig_on_error
    except KeyError as e:
         logging.error(f"Missing expected column during plot generation from {summary_tsv_path}: {e}", exc_info=True)
         return fig_on_error
    except Exception as e:
         logging.error(f"Unexpected error generating plot from {summary_tsv_path}: {e}", exc_info=True)
         return fig_on_error

# --- Helper function to find a file in the rule input list that corresponds to a specific combination of wildcards ---
def find_input_file(pattern_key, wildcards_dict, input_file_list, path_patterns_dict):
    """
    Reconstructs an expected file path using a pattern and wildcards, then finds it in the provided list of 
    input files. Returns the file path string or None if not found or error.
    """
    try:
        # Get the filepath pattern string using the key
        pattern_template = path_patterns_dict[pattern_key]
        
        # Format the pattern string with the current combination of wildcards
        # wildcards_dict looks like {'sanitized_name': 'X', 'tax_id': 'Y', 'genome_basename': 'Z'}
        expected_path_str = pattern_template.format(**wildcards_dict)
        expected_path_obj = Path(expected_path_str)

        # Find the expected path in the list provided by snakemake
        for file_path_str in input_file_list:
            if Path(file_path_str) == expected_path_obj:
                return file_path_str
        
        logging.warning(
            f"Could not find the expected file for '{pattern_key}' matching '{expected_path_str}' in the rule input list for the wildcard combination:\n\
            sanitized_name: '{wildcards_dict.get('sanitized_name', 'N/A')}'\n\
            tax_id: '{wildcards_dict.get('tax_id', 'N/A')}'\n\
            genome_basename: '{wildcards_dict.get('genome_basename', 'N/A')}'.")
        return None
    
    except KeyError: # If pattern_key is missing from path_patterns_dict
        logging.error(f"Error finding input file: Path pattern key '{pattern_key}' not found in params.path_patterns for wildcards {wildcards_dict}.")
        return None
    except Exception as e:
        logging.error(f"Error constructing/finding path for '{pattern_key}' with wildcards {wildcards_dict}: {e}", exc_info=True)
        return None

# --- Main Processing Loop ---
logging.info("Starting per-sample processing loop.")
samples_data = [] # List to hold dictionaries, one per sample (wildcard combination)

if mapping_df.empty:
    logging.error("Mapping specification DataFrame is empty. No samples to process.")
    sys.exit(1)
else:
    for index, map_row in mapping_df.iterrows():
        # map_row contains: tax_id, name, fasta_file, genome_basename, sanitized_name
        sample_wildcards = map_row.to_dict()  # This dict can be used for .format(**sample_wildcards)

        logging.info(f"Processing sample {sample_wildcards['name']}: TaxID={sample_wildcards['tax_id']}, Genome={sample_wildcards['genome_basename']})")

        # Store parsed data for this sample
        current_sample_metrics = {'info': sample_wildcards.copy()}

        # --- MAPPING REAL READS ---
        real_stats_file = find_input_file("real_stats", sample_wildcards, real_stats_files, path_patterns)
        current_sample_metrics['real_mapping_stats'] = parse_bbmap_stats(real_stats_file)

        markdup_log_file = find_input_file("markdup_logs", sample_wildcards, markdup_logs, path_patterns)
        current_sample_metrics['dedup_stats'] = parse_sambamba_markdup_log(markdup_log_file)

        dedup_bam_file = find_input_file("dedup_bams", sample_wildcards, dedup_bams, path_patterns)
        current_sample_metrics['genome_info'] = get_genome_stats_from_bam(dedup_bam_file)

        coverage_file = find_input_file("coverage_files", sample_wildcards, coverage_files, path_patterns)
        current_sample_metrics['real_coverage_depth'] = read_coverage_value(coverage_file)

        # --- SLICING / TILING ---
        current_sample_metrics['tiling_info'] = {
            'window_size': report_config_dict.get('slice_window', 300),
            'step_size': report_config_dict.get('slice_step', 50)
        }

        # --- SIMULATION (Legacy Fallback) ---
        sim_params_file = find_input_file("sim_params_jsons", sample_wildcards, sim_params_jsons, path_patterns)
        current_sample_metrics['sim_params_data'] = parse_sim_params_json(sim_params_file)
        current_sample_metrics['sim_target_coverage'] = report_config_dict.get('sim_cov_target', 'N/A')

        pirs_log_file = find_input_file("pirs_logs", sample_wildcards, pirs_logs, path_patterns)
        current_sample_metrics['pirs_log_stats'] = parse_pirs_log(pirs_log_file)

        # --- CLASSIFICATION (Simulated Reads) ---
        correct_ids_log_file = find_input_file("correct_ids_logs", sample_wildcards, correct_ids_logs, path_patterns)
        current_sample_metrics['sim_classified_reads'] = parse_extract_ids_log(correct_ids_log_file)

        # --- MAPPING SIMULATED READS (Correctly Classified Subset) ---
        sim_stats_file = find_input_file("sim_stats", sample_wildcards, sim_stats_files, path_patterns)
        current_sample_metrics['simulated_mapping_stats'] = parse_bbmap_stats(sim_stats_file)

        # --- CLASSIFIABLE REGIONS ---
        classifiable_bed_file = find_input_file("classifiable_beds", sample_wildcards, classifiable_beds, path_patterns)
        current_sample_metrics['classifiable_regions'] = parse_bed_file(classifiable_bed_file, genome_size=current_sample_metrics['genome_info'].get('genome_size', 0))

        # --- PLOT ---
        summary_tsv_file = find_input_file("summary_tsvs", sample_wildcards, summary_tsvs, path_patterns)
        plot_fig = generate_plotly_fig(summary_tsv_file, sample_wildcards)
        if plot_fig:
            current_sample_metrics['plot_html_div'] = pio.to_html(
                plot_fig,
                full_html=False,
                include_plotlyjs='cdn')
        else:
            current_sample_metrics['plot_html_div'] = "<p>Plot could not be generated.</p>"
        
        samples_data.append(current_sample_metrics)

# --- HTML Generation ---
logging.info("Generating final HTML report.")
html_content = """<html><head><title>Report Generation Error</title></head>
                <body><h1>Error: Report Generation Failed</h1>
                <p>Could not render Jinja2 template. Check logs.</p></body></html>"""

try:
    # Set up Jinja2 environment
    script_dir = Path(__file__).parent
    template_dir = script_dir / "templates"
    if not template_dir.is_dir():
         logging.error(f"Jinja2 template directory not found at expected location: {template_dir}")
         raise FileNotFoundError(f"Template directory not found: {template_dir}")

    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=True)

    # Add custom filter for formatting numbers nicely or showing 'N/A'
    def format_metric(value, precision=2, na_value='N/A'):
        """Jinja filter to format numbers with commas and precision, or return N/A."""
        if value is None or value == 'N/A':
            return f'<span class="value-na">{na_value}</span>'
        try:
            # Attempt to convert to float first for flexible input
            num_value = float(value)
            # Format based on precision
            if precision == 0:
                return f"{int(num_value):,}"          # Comma separators, no decimal
            else:
                return f"{num_value:,.{precision}f}"  # Comma separators, specified decimal places
        except (ValueError, TypeError):
            # If conversion fails, return original value (might be error string) or N/A
             logging.warning(f"format_metric filter could not format value: {value}")
             return f'<span class="value-na">{na_value}</span>'

    # Register the filter
    env.filters['format_metric'] = format_metric

    # Load the template file
    template_file_name = "report_template.html"
    template = env.get_template(template_file_name)
    logging.info(f"Loaded Jinja2 template: {template_file_name}")

    # Prepare the full data context for the template
    template_data = {
         'report_generation_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S %Z'),
         'general_info': general_info, # Dictionary of general run info
         'samples': samples_data       # List of sample dictionaries
    }

    # Render HTML
    html_content = template.render(template_data)
    logging.info(f"Successfully rendered HTML report template.")

except Exception as e:
     logging.error(f"Error generating HTML report with Jinja2: {e}", exc_info=True)


# --- Write Output ---
logging.info(f"Writing HTML report to {report_html_file}")
try:
    with open(report_html_file, "w", encoding='utf-8') as f_out:
        f_out.write(html_content)
    logging.info("Report generation script finished successfully.")
except Exception as e:
     logging.error(f"Failed to write final report {report_html_file}: {e}", exc_info=True)
     sys.exit(1)