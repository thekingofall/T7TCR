#!/usr/bin/env python3

import os
import gzip
import glob
import re
import datetime
import shutil
import itertools
import sys
from collections import defaultdict

# --- Configuration ---
DEFAULT_READS = 1000000 # Set to None to process all reads
UMI_LEN = 8 # Set UMI length based on the NNNNNNNN segments

# Define patterns and their fragments (used for identifying the overall structure type)
# These help broadly categorize the read pair structure first.
PATTERNS_DATA = {
    1: ("TAATACGACTCACTATAGGG", "CGAAACATCGGCCACCCC", "CAACCCTGCGAC"),
    2: ("CGCAATGAAGTCGCAGGGTTG", "GGGGTGGCCGATGTTTCG", "CCCTATAGTGAGTCGTATTA"),
    3: ("TAATACGACTCACTATAGGGGTGGCCGATGTTTCG", "CCC", "CAACCCTGCGACTTCA"),
    4: ("CAATGAAGTCGCAGGGTTG", "GGG", "CGAAACATCGGCCACCCCTATAGTGAGTCGTATTA")
}

# --- FINAL CORRECTED UMI Flanking Sequences ---
# Define the flanking sequences for the UMI segment based *directly* on the user's
# last set of unambiguous examples.
# Format: {pattern_index: (sequence_before_UMI, sequence_after_UMI)}
UMI_FLANKING_SEQUENCES = {
    # Pattern 1 UMI: T7 promoter NNNNNNNN TSO_start
    1: ("TAATACGACTCACTATAGGG", "CGAAACATCGGCCAC"),
    # Pattern 2 UMI: Barcode_end NNNNNNNN P5_adapter_start
    2: ("GTGGCCGATGTTTCG", "CCCTATAGTGAGTCGTATTA"),
    # Pattern 3 UMI: Barcode_end NNNNNNNN CCC_linker
    3: ("GTGGCCGATGTTTCG", "CCC"),
    # Pattern 4 UMI: GGG_linker NNNNNNNN TSO_start
    4: ("GGG", "CGAAACATCGGCCAC"),
}

# --- Setup Directories and Report File ---
HEAD_LABEL = f"head{DEFAULT_READS}" if DEFAULT_READS else "all"
MATCHED_BASE_DIR = f"{HEAD_LABEL}_matched_reads_umi_format_py"
REPORT_FILE = f"{HEAD_LABEL}_analysis_report_merged_umi_py.txt"


# --- Helper Functions ---

def read_fastq_records(file_handle):
    """Generator to yield 4-line FASTQ records from a file handle."""
    while True:
        lines = list(itertools.islice(file_handle, 4))
        if not lines or len(lines) < 4:
            break
        yield tuple(line.strip() for line in lines)

def extract_umi(seq, frag_before, frag_after, umi_len):
    """
    Extracts a UMI of specified length between two flanking fragments using regex.
    Returns UMI sequence or 'N'*umi_len if not found or if flanks are invalid.
    """
    if not frag_before or not frag_after:
        return 'N' * umi_len
    try:
        f_before_esc = re.escape(frag_before)
        f_after_esc = re.escape(frag_after)
        regex = re.compile(f"{f_before_esc}(.{{{umi_len}}}){f_after_esc}")
        match = regex.search(seq)
        if match:
            return match.group(1)
    except re.error as e:
        print(f"  Regex error extracting UMI between '{frag_before}' and '{frag_after}': {e}", file=sys.stderr)
    except TypeError as e:
        print(f"  Type error during UMI extraction (likely None flanks) between '{frag_before}' and '{frag_after}': {e}", file=sys.stderr)

    return 'N' * umi_len

def get_umi_from_pair(r1_seq, r2_seq, frag_before, frag_after, priority='R1', umi_len=8):
    """
    Extracts UMI based on priority/fallback logic between R1 and R2 reads.
    priority: 'R1' (try R1 first), 'R2' (try R2 first)
    """
    primary_seq = r1_seq if priority == 'R1' else r2_seq
    fallback_seq = r2_seq if priority == 'R1' else r1_seq
    n_fill = 'N' * umi_len

    umi = extract_umi(primary_seq, frag_before, frag_after, umi_len)

    if umi == n_fill:
        umi = extract_umi(fallback_seq, frag_before, frag_after, umi_len)

    return umi


# --- Main Processing Logic ---

def main():
    os.makedirs(MATCHED_BASE_DIR, exist_ok=True)

    total_lines_analyzed_per_file = DEFAULT_READS * 4 if DEFAULT_READS else "All"
    reads_analyzed_desc = DEFAULT_READS if DEFAULT_READS else "All"

    try:
        with open(REPORT_FILE, 'w') as rp:
            rp.write("===== Hi-C Samples Analysis Report (Python Version) =====\n")
            rp.write(f"Analysis Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            rp.write(f"Reads analyzed per sample (pairs): {reads_analyzed_desc}\n")
            if isinstance(total_lines_analyzed_per_file, int):
                 rp.write(f"Lines analyzed per R1/R2 file: {total_lines_analyzed_per_file}\n")
            rp.write(f"UMI Length: {UMI_LEN}\n")
            rp.write("UMI Flanking Sequences (Before, After) - FINAL based on user examples:\n") # Indicate final correction
            for idx, flanks in UMI_FLANKING_SEQUENCES.items():
                 rp.write(f"  Pattern {idx}: ('{flanks[0]}', '{flanks[1]}')\n")
            rp.write("\n=== Sample Statistics ===\n\n")
            rp.write("Sample\tTotal_Pairs_Analyzed\tMatched_Pairs\tPercentage (% of pairs)\tPattern_Index\n")
            rp.write("-" * 90 + "\n")

        report_data = []
        r1_files = sorted(glob.glob('*_1.fq.gz'))
        if not r1_files:
             print("No *_1.fq.gz files found in the current directory.")
             return

        for r1_gz in r1_files:
            # --- Base Name Extraction ---
            base_name = None
            if r1_gz.endswith('_1.fq.gz'):
                base_name = r1_gz[:-len('_1.fq.gz')]
            elif r1_gz.endswith('_R1.fq.gz'): # Handle _R1 suffix
                 base_name = r1_gz[:-len('_R1.fq.gz')]
            elif r1_gz.endswith('_R1_001.fq.gz'): # Handle common Illumina suffix
                 base_name = r1_gz[:-len('_R1_001.fq.gz')]
            else:
                 # Fallback attempt: remove common suffixes
                temp_name = r1_gz
                for suffix in ['.fq.gz', '.fastq.gz']:
                    if temp_name.endswith(suffix):
                        temp_name = temp_name[:-len(suffix)]
                        break
                # Try removing common read 1 indicators
                for indicator in ['_1', '_R1', '_R1_001']:
                     if temp_name.endswith(indicator):
                          base_name = temp_name[:-len(indicator)]
                          break
                if base_name is None: # If still no match, use filename without final extension
                    base_name = os.path.splitext(os.path.splitext(r1_gz)[0])[0]
                    print(f"Warning: Could not reliably determine base name from '{r1_gz}'. Using '{base_name}'. Assumed R2 is '{base_name}_2.fq.gz' or similar.", file=sys.stderr)

            # --- Find Corresponding R2 File ---
            r2_gz = None
            possible_r2_suffixes = ['_2.fq.gz', '_R2.fq.gz', '_R2_001.fq.gz']
            for r2_suffix in possible_r2_suffixes:
                potential_r2 = f"{base_name}{r2_suffix}"
                if os.path.exists(potential_r2):
                    r2_gz = potential_r2
                    break

            if r2_gz is None:
                print(f"Warning: Cannot find corresponding R2 file for R1 file '{r1_gz}' (base name: {base_name}). Skipping sample.")
                continue

            print(f"Processing sample: {base_name} (R1: {r1_gz}, R2: {r2_gz})")
            sample_output_dir = os.path.join(MATCHED_BASE_DIR, base_name)
            os.makedirs(sample_output_dir, exist_ok=True)

            match_counts = defaultdict(int)
            pairs_processed = 0
            output_files = {}

            try:
                for i in PATTERNS_DATA.keys():
                    output_filename = os.path.join(sample_output_dir, f"pattern_{i}_matches.fq")
                    output_files[i] = open(output_filename, 'w')

                with gzip.open(r1_gz, 'rt', encoding='utf-8') as f_r1, \
                     gzip.open(r2_gz, 'rt', encoding='utf-8') as f_r2:

                    r1_records = read_fastq_records(f_r1)
                    r2_records = read_fastq_records(f_r2)
                    record_pairs = zip(r1_records, r2_records)

                    if DEFAULT_READS is not None:
                        record_pairs = itertools.islice(record_pairs, DEFAULT_READS)
                        print(f"  Analyzing first {reads_analyzed_desc} read pairs...")
                    else:
                        print(f"  Analyzing all read pairs...")


                    for r1_rec, r2_rec in record_pairs:
                        pairs_processed += 1
                        r1_header, r1_seq, _, r1_qual = r1_rec
                        r2_header, r2_seq, _, r2_qual = r2_rec

                        for pattern_idx, fragments in PATTERNS_DATA.items():
                            # Check if defining fragments are present (initial filter)
                            fragments_found = True
                            for frag in fragments:
                                if frag not in r1_seq and frag not in r2_seq:
                                    fragments_found = False
                                    break

                            if fragments_found:
                                # If broadly matched, attempt UMI extraction using precise flanks
                                umi = 'N' * UMI_LEN # Default

                                if pattern_idx in UMI_FLANKING_SEQUENCES:
                                    f_before, f_after = UMI_FLANKING_SEQUENCES[pattern_idx]

                                    # --- Refined Read Priority Logic ---
                                    # Based on typical locations of these fragments:
                                    # Pat 1 (T7...UMI...TSO): Likely R1
                                    # Pat 2 (BC_end...UMI...P5): Likely R2
                                    # Pat 3 (BC_end...UMI...CCC): Could be R1 based on full example? Revert to R1 priority.
                                    # Pat 4 (GGG...UMI...TSO): Likely R2 based on full example? Revert to R2 priority.
                                    prio = 'R1' # Default
                                    if pattern_idx in [2, 4]:
                                        prio = 'R2'
                                    # Allow Pattern 3 to default to R1

                                    umi = get_umi_from_pair(r1_seq, r2_seq, f_before, f_after, priority=prio, umi_len=UMI_LEN)
                                else:
                                    # This shouldn't happen with current setup
                                    print(f"  Warning: No UMI flanking sequences defined for matched pattern {pattern_idx}. UMI set to Ns.", file=sys.stderr)

                                # --- Prepare Output ---
                                # Only write if UMI was potentially found (or pattern matched)
                                # If we only want reads where UMI extraction works, add: if umi != 'N' * UMI_LEN:
                                clean_header = re.sub(r'/[12]$', '', r1_header.split()[0])
                                umi_tag = f"UMI:{umi}"
                                new_header = f"{clean_header} {umi_tag}"

                                output_files[pattern_idx].write(f"{new_header}\n")
                                output_files[pattern_idx].write(f"{r1_seq}\t{r2_seq}\n")
                                output_files[pattern_idx].write("+\n")
                                output_files[pattern_idx].write(f"{r1_qual}\t{r2_qual}\n")

                                # Increment match count *after* successful processing/writing
                                match_counts[pattern_idx] += 1
                                # Break after finding the first matching pattern and processing it
                                break


                    total_analyzed_for_sample = pairs_processed
                    print(f"  Finished processing {total_analyzed_for_sample} pairs for sample {base_name}.")
                    total_matched_for_sample = sum(match_counts.values())
                    print(f"  Total pairs matched to any pattern: {total_matched_for_sample}")
                    for pattern_idx in PATTERNS_DATA.keys():
                         matched = match_counts[pattern_idx]
                         percentage = (matched / total_analyzed_for_sample * 100) if total_analyzed_for_sample > 0 else 0.0
                         print(f"    Pattern {pattern_idx}: Matched {matched} pairs ({percentage:.4f}%)")
                         report_data.append({
                             "sample": base_name,
                             "total": total_analyzed_for_sample,
                             "matched": matched,
                             "percentage": percentage,
                             "pattern": pattern_idx
                         })

            except FileNotFoundError as e:
                 print(f"Error: Input/Output file error for sample {base_name}: {e}", file=sys.stderr)
            except Exception as e:
                 print(f"An unexpected error occurred while processing sample {base_name}: {e}", file=sys.stderr)
            finally:
                for f_handle in output_files.values():
                    if f_handle and not f_handle.closed:
                         f_handle.close()

        # --- Write Final Report ---
        with open(REPORT_FILE, 'a') as rp:
            # Sort report data for readability
            for entry in sorted(report_data, key=lambda x: (x['sample'], x['pattern'])):
                 rp.write(f"{entry['sample']}\t{entry['total']}\t{entry['matched']}\t{entry['percentage']:.4f}%\tPattern {entry['pattern']}\n")
                 # Add a separator line after the last pattern for each sample
                 if entry['pattern'] == max(PATTERNS_DATA.keys()):
                      rp.write("-" * 90 + "\n")

        print("\nAnalysis finished!")
        print(f"Statistics report saved to: {REPORT_FILE}")
        print(f"Matched read pairs (with extracted UMIs based on FINAL flanks) saved to: {MATCHED_BASE_DIR}")
        print("Output files are in a 4-line-per-pair format (Header+UMI / R1_Seq<TAB>R2_Seq / + / R1_Qual<TAB>R2_Qual).")

    finally:
        pass # No temp dir cleanup needed


if __name__ == "__main__":
    main()
