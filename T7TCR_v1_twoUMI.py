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

# Define patterns and their fragments
PATTERNS_DATA = {
    1: ("TAATACGACTCACTATAGGG", "CGAAACATCGGCCACCCC", "CAACCCTGCGAC"),
    2: ("CGCAATGAAGTCGCAGGGTTG", "GGGGTGGCCGATGTTTCG", "CCCTATAGTGAGTCGTATTA"),
    3: ("TAATACGACTCACTATAGGGGTGGCCGATGTTTCG", "CCC", "CAACCCTGCGACTTCA"),
    4: ("CAATGAAGTCGCAGGGTTG", "GGG", "CGAAACATCGGCCACCCCTATAGTGAGTCGTATTA")
}

# Define where UMIs are *expected* to be found relative to fragments
UMI_LOCATIONS = {
    1: (('R1', 0, 1), ('R2', 1, 2)),
    2: (('R1', 0, 1), ('R2', 1, 2)),
    3: (('R1', 0, 1), ('R2', 1, 2)),
    4: (('R1', 0, 1), ('R2', 1, 2)),
}

HEAD_LABEL = f"head{DEFAULT_READS}" if DEFAULT_READS else "all"
TEMP_DIR = f"{HEAD_LABEL}_temp_fq_py"
MATCHED_BASE_DIR = f"{HEAD_LABEL}_matched_reads_umi_format_py"
REPORT_FILE = f"{HEAD_LABEL}_analysis_report_merged_umi_py.txt"
UMI_LEN = 8

# --- Helper Functions ---

def read_fastq_records(file_handle):
    """Generator to yield 4-line FASTQ records from a file handle."""
    while True:
        lines = list(itertools.islice(file_handle, 4))
        if not lines or len(lines) < 4:
            break
        yield tuple(line.decode('utf-8').strip() if isinstance(line, bytes) else line.strip() for line in lines)

def extract_umi(seq, frag_before, frag_after, umi_len=8):
    """
    Extracts a UMI of specified length between two flanking fragments using regex.
    Returns UMI sequence or 'N'*umi_len if not found.
    """
    try:
        f_before_esc = re.escape(frag_before)
        f_after_esc = re.escape(frag_after)
        regex = re.compile(f"{f_before_esc}(.{{{umi_len}}}){f_after_esc}")
        match = regex.search(seq)
        if match:
            return match.group(1)
    except re.error as e:
        print(f"  Regex error extracting UMI between '{frag_before}' and '{frag_after}': {e}", file=sys.stderr)
    return 'N' * umi_len

def get_umi_from_pair(r1_seq, r2_seq, frag_before, frag_after, priority='R1', umi_len=8):
    """
    Extracts UMI based on priority/fallback logic.
    priority: 'R1' (try R1 first), 'R2' (try R2 first)
    """
    primary_seq = r1_seq if priority == 'R1' else r2_seq
    fallback_seq = r2_seq if priority == 'R1' else r1_seq

    umi = extract_umi(primary_seq, frag_before, frag_after, umi_len)
    # Use 'N'*umi_len for comparison, not just "N"
    if umi == 'N' * umi_len:
        umi = extract_umi(fallback_seq, frag_before, frag_after, umi_len)
    return umi


# --- Main Processing Logic ---

def main():
    os.makedirs(TEMP_DIR, exist_ok=True)
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
            rp.write("\n=== Sample Statistics ===\n\n")
            rp.write("Sample\tTotal_Pairs_Analyzed\tMatched_Pairs\tPercentage (% of pairs)\tPattern_Index\n")
            rp.write("-" * 90 + "\n")

        report_data = []
        r1_files = sorted(glob.glob('*_1.fq.gz'))
        if not r1_files:
             print("No *_1.fq.gz files found in the current directory.")
             return

        for r1_gz in r1_files:
            # --- COMPATIBILITY FIX ---
            suffix = '_1.fq.gz'
            if r1_gz.endswith(suffix):
                base_name = r1_gz[:-len(suffix)]
            else:
                base_name = r1_gz
                print(f"警告: 文件名 '{r1_gz}' 不以 '{suffix}' 结尾。", file=sys.stderr)
            # --- END FIX ---

            r2_gz = f"{base_name}_2.fq.gz"

            if not os.path.exists(r2_gz):
                print(f"警告: 找不到对应的 R2 文件 ({r2_gz}) for {r1_gz}。跳过样本 {base_name}。")
                continue

            print(f"处理样本: {base_name}")
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

                    print(f"  分析 {reads_analyzed_desc} read 对...")
                    for r1_rec, r2_rec in record_pairs:
                        pairs_processed += 1
                        r1_header, r1_seq, _, r1_qual = r1_rec
                        r2_header, r2_seq, _, r2_qual = r2_rec

                        for pattern_idx, fragments in PATTERNS_DATA.items():
                            if len(fragments) < 3:
                                # Removed redundant warning print here, user likely knows
                                continue

                            fragments_found = True
                            for frag in fragments:
                                if frag not in r1_seq and frag not in r2_seq:
                                    fragments_found = False
                                    break

                            if fragments_found:
                                match_counts[pattern_idx] += 1
                                frag1, frag2, frag3 = fragments
                                umi1 = 'N' * UMI_LEN
                                umi2 = 'N' * UMI_LEN

                                if pattern_idx in UMI_LOCATIONS:
                                    (umi1_loc, f1_idx, f2_idx), (umi2_loc, f3_idx, f4_idx) = UMI_LOCATIONS[pattern_idx]
                                    # Check if indices are valid for the fragments list
                                    if max(f1_idx, f2_idx, f3_idx, f4_idx) < len(fragments):
                                        f_before1 = fragments[f1_idx]
                                        f_after1 = fragments[f2_idx]
                                        f_before2 = fragments[f3_idx]
                                        f_after2 = fragments[f4_idx]
                                        prio1 = 'R1' if 'R1' in umi1_loc else 'R2'
                                        prio2 = 'R2' if 'R2' in umi2_loc else 'R1'
                                        umi1 = get_umi_from_pair(r1_seq, r2_seq, f_before1, f_after1, priority=prio1, umi_len=UMI_LEN)
                                        umi2 = get_umi_from_pair(r1_seq, r2_seq, f_before2, f_after2, priority=prio2, umi_len=UMI_LEN)
                                    else:
                                         print(f"  警告: 模式 {pattern_idx} 的 UMI_LOCATIONS 索引超出范围。", file=sys.stderr)
                                else:
                                    umi1 = get_umi_from_pair(r1_seq, r2_seq, frag1, frag2, priority='R1', umi_len=UMI_LEN)
                                    umi2 = get_umi_from_pair(r1_seq, r2_seq, frag2, frag3, priority='R2', umi_len=UMI_LEN)

                                clean_header = re.sub(r'/[12]$', '', r1_header)
                                umi_tag = f"UMI:{umi1}_{umi2}"
                                new_header = f"{clean_header} {umi_tag}"

                                output_files[pattern_idx].write(f"{new_header}\n")
                                output_files[pattern_idx].write(f"{r1_seq}\t{r2_seq}\n")
                                output_files[pattern_idx].write("+\n")
                                output_files[pattern_idx].write(f"{r1_qual}\t{r2_qual}\n")

                    total_analyzed_for_sample = pairs_processed
                    print(f"  处理完成 {total_analyzed_for_sample} 对 reads.")
                    for pattern_idx in PATTERNS_DATA.keys():
                         matched = match_counts[pattern_idx]
                         percentage = (matched / total_analyzed_for_sample * 100) if total_analyzed_for_sample > 0 else 0.0
                         print(f"    模式 {pattern_idx}: 匹配 {matched} 对 ({percentage:.4f}%)")
                         report_data.append({
                             "sample": base_name,
                             "total": total_analyzed_for_sample,
                             "matched": matched,
                             "percentage": percentage,
                             "pattern": pattern_idx
                         })

            except FileNotFoundError as e:
                 print(f"错误: 输入/输出文件错误 for sample {base_name}: {e}", file=sys.stderr)
            except Exception as e:
                 print(f"处理样本 {base_name} 时发生错误: {e}", file=sys.stderr)
                 # Optionally re-raise or handle more specifically
                 # raise
            finally:
                for f_handle in output_files.values():
                    if f_handle and not f_handle.closed:
                         f_handle.close()

        with open(REPORT_FILE, 'a') as rp:
            for entry in sorted(report_data, key=lambda x: (x['sample'], x['pattern'])):
                 rp.write(f"{entry['sample']}\t{entry['total']}\t{entry['matched']}\t{entry['percentage']:.4f}%\tPattern {entry['pattern']}\n")
                 if entry['pattern'] == max(PATTERNS_DATA.keys()):
                      rp.write("-" * 90 + "\n")

        print("\n分析完成！")
        print(f"统计结果已保存到：{REPORT_FILE}")
        print(f"每个样本匹配每个模式、并提取了UMI的Read对已保存到 {MATCHED_BASE_DIR} 下。")
        print("输出文件是 4 行合并格式（非标准 FASTQ），Header包含UMI信息。")

    finally:
        if os.path.exists(TEMP_DIR):
             print(f"执行清理: 删除临时目录 {TEMP_DIR}")
             # Add error handling for directory removal if needed
             try:
                 shutil.rmtree(TEMP_DIR)
             except OSError as e:
                 print(f"错误: 无法删除临时目录 {TEMP_DIR}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
