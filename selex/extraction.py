"""
序列提取模块：从 FASTA/FASTQ 中提取随机区域，统计频次
"""
import os
import csv
import re
from collections import Counter


def process_sequences(input_file, f_linker, r_linker, output):
    """
    提取 FASTA 或 FASTQ 序列中 f_linker 和 r_linker 之间的随机序列，并按出现频次保存 CSV。

    参数:
    - input_file: 输入 FASTA 或 FASTQ 文件路径
    - f_linker: 正向引物或接头序列
    - r_linker: 反向引物或接头序列
    - output: 输出 CSV 文件路径
    """
    f_linker = f_linker.upper()
    r_linker = r_linker.upper()
    seq_counter = Counter()

    file_ext = os.path.splitext(input_file)[1].lower()
    is_fastq = file_ext in [".fq", ".fastq"]

    def _extract(seq):
        f_pos = seq.find(f_linker)
        r_pos = seq.find(r_linker)
        if f_pos != -1 and r_pos != -1 and f_pos + len(f_linker) < r_pos:
            random_seq = seq[f_pos + len(f_linker):r_pos]
            if random_seq and len(random_seq) == 25 and random_seq[9:13] == "GCCG":
                seq_counter[random_seq] += 1

    with open(input_file, "r") as f:
        if is_fastq:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip().upper()
                f.readline()   # +
                f.readline()   # quality
                _extract(seq)
        else:
            seq = ""
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq:
                        _extract(seq)
                    seq = ""
                else:
                    seq += line.upper()
            if seq:
                _extract(seq)

    sorted_seqs = seq_counter.most_common()
    with open(output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Sequence", "Count"])
        for seq, count in sorted_seqs:
            writer.writerow([seq, count])

    print(f"Finished: {os.path.basename(input_file)} -> {output}")
    return output


def process_sequences_regex_only(input_file, pattern, output):
    """
    使用正则表达式提取 FASTA/FASTQ 序列中的随机片段，按频次保存 CSV。
    pattern 须包含一个捕获组，提取其中的核心序列；若无捕获组则取整个匹配。

    参数:
    - input_file: 输入 FASTA 或 FASTQ 文件路径
    - pattern: 正则表达式字符串
    - output: 输出 CSV 文件路径
    """
    regex = re.compile(pattern)
    seq_counter = Counter()

    file_ext = os.path.splitext(input_file)[1].lower()
    is_fastq = file_ext in [".fq", ".fastq"]

    def _get_captured(match):
        """有捕获组取 group(1)，否则取 group(0)（整个匹配）。"""
        return match.group(1).upper() if match.lastindex else match.group(0).upper()

    with open(input_file, "r") as f:
        if is_fastq:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                match = regex.search(seq)
                if match:
                    seq_counter[_get_captured(match)] += 1
        else:
            seq = ""
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq:
                        match = regex.search(seq)
                        if match:
                            seq_counter[_get_captured(match)] += 1
                        seq = ""
                else:
                    seq += line.upper()
            if seq:
                for match in regex.finditer(seq):
                    seq_counter[_get_captured(match)] += 1

    with open(output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Sequence", "Count"])
        for seq, count in seq_counter.most_common():
            writer.writerow([seq, count])

    print(f"Finished: {os.path.basename(input_file)} -> {output}")
    return output
