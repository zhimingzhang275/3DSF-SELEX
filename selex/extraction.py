"""
序列提取模块：从 FASTA/FASTQ 中提取随机区域，统计频次
"""
import os
import csv
import re
from collections import Counter


def _get_captured(match):
    return match.group(1).upper() if match.lastindex else match.group(0).upper()


def _iter_fasta_seqs(f):
    seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq:
                yield seq
            seq = ""
        else:
            seq += line.upper()
    if seq:
        yield seq


def _iter_fastq_seqs(f):
    while True:
        header = f.readline()
        if not header:
            break
        seq = f.readline().strip().upper()
        f.readline()  # +
        f.readline()  # quality
        yield seq


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
    is_fastq = file_ext in (".fq", ".fastq")

    with open(input_file, "r") as f:
        iter_seqs = _iter_fastq_seqs(f) if is_fastq else _iter_fasta_seqs(f)
        for seq in iter_seqs:
            match = regex.search(seq)
            if match:
                seq_counter[_get_captured(match)] += 1

    with open(output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Sequence", "Count"])
        for seq, count in seq_counter.most_common():
            writer.writerow([seq, count])

    print(f"Finished: {os.path.basename(input_file)} -> {output}")
    return output
