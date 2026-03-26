"""
reads预处理模块：去接头、质量筛选、双端拼接
"""
import subprocess


def run_pandaseq(forward_fastq, reverse_fastq, output_fasta,
                 forward_primer=None, reverse_primer=None,
                 threads=16,
                 p7_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                 p5_adapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                 quality_cutoff=20,
                 min_length=50,
                 use_cutadapt=False):
    """
    预处理 + 拼接测序 reads，可选是否使用 Cutadapt 去接头和质量筛选。

    步骤：
      1. 可选 Cutadapt 去除 Illumina 接头 + 质量筛选
      2. 使用 Pandaseq 拼接双端 reads，生成 FASTA

    参数:
    - forward_fastq: R1 文件路径
    - reverse_fastq: R2 文件路径
    - output_fasta: 输出 FASTA 文件路径
    - forward_primer, reverse_primer: 可选引物序列
    - p7_adapter, p5_adapter: 接头序列
    - quality_cutoff: 最低 Q 值 (默认 20)
    - min_length: 最短保留长度 (默认 50)
    - threads: 线程数
    - use_cutadapt: 是否使用 Cutadapt (默认 False)

    返回:
    - output_fasta 文件路径，失败返回 None
    """
    if use_cutadapt:
        cleaned_fwd = forward_fastq.replace(".fq", "_clean.fq")
        cleaned_rev = reverse_fastq.replace(".fq", "_clean.fq")

        print("开始去接头和质量筛选...")
        cutadapt_cmd = [
            "cutadapt",
            "-a", p7_adapter,
            "-A", p5_adapter,
            "-q", str(quality_cutoff),
            "-m", str(min_length),
            "-j", str(threads),
            "--max-n", "0",
            "-o", cleaned_fwd,
            "-p", cleaned_rev,
            forward_fastq, reverse_fastq
        ]

        cutadapt_proc = subprocess.run(cutadapt_cmd, stderr=subprocess.PIPE, text=True)
        if cutadapt_proc.returncode != 0:
            print("Cutadapt 处理失败：")
            print(cutadapt_proc.stderr)
            return None
        print("Cutadapt 处理完成")
    else:
        cleaned_fwd = forward_fastq
        cleaned_rev = reverse_fastq
        print("未使用 Cutadapt，直接使用原始 FASTQ 文件进行拼接")

    print("开始 Pandaseq 拼接双端 reads...")
    cmd = [
        "pandaseq",
        "-f", cleaned_fwd,
        "-r", cleaned_rev,
        "-T", str(threads),
        "-w", output_fasta
    ]

    if forward_primer:
        cmd.extend(["-p", forward_primer])
    if reverse_primer:
        cmd.extend(["-q", reverse_primer])

    process = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)

    if process.returncode != 0:
        print("Pandaseq 执行失败：")
        print(process.stderr)
        return None

    print(f"Pandaseq 拼接完成，结果保存在: {output_fasta}")
    return output_fasta
