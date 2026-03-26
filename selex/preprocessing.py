"""
reads预处理模块：双端拼接
"""
import subprocess


def run_pandaseq(forward_fastq, reverse_fastq, output_fasta, threads=16):
    """
    使用 Pandaseq 拼接双端 reads，生成 FASTA。

    参数:
    - forward_fastq: R1 文件路径
    - reverse_fastq: R2 文件路径
    - output_fasta: 输出 FASTA 文件路径
    - threads: 线程数

    返回:
    - output_fasta 文件路径，失败返回 None
    """
    print("开始 Pandaseq 拼接双端 reads...")
    cmd = [
        "pandaseq",
        "-f", forward_fastq,
        "-r", reverse_fastq,
        "-T", str(threads),
        "-w", output_fasta
    ]

    process = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)

    if process.returncode != 0:
        print("Pandaseq 执行失败：")
        print(process.stderr)
        return None

    print(f"Pandaseq 拼接完成，结果保存在: {output_fasta}")
    return output_fasta
