"""
主流程脚本：针对 liuyang_PTK7_yingxingduizhao 数据集
"""
import os
import glob
import shutil

from selex.preprocessing import run_pandaseq
from selex.extraction import process_sequences_regex_only
from selex.analysis import compute_fold_change_by_round
from selex.structure import add_similarity_column_and_auc  # 取消注释以启用结构分析

# ── 配置 ────────────────────────────────────────────────────────────────────
DATA_DIR = "/home/hanlab/2022_zzm/3DSF-SELEX/data/PTK7"
PATTERN  = r"TTCAGCACTCCACGATCTAAC([ACGT]{9}GCCG[ACGT]{12})CGGTTAGATCGTGCTACCGTGAA"
ROUNDS   = [f"R{i}" for i in range(1, 11)]

P7_ADAPTER = None
P5_ADAPTER = None
# ────────────────────────────────────────────────────────────────────────────


def main():
    output_dir = os.path.join(DATA_DIR, "output")
    os.makedirs(output_dir, exist_ok=True)

    fre_outputs = {}
    for R_path in glob.glob(os.path.join(DATA_DIR, "*")):
        round_name = os.path.basename(R_path)
        if round_name not in ROUNDS:
            continue

        # 双端拼接
        merged_fa = os.path.join(R_path, "merged.fa")
        if not os.path.exists(merged_fa):
            f_fastq = glob.glob(os.path.join(R_path, "*1.f*q"))
            r_fastq = glob.glob(os.path.join(R_path, "*2.f*q"))
            merged_fa = run_pandaseq(
                f_fastq[0], r_fastq[0], merged_fa,
                p5_adapter=P5_ADAPTER, p7_adapter=P7_ADAPTER
            )

        # 序列提取
        fre_output = os.path.join(R_path, "fre.csv")
        fre_outputs[round_name] = process_sequences_regex_only(merged_fa, PATTERN, fre_output)

        # 备份到 output 目录
        shutil.copy(fre_output, os.path.join(output_dir, f"{round_name}_fre.csv"))

    # Fold Change 计算
    compute_fold_change_by_round(DATA_DIR, "fold_change.csv")

    # 结构分析（可选）
    add_similarity_column_and_auc("fold_change.csv", "fold_change_with_similarity.csv")


if __name__ == "__main__":
    main()
