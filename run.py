"""
主流程脚本：针对 liuyang_PTK7_yingxingduizhao 数据集
"""
import os
import glob
import shutil

from selex.preprocessing import run_pandaseq
from selex.extraction import process_sequences_regex_only
from selex.analysis import compute_fold_change_by_round
from selex.structure import add_similarity_column_and_auc

# ── 配置 ────────────────────────────────────────────────────────────────────
DATA_DIR = "/home/hanlab/2022_zzm/3DSF-SELEX/data/IFN-gamma"
PATTERN  = r"TTCAGCACTCCACGATCTAAC([ACGT]{9}GCCG[ACGT]{12})CGGTTAGATCGTGCTACCGTGAA"
# ────────────────────────────────────────────────────────────────────────────


def process_round(r_path, round_name, output_dir):
    merged_fa = os.path.join(r_path, "merged.fa")
    fre_output = os.path.join(r_path, "fre.csv")

    if not os.path.exists(merged_fa):
        f_fastq = glob.glob(os.path.join(r_path, "*1.f*q"))
        r_fastq = glob.glob(os.path.join(r_path, "*2.f*q"))
        if not f_fastq or not r_fastq:
            print(f"Warning: {round_name} 缺少 FASTQ 文件，跳过")
            return None
        merged_fa = run_pandaseq(f_fastq[0], r_fastq[0], merged_fa)
        if merged_fa is None:
            return None

    process_sequences_regex_only(merged_fa, PATTERN, fre_output)
    shutil.copy(fre_output, os.path.join(output_dir, f"{round_name}_fre.csv"))
    return fre_output


def main():
    output_dir = os.path.join(DATA_DIR, "output")
    os.makedirs(output_dir, exist_ok=True)

    for r_path in sorted(glob.glob(os.path.join(DATA_DIR, "*"))):
        round_name = os.path.basename(r_path)
        process_round(r_path, round_name, output_dir)

    fc_csv = os.path.join(output_dir, "fold_change.csv")
    compute_fold_change_by_round(DATA_DIR, fc_csv)

    add_similarity_column_and_auc(fc_csv, os.path.join(output_dir, "fold_change_with_similarity.csv"))


if __name__ == "__main__":
    main()
