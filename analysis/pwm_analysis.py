"""
PWM 计算与可视化脚本
用法:
    python pwm_analysis.py input.fa [--pattern REGEX] [--output_dir DIR] [--top N] [--title TITLE]

流程:
    1. 从 FASTA 用正则提取捕获组序列（如指定 --pattern）
    2. 计算 PFM / PPM / PWM
    3. 输出 CSV 矩阵 + sequence logo + PWM 热图

输出文件（前缀 = output_dir/title）:
    - _extracted.fa      : 提取出的等长序列（FASTA）
    - _pwm_counts.csv    : 计数矩阵 (PFM)
    - _pwm_freq.csv      : 频率矩阵 (PPM)
    - _pwm_score.csv     : 权重矩阵 (PWM, log2 odds)
    - _logo.png          : sequence logo（bits）
    - _heatmap.png       : PWM 热图
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logomaker
from Bio import SeqIO

# ── 正则提取 ──────────────────────────────────────────────────────────────────

def extract_by_pattern(fasta_path, pattern, top_n=None):
    """从 FASTA 每条序列中用正则提取第一个捕获组，返回提取到的序列列表。"""
    regex = re.compile(pattern)
    extracted = []
    total = 0
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper()
        total += 1
        m = regex.search(seq)
        if m:
            # 有捕获组则取第 1 组，否则取整个匹配
            hit = m.group(1) if m.lastindex and m.lastindex >= 1 else m.group(0)
            extracted.append(hit)
        if top_n and total >= top_n:
            break
    print(f"[INFO] 扫描序列数: {total}，正则匹配成功: {len(extracted)}")
    if not extracted:
        sys.exit("[ERROR] 正则未匹配到任何序列，请检查 --pattern 参数")
    return extracted


def save_extracted_fasta(seqs, out_path):
    """将提取的序列写成简单 FASTA 文件。"""
    with open(out_path, "w") as f:
        for i, s in enumerate(seqs):
            f.write(f">seq_{i+1}\n{s}\n")
    print(f"[OK] 提取序列已保存至 {out_path}（共 {len(seqs)} 条）")


# ── 读取 FASTA（无正则时用） ──────────────────────────────────────────────────

def read_fasta(fasta_path, top_n=None):
    seqs = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seqs.append(str(rec.seq).upper())
        if top_n and len(seqs) >= top_n:
            break
    if not seqs:
        sys.exit(f"[ERROR] {fasta_path} 中未读取到任何序列")
    return seqs


# ── 序列预处理 ────────────────────────────────────────────────────────────────

def validate_and_trim(seqs):
    """过滤非 ACGT 序列，截断到最短长度。"""
    valid = [s for s in seqs if set(s).issubset(set("ACGT"))]
    if len(valid) < len(seqs):
        print(f"[WARN] 过滤掉 {len(seqs) - len(valid)} 条含非 ACGT 字符的序列")
    if not valid:
        sys.exit("[ERROR] 过滤后无有效序列")
    lengths = set(len(s) for s in valid)
    if len(lengths) > 1:
        min_len = min(lengths)
        print(f"[WARN] 序列长度不一致（{min(lengths)}~{max(lengths)}），截断至 {min_len}")
        valid = [s[:min_len] for s in valid]
    return valid, len(valid[0])


# ── PWM 计算 ──────────────────────────────────────────────────────────────────

def compute_pfm(seqs, length):
    counts = np.zeros((length, 4), dtype=int)
    base_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    for seq in seqs:
        for i, base in enumerate(seq):
            counts[i, base_idx[base]] += 1
    return pd.DataFrame(counts, columns=["A", "C", "G", "T"])


def pfm_to_ppm(pfm, pseudocount=0.5):
    pfm_pseudo = pfm + pseudocount
    return pfm_pseudo.div(pfm_pseudo.sum(axis=1), axis=0)


def ppm_to_pwm(ppm, bg=None):
    if bg is None:
        bg = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    bg_arr = np.array([bg["A"], bg["C"], bg["G"], bg["T"]])
    return pd.DataFrame(np.log2(ppm.values / bg_arr), columns=["A", "C", "G", "T"])


def ppm_to_ic(ppm):
    entropy = -(ppm * np.log2(ppm + 1e-9)).sum(axis=1)
    ic = 2 - entropy
    return ppm.multiply(ic, axis=0)


# ── 可视化 ────────────────────────────────────────────────────────────────────

def plot_logo(ic_matrix, title, out_path):
    fig, ax = plt.subplots(figsize=(max(8, len(ic_matrix) * 0.45), 3))
    logomaker.Logo(ic_matrix, ax=ax, color_scheme="classic")
    ax.set_xlabel("Position")
    ax.set_ylabel("Information Content (bits)")
    ax.set_title(f"Sequence Logo — {title}")
    ax.set_ylim(0, 2)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[OK] logo 保存至 {out_path}")


def plot_heatmap(pwm, title, out_path):
    fig, ax = plt.subplots(figsize=(max(10, len(pwm) * 0.45), 3))
    im = ax.imshow(pwm.T.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(pwm)))
    ax.set_xticklabels([str(i + 1) for i in range(len(pwm))], fontsize=8)
    ax.set_yticks(range(4))
    ax.set_yticklabels(["A", "C", "G", "T"])
    ax.set_xlabel("Position")
    ax.set_title(f"PWM Heatmap (log2 odds) — {title}")
    plt.colorbar(im, ax=ax, label="log2 odds")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[OK] heatmap 保存至 {out_path}")


# ── 主流程 ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="正则提取序列 → 计算 PWM → 可视化")
    parser.add_argument("fasta",help="输入 FASTA 文件路径")
    parser.add_argument("--pattern", default='TTCAGCACTCCACGATCTAAC([ACGT]{9}GCCG[ACGT]{12})CGGTTAGATCGTGCTACCGTGAA',
                        help="正则表达式，用第一个捕获组 () 提取目标序列。"
                             "例: 'TTCAGCACTCCACGATCTAAC([ACGT]{9}GCCG[ACGT]{12})CGGTTAGATCGTGCTACCGTGAA'")
    parser.add_argument("--output_dir", default=None,
                        help="输出目录（默认与 fasta 同目录）")
    parser.add_argument("--top", type=int, default=None,
                        help="只扫描前 N 条序列（默认全部）")
    parser.add_argument("--title", default=None,
                        help="图标题前缀（默认使用文件名）")
    parser.add_argument("--pseudocount", type=float, default=0.5,
                        help="伪计数（默认 0.5）")
    args = parser.parse_args()

    if not os.path.isfile(args.fasta):
        sys.exit(f"[ERROR] 文件不存在: {args.fasta}")

    out_dir = args.output_dir or os.path.join(os.path.dirname(os.path.abspath(args.fasta)), "pwm_out")
    os.makedirs(out_dir, exist_ok=True)
    title = args.title or os.path.splitext(os.path.basename(args.fasta))[0]
    prefix = os.path.join(out_dir, title.replace(" ", "_"))

    # 1. 提取序列
    if args.pattern:
        print(f"[INFO] 使用正则提取: {args.pattern}")
        seqs = extract_by_pattern(args.fasta, args.pattern, top_n=args.top)
        save_extracted_fasta(seqs, f"{prefix}_extracted.fa")
    else:
        print(f"[INFO] 未指定 --pattern，直接读取全部序列")
        seqs = read_fasta(args.fasta, top_n=args.top)

    # 2. 预处理
    seqs, length = validate_and_trim(seqs)
    print(f"[INFO] 用于 PWM 的序列数: {len(seqs)}，长度: {length}")

    # 3. 计算矩阵
    pfm = compute_pfm(seqs, length)
    ppm = pfm_to_ppm(pfm, pseudocount=args.pseudocount)
    pwm = ppm_to_pwm(ppm)
    ic_matrix = ppm_to_ic(ppm)

    # 4. 保存 CSV（index 从 1 开始）
    for df, suffix in [(pfm, "_pwm_counts.csv"),
                       (ppm, "_pwm_freq.csv"),
                       (pwm, "_pwm_score.csv")]:
        df_out = df.copy()
        df_out.index = range(1, length + 1)
        fmt = "%.4f" if suffix != "_pwm_counts.csv" else None
        df_out.to_csv(f"{prefix}{suffix}", index_label="Position",
                      float_format=fmt)
    print(f"[OK] CSV 矩阵已保存（前缀: {prefix}）")

    # 5. 可视化
    plot_logo(ic_matrix, title, f"{prefix}_logo.png")
    plot_heatmap(pwm, title, f"{prefix}_heatmap.png")

    # 6. 摘要
    print(f"\n── PWM (log2 odds)，共 {length} 位 ──")
    pwm_show = pwm.copy()
    pwm_show.index = range(1, length + 1)
    print(pwm_show.to_string(float_format="%.3f"))


if __name__ == "__main__":
    main()
