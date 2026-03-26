"""
富集分析模块：计算各轮次 Fold Change
"""
import os
import re
import pandas as pd
from functools import reduce


def compute_fold_change_by_round(parent_folder, output_csv,
                                 epsilon=1e-9, outer_join=False,
                                 target_length=None, fre_file_name="fre.csv"):
    """
    读取各轮次子目录中的 fre.csv，计算相对第一轮的 Fold Change。

    参数:
    - parent_folder: 父目录，每个以 R 开头的子目录表示一轮
    - output_csv: 输出 CSV 路径
    - epsilon: 防除零小量（当前未使用，保留接口）
    - outer_join: 是否外连接（True 保留所有序列，False 仅保留每轮都出现的序列）
    - target_length: 若指定，只保留该长度的序列
    - fre_file_name: 各轮次频率文件名（默认 fre.csv）
    """
    def round_key(name):
        m = re.search(r'\d+', name)
        return int(m.group()) if m else float('inf')

    rounds = sorted(
        [d for d in os.listdir(parent_folder)
         if os.path.isdir(os.path.join(parent_folder, d)) and d.startswith('R')],
        key=round_key
    )

    if not rounds:
        print("No round folders found.")
        return

    round_freqs = []
    for r in rounds:
        csv_path = os.path.join(parent_folder, r, fre_file_name)
        if not os.path.exists(csv_path):
            print(f"Warning: {r} has no {fre_file_name}, skipped.")
            continue

        df = pd.read_csv(csv_path)
        if 'Sequence' not in df.columns or 'Count' not in df.columns:
            print(f"Warning: {r}/{fre_file_name} 格式不正确，需包含 'Sequence' 和 'Count'")
            continue

        if target_length is not None:
            before = len(df)
            df = df[df['Sequence'].apply(lambda x: isinstance(x, str) and len(x) == target_length)]
            print(f"{r}: {before} 条 -> {len(df)} 条 (保留长度={target_length})")

        df = df[['Sequence', 'Count']].rename(columns={'Count': f'Freq_{r}'})
        round_freqs.append(df)

    if not round_freqs:
        print("No valid fre.csv files found.")
        return

    how = 'outer' if outer_join else 'inner'
    merged = reduce(lambda left, right: pd.merge(left, right, on='Sequence', how=how), round_freqs)
    merged.fillna(0, inplace=True)

    fc_df = merged.copy()
    base_col = f'Freq_{rounds[0]}'
    for i in range(1, len(rounds)):
        curr_col = f'Freq_{rounds[i]}'
        fc_df[f'FC_{rounds[i]}'] = merged[curr_col] / merged[base_col]

    fc_df.to_csv(output_csv, index=False)
    print(f"Fold change CSV saved to: {output_csv}")
