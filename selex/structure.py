"""
结构分析模块：预测 DNA 二级结构，计算支架相似度和富集 AUC
"""
import pandas as pd
import numpy as np
from scipy.integrate import trapezoid
import nupack


TARGET_SCAFFOLD = ".(((((.*********.(((************))))))))"


def check_conditions(seq, scaffold=TARGET_SCAFFOLD):
    """
    预测序列的 MFE 二级结构，并计算与支架模板的相似度。

    参数:
    - seq: DNA 序列字符串
    - scaffold: 支架模板点括号字符串，'*' 表示不限制位置

    返回:
    - (similarity, predicted_structure): 相似度 (0~1)，预测结构字符串
    """
    model = nupack.Model(material="dna", celsius=25, sodium=0.14, magnesium=0.005)
    result = nupack.mfe(strands=[seq], model=model)
    ss = str(result[0].structure)

    if len(ss) != len(scaffold):
        raise ValueError(f"二级结构长度 {len(ss)} 与支架模板长度 {len(scaffold)} 不一致")

    match_count = 0
    total_fixed = 0
    for s_char, t_char in zip(ss, scaffold):
        if t_char == '*':
            continue
        total_fixed += 1
        if s_char == t_char:
            match_count += 1

    similarity = match_count / total_fixed if total_fixed > 0 else 0
    return similarity, ss


def add_similarity_column_and_auc(fc_csv, output_csv,
                                   added_f_primer="ATCTAAC",
                                   added_r_primer="CGGTTAGA"):
    """
    为 Fold Change 表添加结构相似度列，并计算整体及分段 AUC。

    参数:
    - fc_csv: 输入 FC CSV 文件路径
    - output_csv: 输出 CSV 文件路径
    - added_f_primer: 拼接到随机序列前的引物（用于结构预测）
    - added_r_primer: 拼接到随机序列后的引物（用于结构预测）
    """
    df = pd.read_csv(fc_csv)

    similarities, predicted_structures = [], []
    for seq in df['Sequence']:
        try:
            full_seq = added_f_primer + seq + added_r_primer
            sim, ss = check_conditions(full_seq)
        except Exception as e:
            print(f"序列 {seq} 计算失败: {e}")
            sim, ss = None, None
        similarities.append(sim)
        predicted_structures.append(ss)

    df['Similarity'] = similarities
    df['Predicted_Structure'] = predicted_structures

    fc_cols = sorted(
        [c for c in df.columns if c.startswith("FC")],
        key=lambda x: int("".join(filter(str.isdigit, x)))
    )
    adj_fc_cols = []
    for col in fc_cols:
        new_col = f"Adj_{col}"
        df[new_col] = df[col] * df['Similarity']
        adj_fc_cols.append(new_col)

    # 整体 AUC
    auc_fc, auc_adj = [], []
    for _, row in df.iterrows():
        x = np.arange(len(fc_cols))
        auc_fc.append(trapezoid(row[fc_cols].values.astype(float), x))
        auc_adj.append(trapezoid(row[adj_fc_cols].values.astype(float), x))
    df['AUC_FC'] = auc_fc
    df['AUC_Adj_FC'] = auc_adj

    # 分段 AUC（从第一列开始，逐步延伸到后续每一列）
    for end_idx in range(1, len(fc_cols)):
        start_round = fc_cols[0].replace("FC_", "")
        end_round = fc_cols[end_idx].replace("FC_", "")
        seg_fc = fc_cols[:end_idx + 1]
        seg_adj = adj_fc_cols[:end_idx + 1]
        x = np.arange(end_idx + 1)
        df[f"AUC_FC_{start_round}_{end_round}"] = [
            trapezoid(row[seg_fc].astype(float), x) for _, row in df.iterrows()
        ]
        df[f"AUC_Adj_FC_{start_round}_{end_round}"] = [
            trapezoid(row[seg_adj].astype(float), x) for _, row in df.iterrows()
        ]

    df.to_csv(output_csv, index=False)
    print(f"结果保存到 {output_csv}")
