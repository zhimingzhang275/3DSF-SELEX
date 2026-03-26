# 3DSF-SELEX 数据分析流程

## 项目结构

```
3DSF-SELEX/
├── selex/
│   ├── __init__.py        # 包入口，统一导出所有函数
│   ├── preprocessing.py   # reads 预处理与双端拼接
│   ├── extraction.py      # 随机区域序列提取与频次统计
│   ├── analysis.py        # 各轮次 Fold Change 计算
│   └── structure.py       # DNA 二级结构预测与 AUC 分析
├── run.py                 # 主流程入口脚本
├── process.py             # 原始脚本（已归档）
└── data/                  # 原始测序数据
```

## 模块说明

### 1. `selex/preprocessing.py` — reads 预处理

函数：`run_pandaseq`

- 可选调用 **Cutadapt** 去除 Illumina 接头（P5/P7）并按质量阈值过滤
- 调用 **Pandaseq** 将双端 FASTQ（R1/R2）拼接为单条 FASTA
- 若已存在 `merged.fa`，主流程会跳过此步骤

### 2. `selex/extraction.py` — 序列提取

函数：`process_sequences` / `process_sequences_regex_only`

- 支持 FASTA 和 FASTQ 两种格式
- `process_sequences`：基于正向/反向 linker 的字符串匹配，提取两者之间的随机区域
- `process_sequences_regex_only`：基于正则表达式匹配，灵活度更高
- 输出按出现频次降序排列的 `fre.csv`（含 Sequence、Count 两列）

### 3. `selex/analysis.py` — 富集分析

函数：`compute_fold_change_by_round`

- 读取各轮次子目录（R1、R2、…）中的 `fre.csv`
- 以第一轮（R1）为基准，计算每条序列在后续各轮的 Fold Change
- 支持内连接（只保留每轮都出现的序列）或外连接（保留所有序列，缺失值补 0）
- 输出 `fold_change.csv`，列格式：`Sequence, Freq_R1, Freq_R2, ..., FC_R2, FC_R3, ...`

### 4. `selex/structure.py` — 结构分析

函数：`check_conditions` / `add_similarity_column_and_auc`

- `check_conditions`：调用 **NUPACK** 预测序列的 MFE 二级结构，计算与目标支架模板的相似度（`*` 位置不参与比对）
- `add_similarity_column_and_auc`：
  - 为 FC 表中每条序列拼接引物后预测结构、计算相似度
  - 生成 `Adj_FC`（FC × 相似度）列
  - 计算整体 AUC 及分段 AUC（从 R1 到 R2、R3、…）
  - 输出 `fold_change_with_similarity.csv`

## 运行方式

```bash
cd /home/hanlab/2022_zzm/3DSF-SELEX
python run.py
```

在 `run.py` 顶部修改配置参数：

```python
DATA_DIR = "/path/to/your/data"   # 数据根目录
PATTERN  = r"..."                  # 提取随机区域的正则表达式
ROUNDS   = [f"R{i}" for i in range(1, 11)]  # 轮次列表
```

结构分析默认关闭（计算耗时），取消注释最后一行即可启用：

```python
add_similarity_column_and_auc("fold_change.csv", "fold_change_with_similarity.csv")
```

## 数据目录约定

每个轮次子目录命名须以 `R` 开头（如 `R1`、`R2`），内部结构：

```
R1/
├── *1.fastq / *1.fq    # R1 reads
├── *2.fastq / *2.fq    # R2 reads
├── merged.fa           # 拼接结果（自动生成）
└── fre.csv             # 频次统计（自动生成）
```

## 依赖

```
pandaseq
cutadapt
nupack
pandas
numpy
scipy
```
