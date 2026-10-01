# Step 1 — Sequence QC

## Step 1 完成

### 1. 本步做了什么
从 AlphaFold PDB 的 Chain A 提取氨基酸序列，并把 PDB、UniProt Q925S4 与冻结的 181-aa 项目参考逐位对齐。UniProt 只用于映射和来源注释，未替换项目坐标。

### 2. 使用的算法
Biopython PDBParser 与 Biopython PairwiseAligner（Biopython 1.88）。

### 3. 算法属于什么
传统生物信息学的结构文件解析和全局双序列比对，不是新的深度学习训练。

### 4. 算法原理
通俗地说，它像把两篇文章逐字排齐，记录相同字符、替换、插入和删除。技术上使用全局 PairwiseAligner，match=2、mismatch=-1、gap open=-5、gap extend=-1，并用 PDBParser 读取标准氨基酸 residue 与 N/CA/C 主链原子。

### 5. 为什么本项目需要它
所有后续“第 N 位”分数都依赖同一坐标；此步骤防止 PDB、项目参考与 UniProt 编号被混用。

### 6. 输入
- `inputs/reference_181.fasta`
- `inputs/AF-Q925S4-F1-model_v6.pdb`
- `inputs/uniprotkb_accession_Q925S4_2026_08_20.fasta.gz`

### 7. 输出
- `results/01_sequence_qc/sequence_mapping.csv`
- `results/01_sequence_qc/pdb_sequence.fasta`
- `results/01_sequence_qc/uniprot_mapping.csv`
- `results/01_sequence_qc/sequence_qc.json`
- `reports/step01_sequence_qc.md`

### 8. 关键数值结果
- Project reference 长度：181 aa。
- PDB protein chains：A；Chain A 长度：181 aa。
- Reference↔PDB：matches=181，substitutions=0，insertions=0，deletions=0，reference match fraction=1.000000，PDB residue-number offset(s)=[0]。
- Reference↔UniProt：UniProt 长度=220 aa，matches=180，substitutions=1，insertions=39，deletions=0，offset(s)=[39]，substitution project position(s)=[31]。
- PDB 结构联动允许：True。

### 9. 当前结果怎样解释
本结果只证明序列对象和编号关系。它不证明任何 residue 是实验表位、受体界面或中和位点。

### 10. 是否存在问题
UniProt 与冻结项目参考的差异已完整保留在映射中；项目分析继续使用 1–181 坐标。若 `pdb_structure_linkage_allowed=false`，所有依赖该 PDB 的结构联动必须立即停止。

### 11. 下一步
在 PDB/reference 门槛通过后，Agent B 继续尝试官方 standalone BepiPred-3.0 `vt_pred` 推理；结构 Agent 只可在该门槛通过时继续。

## Runtime provenance
- Python 3.12.13
- Biopython 1.88
