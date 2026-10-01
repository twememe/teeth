# Step 6 — Evolutionary Conservation

## Step 6 完成

### 1. 本步做了什么
从 UniProt 官方 REST API 的 mammalian IL-24 查询结果中，按预注册质量规则选择 one species/one representative 的 25 条序列，与冻结的 181-aa project reference 一起用 MAFFT 比对，并逐位计算保守性、gap fraction 与 100 次固定种子 bootstrap 区间。

### 2. 使用的算法
UniProt REST API；MAFFT `v7.526 (2024/Apr/26)`；normalized Shannon entropy；100-replicate sequence bootstrap（seed=2406）。

### 3. 算法属于什么
UniProt 检索和 MAFFT 属于传统生物信息学；Shannon entropy 和 bootstrap 属于统计方法。没有训练或微调 AI 模型。

### 4. 算法原理
通俗地说，MAFFT 把不同物种的 IL-24 像多篇文章逐字排齐；某列越多物种使用同一种氨基酸，该位置越保守。技术上，去除 gap 后计算 `H=-sum(p*ln p)`，再转换为 `1-H/ln(20)`；gap 单独报告，不伪装成氨基酸。每次 bootstrap 对物种序列有放回抽样，并重新计算 181 个位置。

### 5. 为什么本项目需要它
它提供独立于序列表位、三维表位和表面可接近性的 evolutionary evidence，用于判断哪些 project residues 在哺乳动物中长期稳定。

### 6. 输入
- `inputs/reference_181.fasta`
- `results/06_conservation/raw/uniprot_il24_mammalia.json`

### 7. 输出
- `results/06_conservation/orthologs.fasta`
- `results/06_conservation/ortholog_metadata.csv`
- `results/06_conservation/ortholog_selection_audit.csv`
- `results/06_conservation/mafft_alignment.fasta`
- `results/06_conservation/residue_conservation.csv`
- `results/06_conservation/conservation_profile.png`
- `reports/step06_conservation.md`

### 8. 关键数值结果
- UniProt query records：116；selected species：25；reviewed selected：3。
- Mammalian orders (10)：Artiodactyla, Carnivora, Chiroptera, Dasyuromorphia, Diprotodontia, Primates, Proboscidea, Rodentia, Sirenia, Tubulidentata。
- Conservation min/mean/median/max：0.381188 / 0.798114 / 0.823201 / 1.000000。
- Gap fraction min/mean/max：0.000000 / 0.030055 / 0.120000。
- Exclusion/audit reasons：{'same_taxid_lower_priority_than_Q9JI24': 2, 'eligible_not_selected_target_25': 38, 'explicit_fragment': 5, 'length_below_140': 27, 'same_taxid_lower_priority_than_A0A2K5P2T7': 1, 'same_taxid_lower_priority_than_A0A7J8CR82': 2, 'noncanonical_amino_acid': 2, 'same_taxid_lower_priority_than_A0A2K5IP10': 1, 'same_taxid_lower_priority_than_A0A2K5QSY6': 1, 'same_taxid_lower_priority_than_A0A2K6PT46': 1, 'same_taxid_lower_priority_than_A0A8C9UJS4': 3, 'same_taxid_lower_priority_than_H2Q108': 1, 'same_taxid_lower_priority_than_I3M7Y2': 1, 'same_taxid_lower_priority_than_G1QJB7': 1, 'same_taxid_lower_priority_than_A0A2K6KH42': 1, 'same_taxid_lower_priority_than_A0A671FWG8': 1, 'same_taxid_lower_priority_than_F6YRZ5': 1, 'length_above_280': 1, 'same_taxid_lower_priority_than_A0A2K5UB53': 1, 'isoform_name': 6, 'explicit_low_quality_name': 2}。

### 9. 当前结果怎样解释
高 conservation 表示该对齐位置在所选哺乳动物 IL-24 中氨基酸组成较一致；它不证明实验表位、受体界面或中和功能。高 gap fraction 位置需要降低解释置信度。

### 10. 是否存在问题
只有三条查询结果是 reviewed Swiss-Prot；其余代表来自透明筛选的 unreviewed entries。Bootstrap 反映当前 ortholog set 的抽样不确定性，不覆盖数据库注释错误或系统发育非独立性。

### 11. 下一步
Integration Agent 只读取冻结的 181-row conservation table，并在不访问历史区间的前提下与其他独立证据整合。
