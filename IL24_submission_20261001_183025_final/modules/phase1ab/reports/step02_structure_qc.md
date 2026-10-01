# Step 2：3D Structure QC

## 1. 本步做了什么

对 AlphaFold PDB 做了序列硬门控、链与残基完整性检查，并从每个 residue 的 B-factor 字段提取 pLDDT。所有结果使用 project 1–181 编号。

## 2. 使用的算法

Biopython PDBParser 1.88；AlphaFold per-residue pLDDT 读取。

## 3. 算法属于什么

传统结构生物信息学解析与质量控制，不是新的 AI 推理或训练。

## 4. 算法原理

PDBParser 把 PDB 原子记录组织为 model、chain、residue 和 atom 层级。AlphaFold 将局部结构置信度写入 PDB B-factor 字段；本分析取同一 residue 全部原子的均值。pLDDT 仅表示局部坐标可信度，不表示功能重要性。

## 5. 为什么本项目需要它

后续 RSA、DSSP 和结构表位分数必须严格映射到同一条 181-aa 序列；结构低置信位置也需要被显式标记以限制解释。

## 6. 输入

`inputs/AF-Q925S4-F1-model_v6.pdb`；`inputs/reference_181.fasta`。

## 7. 输出

`results/02_structure_qc/residue_plddt.csv`、`structure_qc.json`、`plddt_profile.png`。

## 8. 关键数值结果

- model 数：1；chain：A；残基数：181。
- 与冻结 reference 完全一致：True；结构分析允许继续：True。
- 缺失主链原子 residue 数：0；residue ID 断点：0；肽键几何断点：0。
- pLDDT minimum/mean/median/maximum：35.88 / 86.99 / 93.69 / 98.69。
- very high (>90)：119；high (70–90)：34；low (50–<70)：21；very low (<50)：7。

## 9. 当前结果怎样解释

pLDDT <70 的位置仅标记结构不确定性。高 pLDDT 不能证明 residue 具有抗原性或功能重要性，pLDDT 不进入最终生物学评分。

## 10. 是否存在问题

序列硬门控状态见 `structure_qc.json`。任何低 pLDDT 位置的结构来源指标应降级解释。

## 11. 下一步

只有序列门控通过后，才使用 FreeSASA 计算 SASA/RSA，并用 DSSP 添加二级结构注释。
