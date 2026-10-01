# Step 4：Structural Surface Accessibility

## 1. 本步做了什么

在通过 181-aa 序列硬门控的 AlphaFold 坐标上，逐 residue 计算 absolute SASA 和 RSA，并用 DSSP 添加二级结构注释。所有表严格使用 project 1–181 编号。

## 2. 使用的算法

FreeSASA 2.2.1（LeeRichards，probe radius 1.40 Å）；mkdssp version 4.6.1，通过 Biopython DSSP 1.88 调用。

## 3. 算法属于什么

FreeSASA 是几何表面积计算；DSSP 是传统结构生物信息学注释。两者都不是 AI 模型。

## 4. 算法原理

FreeSASA 模拟一个水分子尺度的探针沿蛋白表面滚动，计算探针可接触的原子表面积；absolute SASA 是面积，RSA 是按 residue 类型最大可接触面积归一化后的相对值。DSSP 根据主链氢键与几何模式判定 helix、strand 和 coil。

## 5. 为什么本项目需要它

抗体必须物理接触抗原表面，因此 RSA 提供独立的结构可接近性证据。DSSP 只帮助解释局部结构背景，不自动加分或减分。

## 6. 输入

`inputs/AF-Q925S4-F1-model_v6.pdb`；`inputs/reference_181.fasta`，且结构序列必须完全一致后才运行。

## 7. 输出

`results/04_surface/residue_sasa.csv`、`residue_secondary_structure.csv`、`rsa_profile.png`。

## 8. 关键数值结果

- FreeSASA total SASA：11579.18 Å²。
- residue absolute SASA minimum/mean/median/maximum：0.00 / 63.97 / 55.02 / 246.68 Å²。
- RSA minimum/mean/median/maximum：0.0000 / 0.3854 / 0.3752 / 1.2775。
- DSSP helix/strand/coil：115 / 12 / 54 residues。

## 9. 当前结果怎样解释

较高 RSA 表示 AlphaFold 模型中更暴露、更容易被溶剂探针接触；它支持表面可接近性，但不能单独证明抗体结合。DSSP 仅为 annotation，不进入最终 score。

## 10. 是否存在问题

SASA/RSA 依赖当前单体 AlphaFold 构象，并受低 pLDDT 区域坐标不确定性影响；低 pLDDT residue 的结构指标必须降级解释。

## 11. 下一步

Integration Agent 仅在严格 181 行映射下使用 RSA；DSSP 保留为解释性注释。
