# 外部资料归档：IL-24 / 6DF3

归档日期：2026-09-06

本目录只收录与 Phase 2 受体竞争分析直接相关的外部资料：两篇论文、ELISA 原始 CSV，以及从 RCSB PDB 6DF3 下载的结构和序列文件。

## 文件说明

- `Lubkowski_2018_6DF3_crystal_structure.pdf`：6DF3 三元晶体结构论文。论文报告 IL-24–IL-22R1–IL-20R2 复合物，分辨率 2.15 Å。
- `Flexible_regions_IL24_receptors.pdf`：IL-24 与 IL-20R1/IL-22R1 受体结合柔性区域研究论文，包含 T198 相关证据。
- `ELISA_competition_data.csv`：用户提供的 ELISA 抑制率原始数据，未改写数值。
- `6DF3.pdb` / `6DF3.cif`：RCSB 结构坐标（下载地址：`https://www.rcsb.org/structure/6DF3`）。
- `6DF3.fasta`：RCSB 提供的结构实体序列。

## 结构核查

6DF3 的 PDB 坐标包含三条主要链：

- C：IL-24，结构片段对应 UniProt 52–206；
- L：IL-22R1 胞外结构域，结构片段对应 UniProt 24–228；
- H：IL-20R2 胞外结构域，结构片段对应 UniProt 35–224。

因此，6DF3 可以直接用于 IL-24–IL-22R1–IL-20R2 的受体界面分析。它不能直接代表 IL-20R1–IL-20R2 复合物。

## 证据边界

- 论文中的受体接触残基是结构/文献证据，可用于定义分析参考界面；不能单独证明 IA6-13-8 的竞争机制。
- T198 属于文献编号体系。必须先与本项目 IL-24 序列程序化比对，再转换为项目编号；在映射完成前不使用项目残基号 198。
- ELISA 数据可作为外部一致性证据；在包被物、加入顺序、浓度和抑制率公式确认前，不将其命名为受体竞争 IC50 或真实 KD。
- 聊天截图仅包含背景说明和建议，不作为实验原始数据归档。

## 来源

- RCSB PDB 6DF3：`https://www.rcsb.org/structure/6DF3`
- 结构文件：`https://files.rcsb.org/download/6DF3.pdb`、`https://files.rcsb.org/download/6DF3.cif`
- 序列文件：`https://www.rcsb.org/fasta/entry/6DF3`
