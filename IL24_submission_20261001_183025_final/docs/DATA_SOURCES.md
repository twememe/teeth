# 数据来源、用途与划分

历史未记载字段标为“未记录”。本轮未重新清洗标签、重切划分或补造访问日期。文件级原URL、版本、长度和SHA256见`modules/inovate/data/download_manifest.json`；项目来源文件摘要见`logs/source/source_files_sha256.csv`。

| 数据/资产 | 固定版本、来源 | 获取时间记录 | 用途与许可 |
|---|---|---|---|
| 项目IL-24 181-aa、VH/VL与ELISA表 | 项目提供；FASTA与现有PDB保持同版 | 具体提交时间未记录；外部归档2026-09-06 | 区段、结构、解释；本项目授权材料，无新增通用再许可 |
| AlphaFold结构 | AF-Q925S4-F1-model_v6；https://alphafold.ebi.ac.uk/entry/Q925S4 | 2026-08-21输入记录 | 181-aa结构证据；沿用AlphaFold DB条款与CC-BY-4.0来源注释 |
| UniProt/哺乳动物同源序列 | Q925S4、Q13007及`ortholog_metadata.csv` accession/URL | 2026-08-21表内记录 | 比对、同源保守性；UniProt CC-BY-4.0；选择25个代表、去除明显片段/重复，MAFFT比对 |
| 6DF3 | RCSB https://www.rcsb.org/structure/6DF3；PDB/CIF原文件 | 2026-09-06归档 | 外部受体界面；PDB坐标按wwPDB开放数据条款（CC0）使用；论文PDF不放入发布包 |
| AbBiBench | HF AbBibench/Antibody_Binding_Benchmark_Dataset，556fd6913aa231c0d342a8658818be8f963fd582 | 2026-09-22运行记录；逐文件时间未记录 | 215,699条原记录；训练/验证/测试；原README声明CC-BY-4.0 |
| AbAgym | https://github.com/3BioCompBio/AbAgym ，87af82ad68bb90921e15fb8b79c7bcb1d05b23d3 | 2026-09-22运行记录；逐文件时间未记录 | 36,541条interface记录及PDB链；同一固定快照README明确限非商业使用；赛事再分发范围待团队确认，证据见licenses/final_20261001/README.md；不含需另行FoldX许可的衍生结构 |
| IL6 | HF alchemab/il6-binding-prediction，f83d7fa0d08c143aefb4bb5b8a4e55ed0576c34f | 2026-09-22运行记录；逐文件时间未记录 | 1,636条；原README声明CC-BY-NC 4.0，来源AVIDa-hIL6；保留原划分 |
| ESM预训练语料 | 官方UR50D预训练权重；语料快照ID及项目下载语料时间未记录 | 项目未重新收集预训练语料 | 仅使用预训练权重，不宣称重新预训练 |
| BepiPred/DiscoTope/ABodyBuilder2/Boltz训练数据 | 由各第三方模型训练；本项目未重新获取其训练集 | 未记录 | 仅引用官方模型/方法资料，不伪造项目训练集来源或隐藏集检查 |

## 预处理与划分

`modules/inovate/scripts/preprocess.py`是实际处理入口。AbBiBench以普通H+L拼接表示，训练部分中位数7.828436341974用于二值化；等于阈值的记录排除。AbAgym验证残基/链映射后构建实际突变链，WT/site/mutation仅保留为元数据；其中多数为抗原突变，高DMS数值可能意味着免疫逃逸。IL6为单VH序列。无效序列、长度>1022和无效标签剔除，不静默截断；拒绝记录表保留。

精确序列整体分组，seed42。原始253,876条；实际训练199,999，验证20,000，测试33,508。集合间完全相同序列重叠为0；近同源、家族和抗原未严格隔离。隐藏赛事评测集未访问，不能宣称已排除未知隐藏集重叠。完整划分CSV保留原样，SHA256可核对。

小规模ESM示例仅取AbBiBench许可明确的既有验证记录，8个原始行身份在`data/examples/esm2_manifest.json`，没有随机合成数据。阶段一复用真实逐残基表；阶段二示例结构就是主分析的正式MSA30模型及5个代表。

许可补核（2026-10-01）：固定快照的非商业使用条款已找回；原始及加工数据向组委交付的授权范围仍待团队确认。具体证据、署名和唯一待补材料见 `licenses/final_20261001/README.md`。保留原数据、划分及训练成果，不以删除原始文件代替处理加工数据许可。
