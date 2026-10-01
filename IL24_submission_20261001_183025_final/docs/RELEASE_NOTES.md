# RELEASE NOTES

2026-10-01局部更新：外置冻结官方ESM基座。新增prepare_base.py、固定对象版本/源SHA清单、独立A参考与张量指纹；复用原格式转换，不运行smoke训练。568张量和8条真实示例A/B均精确一致；保留全部3份项目checkpoint及所有原始结果，未重跑阶段一/二或完整训练。命令/磁盘/依赖见BASE_MODEL.md。本轮ZIP另命名为_external_base.zip，保留上一版大ZIP。

来源：服务器当前 `/mnt/ws/teeth` 三模块。Phase2 Git提交 `832af53d02341379637544cac3900f840f763e36` 加工作区未提交修复；Phase1/ESM当前目录无Git提交号，以文件摘要记录。`logs/source/source_files_sha256.csv`记录原来源，`release_code_changes.json`记录发布适配前后摘要。未使用旧compact包覆盖正式代码。

本版实际改动：
- 单一顶层入口，区分结果导出、Phase1实际融合排序、Phase2实际结构接触/受体分析、最终ESM加载推理；所有正常输出写新目录。
- 保留全部三模块、各自环境、真实数据划分与训练日志；交付原样LoRA/head、tokenizer与固定官方基座准备材料；大基座外置，运行按独立冻结张量身份校验。
- 修正版6DF3重新生成5代表重叠、48位受体足迹、T198映射、突变说明和SVG；图不读取旧发布表。修复置信度图固定范围漏显高分点、补足轴标签及图例，指标柱图统一使用0–1比例。更正ABodyBuilder2误差指标与当前摘要中的过期说明。
- 保留正式top10构象定义、原排序、种子/MSA/构建、所有主指标及序列；Full181补充分开；ESM不参与结构排序。未修改科学评分公式、训练标签、权重或阈值。
- 模型/数据来源、许可、原环境、完整训练/预测命令、附件5映射集中整理；论文PDF、旧ZIP、虚拟环境、下载缓存与重复审计文档不作为交付资产。

上一轮整包检查（2026-09-30）：导出、181残基融合和5家族、30结构接触（1899对）、5代表修正版6DF3、8真实序列最终ESM推理均已实际运行。新算Phase1共同数值字段与原表容差1e-12，候选身份/顺序不变；Phase2接触/摘要与原表相符；一个同条件代表的PRODIGY真实计算与原值完全相同。GPU/Boltz及训练/报告入口帮助、依赖导入检查通过。详见`logs/validation/prepackage_checks.json`。

上一版ZIP独立解压记录在`logs/validation/acceptance.json`；本次外置版A/B记录在`logs/base_externalization/equivalence.json`，新ZIP最终复核记录在同级`.validation.json`。文档或日志补入后最终ZIP另作一次解压复核，最终文件摘要与结果记录在ZIP同级sidecar；不把历史训练PASS当作本轮训练验证。

本轮未重新执行：完整训练、全量Boltz、上游BepiPred/DiscoTope/DeepSig、远程MSA、全新环境安装。仍需补充的记录：AbAgym原快照未列明确数据许可、BepiPred现有快照未发现LICENSE，提交方需补确认授权；部分第三方调用精确时间/预训练数据版本未记录；序列历史版本差异统一见Model Card。无最终ESM权重或正式入口缺失。代码与材料的技术交付项已具备，许可记录缺口不能视为已确认合规。

| 附件5要求 | 文件位置 |
|---|---|
| 源代码、入口、运行命令 | run.py/run.sh；modules；README；docs/FULL_RUN.md |
| 依赖与系统硬件 | envs；logs/validation/runtime_*.json；hardware.txt |
| 实际示例与最终模型 | data/examples；results；models/base_manifest.json、base_reference.json；tools/prepare_base.py；modules/inovate/output |
| 训练入口/配置/划分/日志 | modules/inovate/scripts/train.py；configs/esm2_train.json；模块data、logs和output |
| 数据与第三方溯源/许可 | docs/DATA_SOURCES.md、THIRD_PARTY.md、licenses；download_manifest.json |
| Model Card和适用限制 | docs/MODEL_CARD.md |
| 标准候选表与结构 | results/results.csv；results/structures |
| 复现与版本校验 | docs/REPRODUCIBILITY.md；logs/source；SHA256SUMS.txt；ZIP独立摘要 |

主入口覆盖真实关键流程，未另造重复Notebook。附件5是通用要求；本地未发现额外赛道字段/候选数量模板，因此不虚构额外门槛。
