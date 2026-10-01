# ESM-2 单序列代理分类：LoRA 可行性验证报告

本报告评估的是单序列二分类代理任务：模型输入只有一条 `model_sequence`，没有抗原条件化输入，也没有同一抗原下的野生型/突变体成对比较，因此不含亲和力或突变效应预测。

真实 ESM-2 650M 已完成前向、反向、LoRA 更新、双卡同步以及保存恢复验证。LoRA 在本验证集上的 MCC 高于冻结基线。
本实验评估的是混合数据上的二分类代理任务；不能据此宣称对未见抗体家族、抗原或真实亲和力改善具有泛化能力。

## 同一验证集的实测对比

主要模型选择指标为 MCC，分类阈值固定为 0.5。下列三组使用完全相同的 20,000 条验证记录。

| 模型 | MCC | F1 | Spearman ρ |
|---|---:|---:|---:|
| 随机预测（实测） | -0.003354 | 0.480650 | -0.006402 |
| 冻结 ESM-2 650M + 线性头 | 0.385864 | 0.710010 | 0.485003 |
| ESM-2 650M + LoRA + 同构线性头 | 0.499083 | 0.754309 | 0.597148 |

LoRA − 冻结基线的 MCC 差值：+0.113219；按相同序列分组的配对 bootstrap 区间为 [+0.103726, +0.124905]（500 次，种子 42）。同一验证集同时参与了模型选择，该区间不能消除选择偏差。
随机基线的 MCC=0、F1=0.5 是类别均衡条件下的理论参考值，不是本次实测值。表内 Spearman 使用预测概率与二分类标签；各来源原始连续分数的相关性另列在 metrics JSON 中，未将不同量纲的分数直接混合。

## 训练与预算

- 开始：2026-09-22T00:37:56.557353+08:00；总截止：2026-09-22T13:37:56.557353+08:00。
- LoRA 已开始 4 轮，完整完成 4 轮；优化步骤 24996。15 轮是上限，不代表全部完成。
- 提前停止：True；达到训练时间上限：False。
- 冻结基线训练 4 轮，学习率 0.001，缓存特征批量 4096。LoRA 从该基线线性头继续训练，分类头结构相同。
- Ab-Tune 固定依赖与现有 CUDA PyTorch 版本冲突，按授权改用 PEFT + Transformers；层 0–15 的 query/value 注入 LoRA，r=8、alpha=16、dropout=0.1。
- 使用 bf16 混合精度、长度分桶、梯度检查点。具体实际批量和设备数见 lora_run.json 与 training_config.yaml。

最佳权重来自第 3 轮的第 14,000 步。早停按每轮末 MCC 的连续未提升轮数计算（patience=3）；每 2,000 步的中途验证也用于保存最佳权重。
最终冻结基线与 LoRA 均重新加载保存权重，使用相同 bf16 推理流程重评完整验证集；逐条行标识、标签和各分组指标已独立复算。
损失图浅色线为每 20 步记录的小批量损失，深色线为连续 50 个记录点的滑动平均；第 10,000 步标示恢复位置。第三面板显示验证集 MCC 与模型选择。

## 中断与一次修复

正式训练在第 10,000 步验证后，因 ESM 旋转位置编码缓存与 inference_mode 的反向传播限制中断。已复现并将验证改为 no_grad；从保留权重继续，已消费批次不重放，原截止时间不变。
原检查点没有优化器和随机数状态，因此恢复了 LoRA/分类头权重，AdamW 与随机数状态重新初始化；这不是优化器状态完全连续的训练。最佳模型仍按全部已完成验证的 MCC 保留。完整错误、修复测试及原始检查点均在 logs/、error_log.txt 和 output/recovery_20260922_0520。

## 数据审计

- 原始数据：AbBiBench 215,699、AbAgym 36,541、IL6 1,636，合计 253,876。
- 训练/验证/测试：199,999 / 20,000 / 33,508。阈值相等记录明确排除；未复制记录补足条数。
- AbBiBench 二值化中位数仅由训练部分估计：7.828436341974。
- 共 62,233 个相同输入序列组在不同实验中标签冲突；相同序列整体分组，跨集合重叠为零。
- AbAgym 的 WT/突变字段实际是残基和 PDB 位点，已校验野生型残基并重建突变链。WT:mutation 格式保留为元数据，模型接收真实氨基酸序列。
- AbAgym 主要记录为抗原突变，且高 DMS 分数可表示免疫逃逸；逐来源及突变对象结果见下表。
- AbBiBench 和 IL6 原始数据没有可靠突变注释，mutation 列保留为空并标注来源，不编造突变信息。
- IL6 原始划分保留；测试集在本次模型选择与主报告中未使用。近同源序列、抗体家族和抗原之间尚未做严格隔离。

## 按来源与突变对象的 MCC

| 分组 | 冻结基线 | LoRA |
|---|---:|---:|
| abagym/antibody (n=194) | 0.000000 | 0.319150 |
| abagym/antigen (n=2823) | 0.000000 | 0.402763 |
| abbench/antibody (n=16820) | 0.406207 | 0.492763 |
| il6/antibody (n=163) | 0.000000 | 0.220834 |

## 交付文件

- 最佳 LoRA 与分类头权重：`output/esm2_650m_lora_best.pt`（预训练主干位于 `../../models/esm2_650m_hf`）。
- 完整配置：`training_config.yaml`；预处理：`scripts/preprocess.py`；下载：`scripts/download_data.py`。
- 对比表：`output/comparison.csv`；详细指标：`output/metrics.json`；逐条预测：`output/*_validation_predictions.csv`。
- 损失曲线：`output/loss_curve.png` / `.pdf`；原始曲线数据：`output/losses.csv`。
- 训练日志：`logs/` 和 `output/run_rank*.log`；进度：`progress.txt`；错误与修复记录：`error_log.txt`。
- 数据版本、哈希：`data/download_manifest.json`；划分审计：`output/data_audit.json`。

全部任务文件位于本目录（项目根由 `scripts/` 的位置解析，不依赖绝对路径）。

## 原始来源

- [AbBiBench](https://huggingface.co/datasets/AbBibench/Antibody_Binding_Benchmark_Dataset)
- [AbAgym](https://github.com/3BioCompBio/AbAgym)
- [IL6](https://huggingface.co/datasets/alchemab/il6-binding-prediction)
- [Ab-Tune 包信息](https://pypi.org/project/Ab-Tune/)
