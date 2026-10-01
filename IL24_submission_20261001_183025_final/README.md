# IL-24 计算项目提交包

本包交付三个实际模块：`phase1ab` 的免疫原区段设计与多证据热点排序，`phase2` 的抗体/复合物结构、接触、PRODIGY与受体界面分析，以及独立的 `inovate` ESM-2单序列代理分类与LoRA微调原型。所属赛道为**赛道一，AI大分子与多肽药物设计**。

阶段一区段选择支持Native155（项目27–181）和Immunogen134（27–160）的构建；热点作为阶段二的后验解释/代表构象选择先验，不是预测时的强制约束。Full181补充计算单独保留其两组15-seed结果。ESM-2不参与抗体排序。本项目未执行抗体序列分子优化；演示中的“设计/训练→优化→筛选”分别对应实际区段选择、多证据选择和结构筛选，以及ESM-2梯度训练/早停/分类，不能解读为已实现实验亲和力优化闭环。

## 最短运行：导出既有结果

在Linux解压后进入本目录（Python 3.10或以上；不需要GPU、联网或第三方Python库）：

```bash
bash run.sh --module export --output ./my_results
```

这会读取正式分析表、执行原有排序、生成UTF-8 `my_results/results.csv`并复制其实际引用的10个结构。**这是结果导出，不是结构重预测或ESM推理。** `my_results`必须是新的空目录。参考清单位于`results/results.csv`。结构路径相对CSV所在目录。

## 安装环境

模块可以使用独立环境。首选安装顺序见`envs/README.md`：CPU后处理使用`envs/postprocess.txt`；ESM先从PyTorch CUDA 12.4官方源安装torch 2.6.0，再安装`envs/esm2.txt`；完整结构预测使用`envs/phase2.txt`。导出仅需标准库。不要将完整GPU依赖当成导出的前置条件。

上一轮发布验证沿用旧环境；本次最终收尾的独立CPU/ESM环境安装和真实测试见ZIP同级验收记录。历史记录仍保留原日期和身份，不冒充本轮测试。

## 冻结基座准备（本版外置）

**官方基座冻结，本项目训练LoRA适配器与分类头。** 项目最终checkpoint仍完整随包，未合并进基座。为减小ZIP，2.6 GB官方基座的运行文件改为按固定清单获取；不重训、不量化、不换模型。

来源为Meta官方 `dl.fbaipublicfiles.com/fair-esm`：`esm2_t33_650M_UR50D.pt`及其contact-regression文件。权重发行使用S3对象版本，分别固定为`2G4typsUSwOKQLH35sqPMJ27ZLX1AEJG`、`cEcpCgy3_xlt3NCWsoXOzeIrh7.EDOBu`；每个文件另有完整SHA256。加载器fair-esm 2.0.0对应源码提交`0b59d87ebef95948c735b1f7aad463dc6dfa991b`。这不是猜测的历史Git权重revision，完整来源见`models/base_manifest.json`。

沿用原转换：Meta参数名映射为Transformers EsmModel，保留float32，rotary模式未使用的绝对位置表置零；原分类器不加载官方masked-LM head。tokenizer与模型配置保留历史版本。**官方文件SHA校验下载；转换后的模型校验由旧正式基座独立冻结的张量、配置、tokenizer指纹。两类SHA不互相比较。** 转换文件SHA仅作记录。

激活`envs/README.md`所列ESM环境后，在本目录执行：

```bash
# 自动：只获取必要的两个官方源文件；已校验缓存可复用
python tools/check_environment.py --module esm2
python tools/prepare_base.py
python run.py --module esm2 --output ./check_esm2

# 离线：将上述两个原始官方文件放入official_files目录
python tools/prepare_base.py --official-dir ./official_files --output ./external_models/esm2_650m_hf

# 离线：若已有转换目录，校验后直接指定它；不改写该目录
python tools/prepare_base.py --converted-dir ./offline_model --verify-only
python run.py --module esm2 --model-dir ./offline_model --output ./offline_esm2
```

下载约2.605 GB；转换输出约2.603 GB，模型准备建议额外6 GB磁盘、16 GB主机内存。转换只用CPU，ESM推理推荐GPU。依赖沿用Python3.10、torch2.6.0+cu124、fair-esm2.0.0、Transformers4.48.3、safetensors0.8.0、numpy1.26.4；安装顺序见`envs/README.md`。默认缓存为`~/.cache/il24-esm2-official`，运行基座为`external_models/esm2_650m_hf`，可用`--cache-dir`、`--output`、`--model-dir`显式指定。缓存与准备产物不属于ZIP静态文件清单，项目训练权重仍在`SHA256SUMS.txt`内。

上一轮外置基座独立A/B验证：568个张量的键、dtype、shape、值全部精确一致，8条真实样例的原始logit与概率最大差均为0，最终CSV与上一版相同。详见`docs/BASE_MODEL.md`和`logs/base_externalization/equivalence.json`。

## 三类真实计算复核

以下命令在对应环境激活后，从发布包根目录执行；每个输出目录须为空。

```bash
# CPU：复用随包官方模型证据，实际重做181残基融合、15/20/25-aa窗口和区域排序
python run.py --module phase1 --output ./check_phase1

# CPU：对30个正式MSA结构重新计算接触，并重新计算5个代表的修正版6DF3重叠和图
python run.py --module phase2 --output ./check_phase2

# GPU推荐：先按上节准备基座，再加载随包最终LoRA及分类头，对8条真实验证集序列推理
python run.py --module esm2 --output ./check_esm2_modules
```

Phase 1复用的BepiPred、DiscoTope、SASA/DSSP和保守性证据在模块`results/`下；本轮没有重跑上游模型或新增bootstrap。Phase 2所有接触来自随包PDB；6DF3是人IL-24受体复合物，项目参考为小鼠IL-24，映射与限定写在`docs/MODEL_CARD.md`。ESM示例来自AbBiBench CC-BY-4.0数据的既有验证划分，身份在`data/examples/esm2_manifest.json`；这8条不构成新增性能实验。

自定义输入与配置：

```bash
python run.py --module export --input modules/phase2/results/14_quality/model_ranking_formal_msa.csv --output ./custom_export
python run.py --module phase1 --input modules/phase1ab --output ./custom_phase1
python run.py --module phase2 --input modules/phase2 --output ./custom_phase2
python run.py --module esm2 --input data/examples/esm2_sequences.csv --config configs/esm2_inference.json --model-dir external_models/esm2_650m_hf --checkpoint modules/inovate/output/esm2_650m_lora_best.pt --output ./custom_esm2
```

`phase1/phase2 --input`接收与对应模块相同相对布局的数据目录；ESM接收CSV，必要列`model_sequence,label,source`（`label`用于现有评估入口，模型本身只读取序列；其他来源/行身份列予以保留）。长度1–1022，不静默截断。ESM输出`lora_validation_predictions.csv`、指标JSON和配置记录；分数是代理任务概率，阈值0.5。完整训练、预处理、上游模型、完整结构预测与评估命令集中在`docs/FULL_RUN.md`。

## 目录和材料

| 位置 | 内容 |
|---|---|
| `modules/phase1ab` | 原计算代码、区段选择分支、证据/中间表、冻结热点结果、原环境规格 |
| `modules/phase2` | 当前正式代码、全部既有Native/Immunogen结构/接触/PRODIGY，独立Full181补充材料 |
| `modules/inovate` | 预处理、固定划分、训练、评估、报告、最终权重及训练记录 |
| `models/base_manifest.json`、`base_reference.json`、`base_config` | 官方源固定版本、参考张量指纹、原配置与tokenizer；大基座不随ZIP |
| `models/esm2_assets.json` | 最终adapter/head与准确基座的SHA256、epoch/step |
| `models/third_party` | 原版BepiPred和DiscoTope小模型代码/资产；大权重获取见models/README |
| `results` | 参考清单、真实示例运行输出、引用的PDB |
| `docs` | Model Card、数据/第三方来源、完整运行、复现说明、发布说明 |
| `logs` | 历史验证记录与源码来源/未提交修复摘要 |

模块内部仍保持已有相对布局，以免引入无意义重构。历史日志中的旧路径只是溯源信息。`run.py`的正常运行输出不会覆盖参考结果。

## 候选与排序

顶层清单沿用正式入口的10条结构构象定义，按`confidence_score`再按`ipTM`降序；10是原提交入口已有选择，不是新设赛事数量门槛。各行共享`antibody_id=IA6-13-8`，`pose_id`区分构建和seed，不代表10个独立设计抗体。Full181没有被强行加进旧主清单；其30模型和独立排序保留在`modules/phase2/full181_supplement/results/full181_supplement/`。

字段：`candidate_id`为原清单行编号；`track`为赛道；`VH_sequence/VL_sequence/antigen_sequence`为实际序列；`construct/msa_condition/seed`定义条件；`source_pdb`是随包结构路径；`model_version/run_version`记录模型和本包版本；`PRODIGY_dG_kcal_mol`单位kcal/mol；confidence、ipTM、Boltz complex_pLDDT（0–1）与hotspot_coverage无量纲。预测能量、结构置信度和ESM代理概率不相互替代；无实测亲和力字段。PDB坐标为Å，A=抗原、H/L=抗体；未经松弛原结构的质子化状态未额外指定，可能有碰撞。

## 耗时与复核状态

导出通常秒级；CPU示例通常数十秒；ESM示例在本机A6000约十余秒（含冻结基座张量校验及加载）。本机两张RTX A6000 48GB、驱动570.133.07。完整LoRA历史运行4轮、24,996步并早停，最佳epoch3/step14000；完整结构预测按构建/条件各15 seeds，可能数小时以上，实际日志保留。上游下载、MSA服务、重新安装需要联网；ESM首次自动准备需要下载，离线准备后示例不需联网。

上一版三模块验收见`logs/validation/acceptance.json`；本次外置基座等价验证见`logs/base_externalization/equivalence.json`，最终新ZIP解压复核记录在ZIP同级`.validation.json`。完整训练、全量Boltz、BepiPred/DiscoTope推理本轮未重新执行，历史成功不替代本轮检查。

```bash
sha256sum -c SHA256SUMS.txt
```

SHA清单不包括自身，ZIP摘要在ZIP同级独立文件中。


## 最终ZIP真实验收与同级记录

本次发布标识为 `IL24_submission_20261001_183025_final`。科研结果中的 `run_version=IL24_submission_20260930_223945` 保留其正式结果版本含义；未修改候选、模型或排序。
与ZIP同时交付、位于其同级目录的文件为：

- `IL24_submission_20261001_183025_final.zip.sha256`
- `IL24_submission_20261001_183025_final.validation.json`
- `IL24_submission_20261001_183025_final.validation.md`
- `IL24_submission_20261001_183025_final.validation_logs.zip`

验收记录包含被实际解压测试的ZIP名称和SHA；它们在ZIP外，避免摘要循环依赖。不得用包内历史PASS代替sidecar的本轮状态。两项许可补核见 `docs/licenses/final_20261001/README.md`；其中唯一待团队补充的赛事授权材料，与技术测试状态分开判断。

在解压目录中，先按envs/README.md建好对应环境，再执行下列验收入口。该入口只调用正式入口并与随包参考比较，不训练、不重跑上游模型。各阶段失败或比较不一致均返回非零。`TEST_ROOT`必须是新的绝对路径；以下目录只是可直接执行的本地示例。

```bash
export TEST_ROOT="$PWD/acceptance_run"
mkdir -p "$TEST_ROOT"
. .env-post/bin/activate
python tools/accept_release.py --stage cpu --test-root "$TEST_ROOT"
deactivate
. .env-esm2/bin/activate
python tools/accept_release.py --stage base --test-root "$TEST_ROOT"
python tools/accept_release.py --stage esm --test-root "$TEST_ROOT"
python tools/accept_release.py --stage local --test-root "$TEST_ROOT"
```

本轮同机A6000、原bf16推理设置，在执行前固定分数比较容差为0（rtol=atol=0），逐行身份、全部CSV字段及标签比较；基座张量/配置/tokenizer指纹要求精确一致。T1–T3比较完整结果内容而非只看计数。CPU环境只需要postprocess.txt；GPU环境需要上述torch及esm2.txt。准备完成后的本地运行不再下载基座，不等于无预置资源也能离线恢复。
