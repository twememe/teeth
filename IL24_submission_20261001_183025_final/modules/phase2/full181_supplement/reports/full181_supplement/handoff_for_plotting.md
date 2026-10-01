# Full181 补跑绘图数据交接

日期：2026-09-21。Full181正式预测、原始产物、接触/CDR、几何QC、ensemble、PRODIGY与共同区对照全部完成；Full30/30，旧60个正式模型只读复用并核验。本轮未生成图、图注或图结论。

本地数据根为 `results/full181_supplement/`。服务器统一归档根为 `/mnt/pxh/teeth/`，项目位于 `phase1ab/` 与 `phase2/`，Full补跑包位于 `phase2/full181_supplement/`。实际执行时使用 `/mnt/il24_full181_supplement_20260921`，原始run_state保留该路径，不为归档改写历史。

关联键固定为 `construct + condition + seed`，不能仅使用seed。condition为empty/formal_msa；聚合行seed=ALL，跨seed配对行另有seed_left/seed_right，不视作生物学配对。原执行、本地和teeth路径对应见artifact_path_map.csv；绘图软件应优先读取selection_manifest_local.csv中的本地代表路径。远端与本地均完成原始结构/序列/置信度/PAE核验；迁移不通过猜测权重路径来绕过验收。

## 坐标、阈值、单位和缺失

- Full181的local=project，范围1–181；Native155范围27–181，Immunogen134范围27–160，旧构建仅在序列和PDB残基身份核对后使用local+26。显式映射在sequence_audit/numbering_map_all_constructs.csv，Full专表为numbering_map_full181.csv。重编号或插入码必须按映射处理。
- 主界面为A对H+L蛋白重原子最小距离严格<4.5Å，扩展阈值严格<5Å。水、氢、非蛋白HETATM排除；H–L内部接触不是抗原界面。
- JSON confidence/pTM/ipTM/complex pLDDT均0–1；PDB B-factor pLDDT为0–100；PAE为Å，Full矩阵406×406。链对ipTM由实际processed records确定A/H/L身份；不要以VH–VL分数替代抗原结合置信度。
- ΔG为kcal/mol，Kd为M，温度为°C。Full30条PRODIGY实际stdout温度均25°C，选择A对H,L。旧60条数值来自原保存表，其实际stdout温度未包含于交付材料，保留NA。
- 支持率、界面区段比例及Jaccard均0–1。构建外位置为NA/不适用；构建内已完成模型无接触才为0。无表位模型、缺失模型和评分失败是不同状态。本次30个Full均有接触及有效PRODIGY评分。
- 每条件独立15个seed，共识≥8/15，不能合并为30-seed共识。Jaccard表保留both_empty等状态，两个空集合的旧式Jaccard=1不能解释为稳定结合。
- 本轮90个结构都有旧阈值下的原子冲突警告，保留全部模型。CA断链、主链缺失和肽键C–N异常均为0，不作事后选择性筛除。

## 核心文件

每张CSV实际字段与行数以table_schemas.json为准；以下列出使用入口。

| 文件/目录 | 内容与粒度 |
|---|---|
| input_audit.json、sequence_audit/ | 181 aa、H180/L181、抗体同一性、构建切片和9项编号测试证据 |
| inputs/、msa/ | 两种正式YAML，同批A/H/L CSV，原A3M、API响应、query、manifest、检索时间与配对key来源 |
| processed_msa_audit.csv | 90条模型×链记录；实际MSA缓存A1251/H8192/L8192与本批CSV解析逐项匹配；内部模型张量深度明确未知 |
| task_manifest.csv、raw/ | 30个唯一任务；每seed原PDB、confidence、PAE、processed输入及原run_state |
| execution_per_seed_audit.csv | 每seed验收、开始/结束、耗时、GPU采样峰值、源路径 |
| full181_model_summary.csv | 30行置信度、表位/paratope/CDR计数、链对ipTM、pLDDT量纲、几何状态与来源 |
| full181_interface_contacts.csv | 2232条<5Å残基对；local/project/PDB编号、插入码、抗体IMGT、CDR/FR、距离和阈值标志 |
| full181_cdr_participation.csv | 420条模型×H/L×CDR/FR区域记录，参与标志和实际局部残基集合 |
| full181_geometry_qc.csv | 30行，既定阈值、断链/肽键/原子冲突/缺失主链；警告不剔除 |
| full181_epitope_support.csv | 362条条件×project残基记录；支持次数、计划/完成分母、支持率及共识状态 |
| full181_epitope_jaccard.csv | 每条件105个seed对，共210行；空集合状态另列 |
| full181_ensemble_summary.csv | 两条件分别15个seed的均值/SD/分布、共识和QC观察数 |
| full181_model_ranking.csv | 30个模型的confidence、consensus、Phase1-informed排序与依据 |
| full181_prodigy.csv | 30条有效评分，ΔG/Kd/温度/接触数、退出码、原始stdout/stderr与命令路径 |
| three_construct_model_comparison.csv | 90个不同构建×条件×seed；Full当前原始分析与旧保存表来源分开；旧原始包核验状态明确 |
| three_construct_comparison.csv | 同MSA条件下三构建的六个ensemble，均15个模型 |
| three_construct_metric_distribution.csv | metric/value长表，附证据来源；用于描述性分布，不能将计算seed视为独立生物重复 |
| three_construct_epitope_support.csv | 1086条构建×条件×project残基；构建外为NA |
| three_construct_cdr_participation.csv | 三构建的H/L CDR与FR参与记录 |
| shared_region_model_epitopes.csv | 先裁剪到project27–160、27–181后的逐模型表位集合 |
| shared_region_epitope_comparison.csv | 1808行；同条件的跨构建225个seed组合及共识比较；三构建共同区27–160，Full/Native扩展共同区27–181 |
| terminal_contact_summary.csv | 90行Full模型×区段，1–26/27–160/161–181接触数及占完整epitope比例 |
| representatives/、selection_manifest_local.csv | 六个代表称号、五个独立PDB，可直接读取本地绝对路径；保留原执行与teeth归档路径 |
| legacy_raw_audit/ | 实际核验旧60套原始文件；60行summary/QC、840行CDR、4459条接触、逐模型validation及舍入差异记录 |
| server_baseline_snapshot/ | 从原服务器只读取得的旧正式输出、表格、prior、原脚本和日志；source_inventory指向原服务器文件 |
| baseline_audit/ | 两条件原表来源、60套原始包核验接口、44项旧汇总回归核对、冻结prior原件来源 |
| input_fingerprints.csv、artifact_fingerprints.csv | 路径、大小、mtime；未计算哈希，不能称为密码学完整性证明 |
| artifact_path_map.csv、local_transfer_validation.json | 30套本地原始产物验收与90项PDB/confidence/PAE路径映射，原run_state未修改 |
| delivery_verification.json | 实际模型/分析/评分/本地产物数量、旧文件元数据核对和远端/本地验收 |

## 代表选择规则

confidence按confidence_score降序；consensus按同条件其他模型的平均表位Jaccard降序，再按confidence；Phase1-informed按原冻结62-residue prior覆盖率，再按confidence。最终并列均按数字seed升序。prior来自原服务器results/12_complex_prediction/phase1_prior_deduplicated.csv，原件保存在inputs/，未按Full结果改变。

empty：confidence seed2，consensus seed6，informed seed11。formal MSA：confidence seed6，consensus和informed均seed9。重复称号共用同一结构，只计一个模型。interim_analysis_20260921/是运行中间快照，不作为最终代表入口。

## 复现和使用边界

原完整GPU命令在逐seed command.json与run_state中，环境和权重路径在environment_fingerprint.json。已完成的30个模型不应因换目录或排序偏好再次运行。严格运行时验收使用原执行根；归档/本地副本通过显式path map读取，不能把POSIX路径缺失误判为预测失败或擅自改写原权重指纹。

现有原始数据的CPU重算入口包括analyze_full181_supplement.py、run_full181_supplement.py的PRODIGY导出模式、full181_ensemble.py、audit_full181_execution.py和audit_full181_legacy_raw.py。原环境和日志在同一服务器保留。PRODIGY已有同产物同参数的评分直接复用，避免覆盖第一次原始日志。

Full181是保留项目N端的计算构建；结构接触、共识、结构置信度、计算亲和力代理分别解读，不称为真实表位、实验KD、中和或受体竞争已经确认。新旧MSA检索批次不同，不声称除抗原长度外所有输入完全一致。

SignalP/DeepSig、阶段一重算、NAG、松弛、6DF3修复及其他模型/实验本轮均未执行。phase1ab仅复制既有项目文件。此前BLOCKED报告包和pre_root记录属于历史状态，最终结论以本报告、Full181_supplement_report.md和delivery_verification.json为准。
