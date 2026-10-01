# IL-24 Full181：阶段二补跑执行提示词

你现在接手已有 IL-24 / IA6-13-8 抗体项目。请在原工程中实际完成 Full181 的阶段二补跑，而不是只提供新计划。顺序是：审计现有工程与输入 → 做必要的增量代码修改和测试 → 补跑 → 下游分析 → 与旧构建对照 → 输出验收报告。本轮不画图，也不重做阶段一。

## 一、已确认事实与本轮范围

1. 用户刚上传的 FASTA 实际为181 aa，不是口头所说的180 aa。四行序列长度为60、60、60、1；第180位为H，第181位为L。禁止删除末端L来凑180，也不能自动裁掉N端。
2. 本次本地检查确认：该序列与旧包 data/IL24_full_1_181.fasta、参考PDB的Chain A完全一致；Native与Immunogen也分别等于它的27–181、27–160截取。在你的执行环境仍须复核，不能只引用这段说明。
3. SignalP和DeepSig后来均未真正完成。本轮不安装、不补跑、不捏造它们的结果。历史1–26边界只按项目既有工程定义/注释处理，不写成这两个软件的实测输出。
4. 阶段一已在181-aa参考上完成抗原侧分析，本轮不重跑BepiPred、DiscoTope、RSA、DSSP、保守性或免疫原筛选。阶段二继续使用同一个IA6-13-8的VH119/VL106，不设计新抗体、不改序列、不训练新模型。
5. 本轮必要闭环：Full181复合物预测 → 接触表位/paratope → 置信度与ensemble → PRODIGY → 三构建对照。本轮新增正式模型30个：Full181 × {empty MSA, formal ColabFold MSA} × seeds 1–15。旧Native/Immunogen只读复用，不重跑旧GPU预测。
6. 本轮不顺带启动N74单NAG重预测、OpenMM松弛、6DF3修复、SignalP/DeepSig、FoldX/Rosetta或其他新模型。可以复用既有几何QC代码作低成本检查。上述扩展保留为未补齐项，不能把主线完工写成所有阶段二扩展均已完成。

## 二、唯一序列与坐标契约

以下是已核验的Full181序列；FASTA换行不计入长度：

>IL24_Full181_project_1_181
MSWGLQILPCLSLILLLWNQVPGLEGQEFRSGSCQVTGVVLPELWEAFWTVKNTVQTQDD
ITSIRLLKPQVLRNVSGAESCYLAHSLLKFYLNTVFKNYHSKIAKFKVLRSFSTLANNFI
VIMSQLQPSKDNSMLPISESAHQRFLLFRRAFKQLDTEVALVKAFGEVDILLTWMQKFYH
L

规范化序列SHA256（只对拼接后的大写氨基酸字符串计算，不含header、空白和换行）：
3642351ae1bb5ca0788598425de59f455672266ed7f434d21ae1bfe654992374

固定三种构建：
- Full181：project 1–181，length=181，offset=0。
- Native155：project 27–181，length=155，offset=26。
- Immunogen134：project 27–160，length=134，offset=26。

若模型抗原链使用从1开始的连续局部编号，则project_position=local_position+offset。先核对输出链序列与编号；遇到缺失/插入/重编号时建立显式映射，不能机械套公式。

必须加入单元测试：Full local1→project1=M，local27→project27=Q，local74→project74=N，local160→project160=A，local181→project181=L；Native/Immunogen local1→project27=Q，local48→project74=N，local134→project160=A；Native local155→project181=L。验证project74–76为NVS。旧糖基化输入的local48不得移用于Full，尽管本轮不生成NAG模型。

FASTA标题Q925S4只作为原文件元数据，不据此自动下载/替换成数据库另一版本的220-aa序列。发现当前工程序列与上述SHA256不符，先输出逐位差异并停止受影响的GPU任务，禁止静默修正。

## 三、先做工程审计，保护旧结果

先阅读工程中的AGENTS.md（若有）、任务记录、阶段总结、MODEL_CARD.md、REPRODUCIBILITY.md，并定位真正的项目根目录。不要假定当前路径是完整服务器目录，也不要依赖ChatGPT沙箱绝对路径。

重点检查实际存在的：
- data/IL24_full_1_181.fasta、IL24_native_27_181.fasta、IL24_immunogen_27_160.fasta、VH_VL_clean.fasta；
- data/antibody_imgt_numbering.csv、cdr_annotation.csv、il24_numbering_map.csv；
- src/generate_boltz_inputs.py、generate_colabfold_msa.py；
- src/run_boltz_formal.sh、run_boltz_formal_msa.sh；
- src/analyze_boltz_complexes.py、run_prodigy_formal.py、rank_formal_models.py、rank_formal_msa_models.py、summarize_formal_msa.py；
- 原始PDB、confidence/PAE、旧MSA、逐模型接触表、PRODIGY表、冻结Phase1 prior；
- 模型checkpoint、实际虚拟环境、GPU和磁盘资源。

已知当前精简release并不包含完整60套原始模型；部分表格和代表结构存在，logs为空，权重没有打包。区分“历史报告记录完成”“当前找到原始文件”“只有汇总表”三种状态。不能用旧results.csv的10个候选代替完整ensemble，也不能把results_immunogen.csv为空解释为该构建没运行。

建立只读输入/旧结果SHA256清单。新增目录建议为results/full181_supplement/、logs/full181_supplement/及reports/full181_supplement/，或在原结构中建立同等隔离分支。所有新增表都携带construct、condition、seed、source_path；不要覆盖旧results.csv或旧统计表，不随意改变工程地址。

优先新增参数化入口或小型补跑脚本；如修改共用代码，保留diff并验证旧Native/Immunogen映射及指标不意外变化。不要直接运行会覆盖旧输入的prepare_phase2_inputs.py，也不要把仅导出旧结果的run.sh当作Boltz补跑入口。

旧原始模型找不到时，继续可执行的Full主任务；旧构建只计算现有材料确实支持的对照项，其余标记缺失。不得为了凑齐90个模型擅自重跑旧60个，也不得声称90套原始结果均已验收。

## 四、生成Full输入与MSA，保持旧协议

链固定为A=Full181抗原、H=原VH119、L=原VL106。核验抗体FASTA、IMGT/CDR注释与旧运行使用的是同一版本。独立ABodyBuilder2结构不是本轮Boltz共折叠的前置条件，不能为等它而另开抗体预测路线。

新增两种YAML：
- full_blind.yaml：A/H/L三条链均明确msa: empty。
- full_blind_msa.yaml：引用Full这一批次对应的A/H/L正式MSA文件。

正式MSA沿用已审计的ColabFold/MMseqs2流程与paired/unpaired组合方法；扩展generate_target以支持Full，不重新覆盖Native/Immunogen的MSA。

特别注意：现有MSA CSV含有配对key，不能只为Full生成一个新A.csv，再把旧H/L.csv随意拼上；不同批次相同数字key不自动代表同一配对关系。默认将Full、VH、VL作为一组按原生成流程处理。只有能证明原始单链搜索缓存可复用、且正确重建本组配对键时才复用缓存。不要把Native抗原MSA简单补26列后冒充Full检索结果。

保存原始检索产物、三条链最终CSV、query核对结果、搜索时间、工具版本、MSA原始深度、配对信息和SHA256。MSA含小写插入时，按A3M规则核对query/比对列；不能用原始字符串长度直接判定所有同源序列错误。区分原始MSA深度、CLI允许上限和实际输入截取量；无法核实实际截取量时明确记为未知，不编造。

新旧MSA来自不同检索时间/数据库快照时，在可比性限制中记录；不能说除抗原长度外一切完全相同。API失败时保留错误并有限重试；可以先完成empty分支，但不得把MSA失败后的empty结果标成formal MSA。

## 五、正式GPU任务：固定30个，不因结果好坏扩增

先核查实际环境。旧记录为Boltz软件2.2.1，使用的是boltz1模型。优先复用原环境及同一boltz1_conf.ckpt；记录版本、checkpoint SHA256和完整命令，禁止静默升级为Boltz-2或改换模型。

沿用原正式参数：
--model boltz1
--checkpoint <核验后的原boltz1_conf.ckpt>
--cache <核验后的原模型缓存>
--accelerator gpu --devices 1
--diffusion_samples 1
--recycling_steps 3
--sampling_steps 200
--no_kernels
--num_workers 0
--output_format pdb
--write_full_pae
--seed 1至15

具体CLI兼容性先在当前冻结环境中检查，不凭记忆编造新参数。模型或环境缺失时先定位已授权可用资源；不要擅自换检查点、创建付费云资源或绕过权限。

两组任务分别为：
Full181 / empty / seed1…15；
Full181 / formal_msa / seed1…15。

先用每组seed1按正式200-step参数验证输入—结构—置信度—分析链路，成功后该seed1计入正式15个，不额外重复跑pilot。后续按固定计划完成，不能因seed1表位或分数不理想就改变参数、提示或筛seed。禁止引入Phase1表位约束或受体信息；目前只有blind预测后的Phase1-informed排序。

检查空闲GPU后调度，若确有两张可用A6000，可各运行一条分支；否则单卡串行。不要照搬旧脚本强制CUDA_VISIBLE_DEVICES=0导致两任务挤同卡，不杀死其他项目进程。每个任务独立目录、独立日志，记录开始/结束、退出码及可测的显存和耗时。

断点续跑不能仅凭一个confidence JSON存在就跳过。至少确认本次输入/权重/参数指纹一致、PDB可解析、A/H/L序列及长度181/119/106正确、confidence字段可读取、PAE产物有效，再标记完成。记录完整任务清单，避免重复seed与重复模型。

## 六、补齐下游分析，先修编号再算结果

现有接触脚本把抗原编号硬编码为local+26，且部分分析/排序脚本只识别native和immunogen。必须将构建元数据/映射贯穿接触、consensus、Phase1重叠、几何QC、代表结构、亲和力和导出；不能只修其中一张表。

1. 复合物接触与表位
沿用抗原A对抗体H+L的重原子最小距离定义，<4.5 Å为主分析，<5.0 Å为扩展记录。过滤水、非蛋白HETATM和氢；不要把VH–VL内部接触当作抗原结合。
输出逐残基对距离、抗原local/project编号、抗体chain/local/IMGT/CDR或FR注释，以及每seed的epitope、paratope、CDR参与情况。保留无接触模型，不为提高稳定性将其删除。

2. 置信度与几何
每seed读取confidence_score、pTM、ipTM、complex pLDDT，并保存原始JSON/PAE。可由实际输出确定链对应关系时，另报A–H/A–L的链对置信信息；不能把VH–VL高置信度当作抗原–抗体结合已可靠。额外字段缺失时写NA，不反推虚构。
复用现有geometry QC规则并保留阈值，记录断链/严重冲突等警告；不通过QC的模型保留并标注，不在看到结果后选择性剔除。PDB与JSON里的pLDDT量纲分别核实，不混用0–1与0–100。

3. 多seed共识及代表模型
每个MSA条件单独分析15个seed，不将两种MSA混成一个30-seed共识。输出逐残基支持次数/支持率、seed间epitope Jaccard、均值/标准差，以及>=50%共识（完整15个seed时至少8个支持）。
支持率分母明确区分计划、完成、缺失；空表位模型是零接触观察，不是运行失败。两个空集合的Jaccard即使按旧代码记1，也不能解释成稳定结合，应单列无接触状态。
复用原三套事后排序：confidence；consensus；Phase1-informed。对Full两种条件各自选confidence和consensus代表，informed代表单独标明选取依据；并列时固定按seed编号处理。同一seed获多个称号可复用结构，不虚构多个独立模型。
Phase1 prior必须读取冻结文件并记录其来源/哈希，不能根据Full结果改范围。原冻结prior找不到时，标记对应排序阻塞，其余主分析继续；不能从报告大致区间临时拼出一个prior当成原件。

4. PRODIGY亲和力代理
沿用原PRODIGY版本及参数，对30个Full原始未松弛复合物计算A对H,L的评分。确认选择为抗原–完整Fv两侧，不能误算A–H、A–L平均或VH–VL结合。
保存每seed stdout/stderr、退出码、DeltaG、预测Kd、温度和接触数；旧脚本列名为prodigy_Kd_M_25C，实际温度必须核实。无有效界面或解析失败时保留记录并写NA/原因，禁止填0或虚构数值。
主比较保持原始结构对原始结构；本轮不新增松弛。所有能量/Kd标为计算代理，不写成实测亲和力，不据此宣称Full更强或已证明中和。

## 七、与旧Native155、Immunogen134进行三构建比较

所有比较在相同MSA条件内进行；新旧结果分源保存，不修改旧基线。明确哪些来自原始重算、哪些来自已保存逐模型表、哪些只有历史汇总。

统一输出模型数、置信度、几何警告、CDR参与、表位支持率、Jaccard和PRODIGY分布。不要拿Full最好一个seed对旧构建平均值，也不要把15个计算seed解释为15个独立生物学样本。若旧原始表支持，复算旧汇总作回归校验；不支持的指标标NA，不用代表PDB外推15个seed。

比较范围固定：
- 三构建共同区域：project27–160。
- Full与Native的扩展共同区域：project27–181。
- 各自完整构建范围的结果另列，不混入共同区指标。

对每个Full模型单列1–26、27–160、161–181三个区段的接触残基数及其在完整epitope中的占比，并保存共同区的表位集合。用于回答新增N端或C端是否改变预测界面；不能预先规定Full必须重复旧132–150共识。

跨构建逐残基表中，未包含在某构建的残基必须为NA/不适用，而不是0接触支持；构建内确无接触才记0。共有区Jaccard先取交集范围再算。无表位状态另报，避免将空集合一致当作结合证据。

若做低成本结构叠合，只按共同project残基对齐：Full的project27–160对应local27–160，旧两构建对应local1–134。缺原始结构则只做可验收代表比较；不能把Full local1–134与Native local1–134直接配对。

结论必须区分：构建敏感性、预测界面共识、结构置信度与计算亲和力代理。Full181是保留项目N端的计算构建，不直接等同于生理成熟蛋白的实验结构。不得把新预测当作真实表位、实测KD或受体竞争机制已经确认。

## 八、交付、验收与停止条件

所有新增文件集中在补跑分支，建议至少包括：
- input_audit.json、input_sha256.csv、numbering_map_full181.csv、mapping_tests.log；
- Full的两种YAML、正式MSA及manifest、checkpoint/environment指纹；
- task_manifest.csv：30个任务的条件、seed、输入hash、输出路径、状态、失败原因；
- 30套可验收原始PDB、confidence、PAE和日志；未完成必须按实际数量报告；
- full181_model_summary.csv、full181_interface_contacts.csv、full181_epitope_support.csv；
- full181_model_ranking.csv、full181_prodigy.csv、full181_geometry_qc.csv；
- three_construct_comparison.csv、shared_region_epitope_comparison.csv、terminal_contact_summary.csv；
- representatives/及selection_manifest.csv；
- Full181_supplement_report.md、handoff_for_plotting.md；
- 本次新增/修改脚本、可复现入口、参数记录及最终SHA256清单。

预留绘图用长表和可追溯PDB即可，不渲染结构图、不生成统计图、不润色图注、不开展全项目大重构。交接文件说明每张未来可用表的字段、坐标、条件、单位、缺失和代表模型选择依据。

验收必须确认：181长度与末端L保留；抗体未变；30个任务无重复；两种MSA不混淆；所有表位project编号正确；预测、接触、能量与排序逐seed能关联；旧文件hash不变；结果数不虚报。旧脚本先跨构建取Top10再拆分的逻辑不能用于本轮交付，各构建的完整表必须分别保留。

开始时简报实际找到的根目录、原始输入/权重/旧结果可用性、将修改的最小脚本集合、两组15-seed任务安排，然后继续实际执行，不要只写计划。运行中维护progress.md和错误记录；硬阻塞时停受影响环节，说明具体缺失资源及已完成内容，不用新模型/假结果绕过。

最终报告用实际结果回答：Full完成多少/30，置信度和表位分布如何，与旧构建共同区是否一致，N端1–26是否参与界面，PRODIGY是否产生可用代理，哪些代表结构可供下一步作图，哪些扩展仍未做。结果差、无共识或无有效接口都照实交付，不追加seed追逐好结果。
