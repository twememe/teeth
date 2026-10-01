# Full181 阶段二补跑验收报告

日期：2026-09-21。**本轮新增正式模型 30/30，全部通过原始产物与本地回收验收；无失败 seed、无重复预测。** empty 与 formal MSA 各 seeds 1–15，seed1 为正式 200-step 任务并计入15。旧 Native155/Immunogen134 的60个原始模型只读核验，旧GPU任务重跑数为0。

## 工程、协议与输入

- 服务器保持 172.16.149.161（xinlab-Super-Server）。用户授权 root 后找到原工程 `/home/xinlab/home/ws/phase2`。
- 复用原 Python 3.10.12、Boltz 软件2.2.1、boltz1 模型、torch2.6.0+cu124、PRODIGY2.4.0；未替换或升级模型。原 checkpoint 为3,595,352,714字节，CCD为345,859,128字节，路径/大小/时间戳在环境与逐seed指纹中。
- `/home` 当时仅剩770MiB，使用同机 `/mnt/il24_full181_supplement_20260921` 实际执行，两张A6000分别运行empty与formal MSA。按用户追加要求，统一归档根为 `/mnt/pxh/teeth/`，含 `phase1ab` 和 `phase2`；补跑工作包位于 `phase2/full181_supplement`。原目录保留以维持历史追溯。
- Full为181 aa，H180/L181保留，上传FASTA、项目FASTA及参考链A逐位一致；VH119/VL106未变。Native/Immunogen分别为project27–181与27–160。Full使用offset0；旧构建使用经序列验证的显式映射，不能沿用Full local+26。9项必需编号测试通过。
- 正式参数固定：boltz1；1个diffusion sample；3 recycling；200 sampling steps；no_kernels；num_workers0；PDB；完整PAE；seeds1–15。每任务有独立原始输出、完整命令、stdout/stderr、时间和GPU采样。峰值采样显存12,859MiB。
- 仅增加补跑、显式坐标、分析、原始文件核验和交付入口；旧工程脚本和结果不覆盖。最初pxh权限阻塞、传输暂停及恢复记录保留，最终任务失败数为0。

## MSA与溯源

Full/A、H、L按同一批ColabFold/MMseqs2流程检索，UTC 2026-09-21 08:53:17–08:54:35。三链paired原始深度均977；unpaired分别712、9978、10335；最终CSV分别1588、10953、11310。每条CSV的key/sequence均从本批A3M直接重建核对，未拼入旧H/L批次。

实际15个MSA任务的processed数组均与原Boltz CSV解析器对本批输入的结果逐项一致：A1251、H8192、L8192条。与最终CSV的差异来自原解析器去重和CLI8192上限。网络内部联合/抽样张量深度未记录，明确为未知，不等同于上述缓存深度。原始A3M、API响应、query、检索记录和processed文件均保留。新旧MSA检索时间不同，数据库快照不保证相同，不能把差异全部归因于抗原长度。

按用户直接要求未计算哈希，也未做冒烟/pilot预测。指纹采用路径、大小、mtime及输入/数组直接比较；原本地96个文件、原服务器1166个文件元数据未变，不作密码学内容不变声明。

## Full181结果

| 条件 | 模型数 | confidence均值 | ipTM均值 | complex pLDDT均值(0–1) | 表位残基数均值 | seed间Jaccard均值 |
|---|---:|---:|---:|---:|---:|---:|
| empty | 15 | 0.6520 | 0.4704 | 0.6975 | 26.0000 | 0.3852 |
| formal_msa | 15 | 0.7538 | 0.4861 | 0.8207 | 16.9333 | 0.4278 |

主接触定义为A对H+L蛋白重原子最小距离严格<4.5Å，<5.0Å另存扩展记录；VH–VL内部接触不计入抗原界面。30个模型均有接触观察，2232条<5Å残基对、420条CDR/FR区域记录已导出。

每个条件独立15个seed，≥50%共识阈值为至少8个支持。共识project位置如下：

- empty：56;57;58;59;60;62;63;64;65;98;99;100;101;103;104;107;108;110;111;114;118。
- formal_msa：73;132;133;135;136;139;140;142;143;144;146;147;149;150。

A–H/A–L链对置信度来自实际processed链名到chain_id的映射，均值如下；整体ipTM或Fv内部置信度不能替代抗原–抗体界面的置信度。

- empty：A→H ipTM 0.1750，A→L ipTM 0.1712。这些链对数值仍较低，不能据整体结构分数确认结合界面可靠。
- formal_msa：A→H ipTM 0.1979，A→L ipTM 0.1912。这些链对数值仍较低，不能据整体结构分数确认结合界面可靠。

PDB B-factor的pLDDT为0–100，JSON confidence/pTM/ipTM/complex pLDDT为0–1；PAE均为有效406×406、有限且非负数组。

## 末端与CDR参与

| 条件 | N端1–26有接触seed | C端161–181有接触seed | N端接触数范围 | C端接触数范围 |
|---|---:|---:|---|---|
| empty | 7/15 | 8/15 | 0–7 | 0–2 |
| formal_msa | 1/15 | 1/15 | 0–1 | 0–1 |

N端在部分模型参与预测界面，两种条件均没有N端残基达到8/15共识。每个seed的1–26、27–160、161–181计数及占完整表位比例保存在terminal_contact_summary.csv。

empty中H的3个CDR及L的3个CDR均在15/15模型参与；formal MSA中H的3个CDR和L-CDR3为15/15，L-CDR1为6/15，L-CDR2为8/15。FR接触另列，未把FR归入CDR。

## 三构建共同区

所有比较在相同MSA条件内进行，并先裁剪到共同project坐标。以下为≥50%共识集合的Jaccard；完整表另含225个跨seed组合，不将seed当作生物学配对或独立实验样本。

| 条件 | 构建对 | 共同区 | 共识Jaccard |
|---|---|---|---:|
| empty | Full181 / Native155 | project27-160 | 0.7083 |
| empty | Full181 / Immunogen134 | project27-160 | 0.0312 |
| empty | Native155 / Immunogen134 | project27-160 | 0.0323 |
| empty | Full181 / Native155 | project27-181 | 0.6800 |
| formal_msa | Full181 / Native155 | project27-160 | 0.8571 |
| formal_msa | Full181 / Immunogen134 | project27-160 | 0.5714 |
| formal_msa | Native155 / Immunogen134 | project27-160 | 0.6667 |
| formal_msa | Full181 / Native155 | project27-181 | 0.8571 |

Full与Native的共同区共识在两种条件均有较大重叠；与Immunogen的重叠明显依赖MSA条件。三构建与MSA处理均体现构建/输入敏感性，不能把任一共识称为已验证真实表位。构建外残基为NA，不按0接触计算。

## PRODIGY与几何QC

30/30原始未松弛结构均以A对H,L获得PRODIGY2.4.0评分，实际stdout温度均25°C；退出码、接触数、ΔG、Kd、stdout/stderr和命令保留。

| 条件 | ΔG均值(kcal/mol) | ΔG范围 | Kd中位数(M) | Kd范围(M) |
|---|---:|---|---:|---|
| empty | -14.980 | -17.8–-12.2 | 5.80e-12 | 8.60e-14–1.20e-09 |
| formal_msa | -12.413 | -15.3–-10.2 | 1.90e-09 | 5.70e-12–3.20e-08 |

这些是计算亲和力代理，不是实测KD，也不据此判定Full更强或证明中和。旧60个PRODIGY数值复用原逐模型表；旧表未包含实际stdout温度，温度字段保留未知，不从25C列名反推实测输出。

30个Full和60个旧结构在既定阈值下均有原子冲突警告，全部保留。所有90个均无CA断链、主链缺失或肽键C–N距离异常。Full empty的跨链/非相邻链内冲突合计132/361；formal MSA为60/112。阈值为CA距离>4.5Å，C–N<1.1或>1.6Å，冲突半径1.8Å；没有以QC警告进行事后筛除或新增松弛。

## 代表结构

| 条件 | 选择依据 | seed | 本地PDB文件 |
|---|---|---:|---|
| empty | confidence | 2 | Full181_empty_seed2.pdb |
| empty | consensus | 6 | Full181_empty_seed6.pdb |
| empty | phase1_informed | 11 | Full181_empty_seed11.pdb |
| formal_msa | confidence | 6 | Full181_formal_msa_seed6.pdb |
| formal_msa | consensus | 9 | Full181_formal_msa_seed9.pdb |
| formal_msa | phase1_informed | 9 | Full181_formal_msa_seed9.pdb |

六个称号对应五个独立PDB；formal MSA seed9同时为consensus和Phase1-informed代表。confidence按分数；consensus按与同条件其他seed的平均Jaccard再按confidence；informed按原冻结62残基prior coverage再按confidence；最后均以数字seed升序打破并列。prior原件、来源和指纹已保存，未由Full结果重建。

## 验收与交付

- 30个唯一Full任务、30套PDB/confidence/PAE、30条PRODIGY有效记录全部验收；A181/H119/L106逐位正确，H180/L181保留，MSA状态未混用。
- 旧60套PDB/confidence/PAE已实际核验，4459条旧接触判定与保存表一致；2条距离有最后一位舍入差异，最大0.00050029Å，记录在legacy_raw_audit，未影响阈值/表位/CDR。44项旧汇总回归核对全部通过。
- prediction、contacts、QC、PRODIGY、ranking按construct+condition+seed关联。原run_state未因迁移改写；artifact_path_map.csv和selection_manifest_local.csv提供原执行路径、本地路径及teeth归档路径。
- 最终数据用途、单位、NA规则、代表选择与文件入口见handoff_for_plotting.md；逐表字段/行数见table_schemas.json；逐seed与旧文件验收见delivery_verification.json。
- 未生成任何图、图注或图结论。未补跑SignalP/DeepSig、阶段一、NAG、松弛、6DF3修复或其他模型/实验。这些扩展仍未补齐，主线完成不代表全部扩展完成。

Full181保留项目N端，是本项目计算构建，不直接等同于生理成熟蛋白的实验结构。
