# 第三方工具与模型

| 工具 | 本项目实际版本/获取 | 使用与许可材料 |
|---|---|---|
| ESM-2/fair-esm | esm2_t33_650M_UR50D；官方源文件版本/SHA及转换后张量指纹见models/base_manifest.json；fair-esm 2.0.0用于ESM模块 | Meta官方https://github.com/facebookresearch/esm ，MIT；UR50D预训练，无项目预训练 |
| BepiPred-3.0 | 98bedf2471e04766c981585180c20da7961c6804；bp3 0.0.12.7；原fair-esm 1.0.3 | https://github.com/UberClifford/BepiPred-3.0；五折官方权重与默认vt_pred/阈值0.1512/9-aa平滑；固定CLI README为学术非商业使用；bp3 0.0.12.7发行包自带MIT，十份随包权重逐字节匹配；CLI赛事交付范围待确认。证据见licenses/final_20261001/README.md |
| DiscoTope-3.0 | 35d9f2e55f97eaba2a7acefbc394db58fb9670bc；100 XGBoost+2 GAM；ESM-IF1 | https://github.com/DTU/DiscoTope-3.0；原README/LICENSE随包；AlphaFold模式，官方0.90阈值，排序用未校准连续分数 |
| ImmuneBuilder ABodyBuilder2 | 1.2；4模型Zenodo7258553 | https://github.com/oxpig/ImmuneBuilder；原MIT许可及Zenodo模型来源，模型SHA独立记录 |
| Boltz | Python包2.2.1，Boltz-1 boltz1_conf.ckpt | https://github.com/jwohlwend/boltz；代码MIT；模型按官方资产卡许可，不统一重新许可；--model boltz1 |
| PRODIGY | prodigy-prot 2.4.0 | https://github.com/haddocking/prodigy；原Apache-2.0许可，A与H,L链，25°C，结果相对比较 |
| Transformers/PEFT | 4.48.3 / 0.14.0 | Hugging Face Apache-2.0；ESM-LoRA实现，未使用弃用的Ab-Tune依赖方案 |
| FreeSASA / DSSP / MAFFT | 2.2.1 / 4.6.1 / 原phase1 env锁定版本 | 表面/二级结构/序列比对；各自原许可与来源，环境规格保留 |
| DeepSig | 0.9；Python3.8.20/TF2.2.0/Keras2.4.3 | https://github.com/BolognaBiocomp/deepsig；仅既有N端证据和实际读取/分析入口，获取说明见FULL_RUN |
| ColabFold MMseqs2服务 | Boltz2.2.1集成调用；服务器数据库精确快照未记录 | https://api.colabfold.com；真实MSA CSV与请求元数据随包，重请求可能受服务器数据库变化影响 |

许可文本收集于`docs/licenses`与对应vendor/third_party目录。部分数据与模型卡的许可未在历史记录中完整抄录时，按DATA_SOURCES中的缺口明确披露；不向第三方资产统一添加项目自选许可证。未使用商业LLM/API提示词生成本项目候选；ColabFold请求序列/参数记录在MSA manifests中。模型调用时间见原日志（阶段一2026-08-21/22、区段分支2026-08-27、阶段二2026-08/09、Full181 2026-09-21、ESM 2026-09-22）；缺少精确调用时刻的资产标记未记录。

2026-10-01仅补核上述BepiPred与AbAgym两项，具体来源、资产范围、版权声明及待补材料见 `licenses/final_20261001/README.md`。
