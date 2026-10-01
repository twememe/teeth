# 两项许可补核（2026-10-01）

本记录只核对当前包实际使用的 AbAgym 与 BepiPred；没有为第三方资产重新许可，也未联系作者或编造授权。以下核查日期是本轮查阅日期，不是历史下载日期。

## AbAgym 数据

固定快照：87af82ad68bb90921e15fb8b79c7bcb1d05b23d3。证据为同目录 AbAgym_README.md、AbAgym_tree.json，原始链接与完整 SHA 见 sources.json。
原始仓库：https://github.com/3BioCompBio/AbAgym/tree/87af82ad68bb90921e15fb8b79c7bcb1d05b23d3 。该快照 README 的 License 明确限非商业使用；商业应用（包括出售基于该数据训练的预测工具访问权）须联系 Fabrizio.Pucci@ulb.be。不是“没有任何许可说明”，也不是通用开源许可证。
适用资产为原始 interface 数据、metadata/PDB，以及包含该数据的项目加工记录和 train/validation/test 划分中的 AbAgym 行。这里没有单独的软件或预训练权重许可；本项目训练权重不能被描述为官方模型。项目没有打包该 README 另行链接的 FoldX benchmark 衍生结构，故不能把 FoldX 的额外授权视为已获得。
署名：G. Cia, D. Li, S. Poblete, M. Rooman, F. Pucci, “AbAgym: a well-curated dataset for the mutational analysis of antibody–antigen complexes”, mAbs 17(1), 2025, https://doi.org/10.1080/19420862.2025.2592421 。论文开放获取许可不替代数据条款。
该文字明确使用范围，但未明确授权向赛事组织方再分发原始与加工数据。未发现项目内额外授权证明。本包保持历史数据、划分、权重，未移走原始文件后假称加工文件不受影响。评测材料现存且可运行；对外提交的授权状态仍待团队确认。

## BepiPred 源码、模型与数据分别处理

CLI 固定快照：98bedf2471e04766c981585180c20da7961c6804，https://github.com/UberClifford/BepiPred-3.0/tree/98bedf2471e04766c981585180c20da7961c6804 。BepiPred_CLI_README.md 明确代码和数据可由学术团体作非商业用途；营利应用需另获许可（Morten Nielsen, morni@dtu.dk）。保留原说明和署名，不将其改标 MIT；其文字未明确赛事再分发权。

实际 bp3 库为 PyPI 的 0.0.12.7，而非仅凭 CLI 仓库名称判断。官方该版本 wheel 内确有完整 MIT License（Copyright 2023 Joakim Clifford），见 bp3-0.0.12.7_LICENSE。wheel SHA256 为 6735063426e09db34136f008c7eb76ae243ff1cdc2a63cd984cea6d8e12e880d。原始分发见 https://pypi.org/project/bp3/0.0.12.7/ 。本轮下载该 wheel 后核对随包 bp3 文件，完整逐文件记录见 bp3_package_matches.json；两套各五折共十个 .pt 权重均与该 MIT 分发包逐字节一致，见 weight_matches.json。该发行包中的代码和模型按随发行版的 MIT 文字保留版权和许可声明；这不是将新仓库许可证倒套到旧快照。CLI 仓库说明与 bp3 分发许可分别保留，不将 bp3 的 MIT 扩张到 CLI 或上游训练数据。

BepiPred 上游训练数据并未作为该工具训练集在本项目重新获取或打包；项目包含的是 IL-24 上的历史预测证据，仍标明来源与学术非商业用途。本轮不重跑 BepiPred。引用：Clifford et al., “BepiPred-3.0: Improved B-cell epitope prediction using protein language models”, 2022, https://doi.org/10.1002/pro.4497 。
另查 DTU 官方下载渠道 https://services.healthtech.dtu.dk/cgi-bin/sw_request?packageversion=3.0b&platform=src&software=bepipred&version=3.0 ，该渠道有独立学术软件协议并限制向第三方分发。它不是本包已证明来源的 PyPI wheel/固定 GitHub 快照，未提交表单或接受该协议，不能未经证明将其当成历史获取协议，也不能忽略不同渠道条款的差异。

## 唯一待团队补充的材料项

请提供一份本次赛事交付的许可/授权确认材料：说明团队与用途符合学术、非商业条件，并确认 AbAgym 原始及加工数据和 BepiPred CLI/相关预测证据可向组委交付；如现有协议未涵盖该范围，需来源方书面授权或可核查的适用条款。bp3 0.0.12.7 的 MIT 文件与十份权重证据已补齐，不将其再次列为无许可权重。
在该材料补齐前，本包可报告技术运行测试的真实结果，但不宣称“完整附件5材料全部通过”或已获组委验收。下载链接本身不是授权。无明确证据表明本包实际 PyPI/GitHub 来源资产禁止本次再分发，因此本轮未擅删资产、重切数据或改权重；待确认状态明确随包交付。
