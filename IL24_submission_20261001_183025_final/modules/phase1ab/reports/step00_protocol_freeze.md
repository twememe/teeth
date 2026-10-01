# Step 0 — Protocol Freeze / Pre-registration

## Step 0 完成

### 1. 本步做了什么

冻结唯一 181-aa 项目参考、坐标系、允许/禁止数据、算法和候选提取规则；复制允许的 PDB 与 UniProt gzip 到 `inputs/` 并核验 SHA-256；隔离未授权 JPG；建立 WSL2/Conda 核心环境、依赖锁定和方法注册。

### 2. 使用的算法/工具

SHA-256 文件完整性校验；WSL2 2.6.3.0；Miniconda/Conda 26.5.3；Python 3.12.13；环境 YAML、build-pinned 与 explicit package export。

### 3. 算法属于什么

这是可复现性工程、文件完整性校验和环境管理，不是深度学习或生物学打分。

### 4. 原理

通俗地说，先把“研究对象、尺子和规则”封存并盖章，后续任何结果都只能按同一把尺子重算。技术上用序列/文件哈希锁定输入，用 Conda 包清单和精确版本锁定运行时，用协议文件锁定权重、窗口和停止条件。

### 5. 为什么本项目需要它

防止坐标漂移到 UniProt 220-aa 编号、防止历史答案泄漏、防止结果出来后调整阈值/权重，也保证结构、序列、进化与模型证据可以追溯重算。

### 6. 输入

- `AF-Q925S4-F1-model_v6.pdb`
- `uniprotkb_accession_Q925S4_2026_08_20.fasta(1).gz`
- 用户冻结的 181-aa reference

### 7. 输出

- `protocol.yaml`
- `inputs/manifest.json`
- `inputs/reference_181.fasta`
- `README.md`
- `logs/environment.txt`
- `methods_registry.csv`
- `envs/il24-core*.yml` 与 explicit/package exports

### 8. 关键数值结果

- Frozen reference：181 aa；sequence SHA-256 `3642351ae1bb5ca0788598425de59f455672266ed7f434d21ae1bfe654992374`。
- PDB SHA-256：`9d722f2c1fdb78c6534cad517a66fe6f7cdd142903d6364e05957ca7ca42a224`。
- UniProt gzip SHA-256：`45d2ce6b94f61c8c202f7b8a84631a4b02e5b87c6da066698968bdc8b1e721d8`。
- WSL2 Ubuntu 26.04、kernel 6.6.87.2；RTX 4050 6141 MiB、driver 581.08。
- Core：Python 3.12.13、Biopython 1.88、FreeSASA 2.2.1、mkdssp 4.6.1、MAFFT 7.526。

### 9. 当前结果怎样解释

这一步只证明输入、坐标、运行环境与方法边界被冻结，不提供任何抗原区或表位结论。

### 10. 是否存在问题

官方 Miniforge 下载受 GitHub release-asset DNS 阻塞，已改用官方校验 Miniconda；NTFS 解包大小写冲突使运行时改到 WSL ext4；DSSP ABI 通过精确 `libmcfp=2.0.1=h171cf75_0` 修复。所有失败和修正均保留日志。

### 11. 下一步

按冻结协议运行序列 QC、结构 QC、官方 BepiPred/DiscoTope、FreeSASA/DSSP 和进化保守性；历史抗原区仍保持隔离。

