# Step V2.1-1 — N-terminal boundary

## SignalP 6
未执行；官方 portable package 和权重不可用，状态为 `UNAVAILABLE_EXTERNAL_LICENSED_PACKAGE`。

## Reviewed UniProt annotation
数据库专家人工审阅的加工注释映射到 project signal peptide 1-26、cleavage 26|27。

## DeepSig
BolognaBiocomp 的深度学习 signal-peptide predictor；正式 euk inference 预测 SP 1-26、confidence 1.0、Chain 27-181。

## N-terminal Kyte-Doolittle
冻结的多尺度曲线显示 N 端疏水核心在 24-28 附近转为亲水，作为 supporting evidence。

## 输入
181-aa FASTA/PDB and UniProt Q925S4 annotation.

## 输出
`n_terminal_boundary.csv` and `n_terminal_hydropathy.png`.

## 实际数值
Construct start=27; UniProt=26|27; DeepSig=26|27; evidence=STRONGLY_CONCORDANT.

## 本步起什么作用
Excludes the secretion signal without optimization.

## 是否有异常
SignalP 6 unavailable，但不再阻塞；DeepSig 与 UniProt 精确一致。Project position 31 的差异位于 signal segment 之外。
