# Step 3 — BepiPred-3.0

## Step 3 完成

### 1. 本步做了什么
对冻结的 181-aa 项目参考运行官方 standalone BepiPred-3.0，并保存逐 residue 正类概率、官方默认阈值分类和官方 9-aa rolling mean。

### 2. 使用的算法
BepiPred-3.0 official standalone，`vt_pred`，默认阈值 0.1512，默认 rolling window 9。

### 3. 算法属于什么
基于预训练 ESM-2 protein language model 表征的深度学习推理；没有训练或微调模型。

### 4. 算法原理
通俗地说，ESM-2 从大量蛋白序列中学习氨基酸“语言规律”，BepiPred 再判断每个位置多像已知抗体接触 residue。技术上，官方五折前馈网络 ensemble 对 ESM-2 residue embedding 输出正类概率，`vt_pred` 在平均正类概率上应用固定阈值。

### 5. 为什么本项目需要它
它提供独立的 1D sequence-based epitope evidence；后续融合应优先使用连续概率，而不是只使用 0/1 分类。

### 6. 输入
- `inputs/reference_181.fasta`

### 7. 输出
- `results/03_bepipred/raw/`
- `results/03_bepipred/bepipred_residue_scores.csv`
- `results/03_bepipred/bepipred_profile.png`
- `reports/step03_bepipred.md`

### 8. 关键数值结果
- Residue rows：181。
- Positive probability：min=0.018547，mean=0.131211，max=0.313823。
- 官方默认阈值阳性 residue：71/181。
- 最高连续概率位置：project 53 (N)，score=0.313823。
- 高分 residue（Top 10）：53N (0.3138), 56Q (0.3062), 58Q (0.2955), 150R (0.2833), 158E (0.2752), 60D (0.2662), 49W (0.2596), 28E (0.2582), 30R (0.2565), 125Q (0.2543)。
- 官方阈值连续阳性区域：150–150 (n=1, mean=0.2833); 49–49 (n=1, mean=0.2596); 53–54 (n=2, mean=0.2521); 124–125 (n=2, mean=0.2466); 56–63 (n=8, mean=0.2435); 73–73 (n=1, mean=0.2361); 153–154 (n=2, mean=0.2316); 27–28 (n=2, mean=0.2281); 127–132 (n=6, mean=0.2186); 156–159 (n=4, mean=0.2014)。
- Runtime device result：cuda。
- Inference time：84.50 seconds。

### 9. 当前结果怎样解释
高分位置是 sequence-based computational evidence，不等于实验验证表位、受体界面或中和位点。9-aa rolling mean 仅用于线性连续性辅助展示。

### 10. 是否存在问题
模型默认阈值不是对本蛋白重新校准的阈值，不能用历史结果调参。完整官方运行日志位于 `logs/bepipred_inference.log`。

### 11. 下一步
Integration 将连续 BepiPred 概率与当前可用的 RSA 表面证据和 conservation 证据透明融合。DiscoTope 分支标记为缺失，不引入替代模型。
