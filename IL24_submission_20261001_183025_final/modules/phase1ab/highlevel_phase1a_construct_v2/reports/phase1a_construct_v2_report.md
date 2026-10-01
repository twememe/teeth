# IL-24 High-level Phase 1A v2.1
## Recombinant Immunogen Construct Selection

## 1. Executive Summary

Recommended construct: **27-160** (134 aa). Start 27 is supported by reviewed UniProt, DeepSig SP 1-26 (confidence 1.0), and N-terminal hydropathy. End 160 is the rightmost q75-acceptable endpoint in the stable hydrophobic band 154-160. Endpoint 160 sensitivity is **ROBUST** (1.000 acceptable).

## 2. Input integrity

FASTA and PDB Chain A are 181 residues with exact 181/181 identity and project numbering 1-181.

## 3. N-terminal signal peptide

- SignalP 6: not executed; `UNAVAILABLE_EXTERNAL_LICENSED_PACKAGE` (official portable package and weights absent).
- Reviewed UniProt Q925S4: project signal peptide 1-26; cleavage 26|27.
- DeepSig 0.9: euk mode; SP 1-26; cleavage 26|27; confidence 1.0; mature Chain 27-181.
- Concordance: STRONGLY_CONCORDANT; final start=27.
- N-terminal hydropathy: hydrophobic transport-label core followed by a polar transition around residues 24-28.

UniProt annotation 是数据库专家结合实验、序列和文献整理的人工审阅加工注释。DeepSig 是识别蛋白 N 端‘分泌运输标签’的深度学习网络，在本项目中只用于独立验证 signal region，不是 SignalP。

## 4. C-terminal hydropathy

All primary endpoints used complete symmetric 7/9/11/15/21-aa windows. The primary stable band is [154, 160] with peak 158 and peak HydroConsensus 0.991228.

| window | DeltaH_at_peak_158 | percentile |
| --- | --- | --- |
| 7 | 3.3143 | 1.0000 |
| 9 | 2.4444 | 0.9825 |
| 11 | 2.2182 | 1.0000 |
| 15 | 1.9933 | 0.9912 |
| 21 | 1.6667 | 0.9912 |

Kyte-Doolittle 给每个氨基酸一个偏水/偏油数值，再用滑动窗口寻找性质变化。完整对称窗口要求切点两侧都有完整数据，避免靠近右端时因缺失窗口产生人为右移。

## 5. 3D structural boundary

Fixed flanks are e-10..e and e+1..e+10; heavy-atom contacts use <5.0 Å. HIGH requires at least two of: density >= q75, continuous helix/strand, or mean RSA <0.10. pLDDT is confidence-only.

| endpoint | local_cross_contact_density | dssp_status | boundary_mean_RSA | boundary_median_pLDDT | structural_risk |
| --- | --- | --- | --- | --- | --- |
| 154 | 0.05 | COIL_BOUNDARY | 0.5118 | 96.88 | ACCEPTABLE_STRUCTURAL_RISK |
| 155 | 0.07 | COIL_BOUNDARY | 0.6261 | 96.62 | ACCEPTABLE_STRUCTURAL_RISK |
| 156 | 0.08 | STRUCTURE_TRANSITION | 0.5834 | 96.0 | ACCEPTABLE_STRUCTURAL_RISK |
| 157 | 0.07 | CONTINUOUS_HELIX | 0.4782 | 95.38 | HIGH_STRUCTURAL_RISK |
| 158 | 0.07 | CONTINUOUS_HELIX | 0.4518 | 95.38 | HIGH_STRUCTURAL_RISK |
| 159 | 0.07 | CONTINUOUS_HELIX | 0.2991 | 95.38 | HIGH_STRUCTURAL_RISK |
| 160 | 0.06 | CONTINUOUS_HELIX | 0.1761 | 95.44 | ACCEPTABLE_STRUCTURAL_RISK |

Contact density 看切点两侧在三维空间中还有多少‘搭桥’；DSSP 判断是否剪断完整 α 螺旋/β 折叠；RSA 判断局部是否深埋核心；pLDDT 只表示 AlphaFold 对该局部结构的把握。

## 6. Contact threshold sensitivity

| endpoint | q60 | q65 | q70 | q75 | q80 | q85 | q90 | acceptable_fraction | sensitivity_level |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 154 | OK | OK | OK | OK | OK | OK | OK | 1.0 | ROBUST |
| 155 | OK | OK | OK | OK | OK | OK | OK | 1.0 | ROBUST |
| 156 | OK | OK | OK | OK | OK | OK | OK | 1.0 | ROBUST |
| 157 | HIGH | HIGH | HIGH | HIGH | HIGH | OK | OK | 0.2857142857142857 | FRAGILE |
| 158 | HIGH | HIGH | HIGH | HIGH | HIGH | OK | OK | 0.2857142857142857 | FRAGILE |
| 159 | HIGH | HIGH | HIGH | HIGH | HIGH | OK | OK | 0.2857142857142857 | FRAGILE |
| 160 | OK | OK | OK | OK | OK | OK | OK | 1.0 | ROBUST |

Endpoint 160 is acceptable in 7/7 settings, therefore sensitivity is **ROBUST**. This analysis changes confidence only; q75 remains the primary veto.

## 7. Final engineering decision

Reviewed UniProt + DeepSig + N-terminal hydropathy → start 27; complete-window multiscale hydropathy → band 154-160; primary q75 3D veto → acceptable 154,155,156,160; maximal retention → end 160.

Maximal retention 的意思是：多个相邻位置都安全时选择靠后的一个，以保留更多 IL-24 抗原序列。

## 8. Why not shorter?

Structurally acceptable endpoints inside the primary band: 154, 155, 156, 160. Length was used only here as the final tie-break; the selected endpoint is 160.

154、155、156 和 160 都在同一疏水边界带且 q75 可接受；选择最右端 160 是为了保留更多抗原表面。

## 9. Why not longer?

161 以后无法同时满足最大 21-aa 完整对称窗口，因此属于 edge exploratory，不与 154-160 在同一证据质量下比较。

## 10. Relationship to Phase 1B

Phase 1A 决定宽免疫原 construct；Phase 1B 在蛋白内部寻找 hotspot。现有 Phase 1B 热点结果未修改。

## 11. Limitations

结构来自 AlphaFold；DeepSig 是计算预测；constructability 不等于实验 expression yield；局部 3D risk 不等于真实 folding stability。最终仍需表达、纯化和免疫实验验证。
