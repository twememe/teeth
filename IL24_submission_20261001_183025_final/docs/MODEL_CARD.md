# Model Card

适用范围：IL-24计算区段/表位候选、已给定抗体的结构构象与界面比较，以及独立单序列二分类原型。所有结果是计算证据，不是实验效能验证。

## 阶段一与构建

参考181-aa序列固定，项目1–181对应小鼠UniProt Q925S4的40–220。第31位项目S与UniProt F的既有差异保留。N端切点26|27结合UniProt、已有DeepSig和疏水性；C端构建边界来自冻结的疏水/结构规则。BepiPred-3.0与DiscoTope-3.0为预训练推理；percentile平均并列秩，`FinalB=(((B+D)/2)+R+C)/3`。pLDDT仅为结构置信注释。15/20/25窗口、重叠>50%和既有Top5规则不变。历史100次保守性bootstrap保留，本轮未新增。

## 阶段二

ImmuneBuilder 1.2的ABodyBuilder2预测Fv；Boltz包2.2.1调用Boltz-1权重。正式Native155、Immunogen134各有empty-MSA和formal-MSA条件各15 seeds；Full181补充同样分两组各15。参数为3 recycling、200 sampling、1 diffusion sample。MSA敏感性3 seeds/50步和N74单NAG敏感性3 seeds/200步单独记录，不与正式组混用。主清单使用正式MSA、未松弛PDB与同一组PRODIGY 2.4.0预测dG；糖基化分支只属探索。

**序列版本限定（集中记录）**：当前结构和输入采用VH119/VL106；历史测序表120/107差异未定案。本版不改序列、不混用结构、不发起重测序或重预测。阶段一/二保持既有冻结先验，ESM不参与排序。

ABodyBuilder2 `error_estimates.npy`是ensemble平均平方位置偏差（Å²），平方根为预测误差（Å）。ImmuneBuilder 1.2在细化输出中将平方根误差写入PDB B列。随包原始rank0–3未细化PDB的B列为0，不能解释为pLDDT，也不使用50/70阈值。新增误差注释CSV仅转换上述原数组单位，未改变坐标或排序。

6DF3的C链是人IL-24，L链是IL-22RA1，H链是IL-20RB；这里L/H不是抗体。修正版只用氨基酸蛋白重原子接触（5Å），将IL-24侧受体足迹R和预测抗体接触表位E映射到同一项目编号，计算`|E∩R|/|E∪R|`。155位对齐中106位相同，约68.4%，48个受体接触位点可映射；T198对应项目173。跨物种位点转移不意味着相同侧链作用；6DF3不能证明IL-20R1/IL-20R2机制。相关报告、CSV、JSON、SVG来自同一次修正版计算。

## 独立ESM-2原型

实际输入只有`model_sequence`；没有抗原条件化、WT/突变成对接口、亲和力/突变效应预测或IL-24实验反馈闭环。ESM-2 650M，33层，UR50D；0–15层query/value注入LoRA（r8、alpha16、dropout0.1），均值池化+线性分类头，阈值0.5。官方基座冻结且未合并适配器；本项目训练LoRA与分类头。包中保留原样项目checkpoint、tokenizer、配置与`models/base_manifest.json`，通过固定官方文件重建并核对基座张量身份；大基座不随ZIP，详见BASE_MODEL.md。

冻结基线MCC=0.385864，LoRA MCC=0.499083；同一20,000条验证集也参与checkpoint选点，因此不是独立确认结果。原500次分组bootstrap区间只是既有条件性不确定性，不能消除选择偏差。测试33,508条未用于该模型选择/主报告。完全相同序列跨划分隔离，但未完成严格家族/抗原/近同源隔离。混合来源标签高值不代表一致生物学方向；62,233个相同输入序列组存在实验依赖标签冲突。

历史训练在step10000后因rotary缓存问题中断；当前正式代码使用no_grad修复。恢复保留权重和已消费批次，原checkpoint未含AdamW/RNG状态，二者重置，因此重新从头训练不保证逐位重现原权重。交付的最佳权重原样冻结；本轮只验证加载和8序列推理。
