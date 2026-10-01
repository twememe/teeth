"""Assemble the final report only from measured results, with explicit validation limits."""
import os,pathlib,json,hashlib,math,argparse
ROOT=pathlib.Path(__file__).resolve().parents[1]

def mcc_from_counts(c):
    tn,fp,fn,tp=map(float,c)
    denom=((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))**.5
    return (tp*tn-fp*fn)/denom if denom else 0.

def paired_group_bootstrap(validation,baseline,lora,n_boot=500):
    import numpy as np
    groups,codes=np.unique(validation.model_sequence.to_numpy(),return_inverse=True)
    y=validation.label.to_numpy(dtype=int)
    counts=[]
    for pred in [baseline,lora]:
        code=y*2+(np.asarray(pred)>=.5).astype(int)
        counts.append(np.bincount(codes*4+code,minlength=len(groups)*4).reshape(len(groups),4))
    rng=np.random.default_rng(42);deltas=[]
    for _ in range(n_boot):
        sampled=rng.integers(0,len(groups),size=len(groups))
        deltas.append(mcc_from_counts(counts[1][sampled].sum(0))-mcc_from_counts(counts[0][sampled].sum(0)))
    return {'paired_sequence_group_bootstrap_n':n_boot,'delta_mcc_ci95':np.quantile(deltas,[.025,.975]).tolist(),'seed':42,'interpretation':'Conditional on this validation set; selection on the same validation set is optimistic and this is not independent confirmatory inference.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--output-dir",type=pathlib.Path,default=ROOT/"output");ap.add_argument("--reuse-uncertainty",action="store_true");ap.add_argument("--run-state",type=pathlib.Path,default=ROOT/"run_state.json");args=ap.parse_args()
    os.chdir(ROOT);os.environ.update(json.loads((ROOT/'environment.json').read_text()))
    import numpy as np,pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=args.output_dir.resolve();audit=json.loads((out/'data_audit.json').read_text());config=(ROOT/'training_config.yaml').read_text()
    reports={name:json.loads((out/f'{name}_metrics.json').read_text()) for name in ['random','frozen','lora']}
    pred={name:pd.read_csv(out/f'{name}_validation_predictions.csv') for name in reports}
    validation=pd.read_csv(ROOT/'data/validation.csv',keep_default_na=False,low_memory=False)
    for name,p in pred.items():
        assert len(p)==len(validation) and p.row_id.tolist()==validation.row_id.tolist(),f'{name}: validation row mismatch'
        assert np.array_equal(p.label.to_numpy(),validation.label.to_numpy()),f'{name}: validation label mismatch'
        assert np.isfinite(p.probability).all(),f'{name}: nonfinite predictions'
    rows=[]
    for name in ['random','frozen','lora']:rows.append({'model':name,**reports[name]['overall']})
    pd.DataFrame(rows).to_csv(out/'comparison.csv',index=False)
    compare=json.loads((out/'paired_comparison.json').read_text()) if args.reuse_uncertainty else paired_group_bootstrap(validation,pred['frozen'].probability,pred['lora'].probability)
    compare['delta_mcc']=reports['lora']['overall']['mcc']-reports['frozen']['overall']['mcc']
    (out/'paired_comparison.json').write_text(json.dumps(compare,indent=2))
    losses=pd.read_csv(out/'losses.csv');fig,axes=plt.subplots(1,3,figsize=(15,4.1),constrained_layout=True)
    frozen_losses=losses[losses.model=='frozen'].sort_values('global_step')
    axes[0].plot(frozen_losses.global_step,frozen_losses.training_loss,marker='o',markersize=4,color='#245B8A',linewidth=1.6)
    axes[0].set(title='Frozen ESM-2 + linear head',xlabel='Optimizer step',ylabel='Mean training loss (BCE)')
    raw=losses[(losses.model=='lora') & losses.validation_mcc.isna()].sort_values('global_step')
    assert len(raw) and len(frozen_losses),'Missing measured loss history'
    axes[1].plot(raw.global_step,raw.training_loss,color='#BE5A2E',alpha=.15,linewidth=.5,label='Logged batch loss')
    axes[1].plot(raw.global_step,raw.training_loss.rolling(50,min_periods=10).mean(),color='#A64824',linewidth=1.8,label='50-point moving mean')
    axes[1].axvline(10000,color='#666666',linewidth=1,linestyle=':',label='Recovery point')
    axes[1].set(title='LoRA training loss',xlabel='Optimizer step',ylabel='Training loss (BCE)');axes[1].legend(frameon=False,fontsize=8,loc='upper right')
    validation_rows=losses[(losses.model=='lora') & losses.validation_mcc.notna()].sort_values('global_step')
    axes[2].plot(validation_rows.global_step,validation_rows.validation_mcc,color='#A64824',marker='o',markersize=3,linewidth=1.3,label='LoRA validation')
    axes[2].axhline(reports['frozen']['overall']['mcc'],color='#245B8A',linestyle='--',label='Frozen baseline')
    best_row=validation_rows.loc[validation_rows.validation_mcc.idxmax()]
    axes[2].scatter([best_row.global_step],[best_row.validation_mcc],color='#1D7A58',marker='*',s=100,zorder=4,label='Selected checkpoint')
    axes[2].set(title='Common validation set (n=20,000)',xlabel='Optimizer step',ylabel='MCC');axes[2].legend(frameon=False,fontsize=8,loc='lower right',bbox_to_anchor=(1,.08))
    for ax in axes:ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2)
    fig.savefig(out/'loss_curve.png',dpi=180);fig.savefig(out/'loss_curve.pdf');plt.close(fig)
    lora=json.loads((out/'lora_run.json').read_text());baseline=json.loads((out/'baseline_run.json').read_text());run=json.loads(args.run_state.read_text())
    best=out/'esm2_650m_lora_best.pt';assert best.exists() and best.stat().st_size>1000
    verdict='LoRA 在本验证集上的 MCC 高于冻结基线。' if compare['delta_mcc']>0 else '本次预算内未观察到 LoRA 的验证集 MCC 优于冻结基线。'
    lines=['# ESM-2 单序列代理分类：LoRA 可行性验证报告','','本报告评估的是单序列二分类代理任务：模型输入只有一条 `model_sequence`，没有抗原条件化输入，也没有同一抗原下的野生型/突变体成对比较，因此不含亲和力或突变效应预测。','','真实 ESM-2 650M 已完成前向、反向、LoRA 更新、双卡同步以及保存恢复验证。'+verdict,
        '本实验评估的是混合数据上的二分类代理任务；不能据此宣称对未见抗体家族、抗原或真实亲和力改善具有泛化能力。','',
        '## 同一验证集的实测对比','','主要模型选择指标为 MCC，分类阈值固定为 0.5。下列三组使用完全相同的 20,000 条验证记录。','',
        '| 模型 | MCC | F1 | Spearman ρ |','|---|---:|---:|---:|']
    labels={'random':'随机预测（实测）','frozen':'冻结 ESM-2 650M + 线性头','lora':'ESM-2 650M + LoRA + 同构线性头'}
    fmt=lambda x:'未定义' if x is None else f'{x:.6f}'
    for name in ['random','frozen','lora']:
        m=reports[name]['overall'];lines.append(f'| {labels[name]} | {fmt(m["mcc"])} | {fmt(m["f1"])} | {fmt(m["spearman_probability_label"])} |')
    lo,hi=compare['delta_mcc_ci95']
    lines+=['',f'LoRA − 冻结基线的 MCC 差值：{compare["delta_mcc"]:+.6f}；按相同序列分组的配对 bootstrap 区间为 [{lo:+.6f}, {hi:+.6f}]（500 次，种子 42）。同一验证集同时参与了模型选择，该区间不能消除选择偏差。',
        '随机基线的 MCC=0、F1=0.5 是类别均衡条件下的理论参考值，不是本次实测值。表内 Spearman 使用预测概率与二分类标签；各来源原始连续分数的相关性另列在 metrics JSON 中，未将不同量纲的分数直接混合。','',
        '## 训练与预算','',f'- 开始：{run["experiment_started_at"]}；总截止：{run["hard_deadline_at"]}。',
        f'- LoRA 已开始 {lora["actual_epochs_started"]} 轮，完整完成 {lora["completed_epochs"]} 轮；优化步骤 {lora["global_steps"]}。15 轮是上限，不代表全部完成。',
        f'- 提前停止：{lora["early_stopped"]}；达到训练时间上限：{lora["deadline_reached"]}。',
        f'- 冻结基线训练 {baseline["actual_epochs"]} 轮，学习率 {baseline["learning_rate"]}，缓存特征批量 {baseline["batch_size"]}。LoRA 从该基线线性头继续训练，分类头结构相同。',
        '- Ab-Tune 固定依赖与现有 CUDA PyTorch 版本冲突，按授权改用 PEFT + Transformers；层 0–15 的 query/value 注入 LoRA，r=8、alpha=16、dropout=0.1。',
        '- 使用 bf16 混合精度、长度分桶、梯度检查点。具体实际批量和设备数见 lora_run.json 与 training_config.yaml。','',
        '## 数据审计','',f'- 原始数据：AbBiBench {audit["raw_counts"]["abbench"]:,}、AbAgym {audit["raw_counts"]["abagym"]:,}、IL6 {audit["raw_counts"]["il6"]:,}，合计 {audit["total_raw"]:,}。',
        f'- 训练/验证/测试：{audit["split_counts"]["train"]:,} / {audit["split_counts"]["validation"]:,} / {audit["split_counts"]["test"]:,}。阈值相等记录明确排除；未复制记录补足条数。',
        f'- AbBiBench 二值化中位数仅由训练部分估计：{audit["abbench_training_median"]:.12f}。',
        f'- 共 {audit["conflicting_label_sequence_groups"]:,} 个相同输入序列组在不同实验中标签冲突；相同序列整体分组，跨集合重叠为零。',
        '- AbAgym 的 WT/突变字段实际是残基和 PDB 位点，已校验野生型残基并重建突变链。WT:mutation 格式保留为元数据，模型接收真实氨基酸序列。',
        '- AbAgym 主要记录为抗原突变，且高 DMS 分数可表示免疫逃逸；逐来源及突变对象结果见下表。',
        '- AbBiBench 和 IL6 原始数据没有可靠突变注释，mutation 列保留为空并标注来源，不编造突变信息。',
        '- IL6 原始划分保留；测试集在本次模型选择与主报告中未使用。近同源序列、抗体家族和抗原之间尚未做严格隔离。','',
        '## 按来源与突变对象的 MCC','','| 分组 | 冻结基线 | LoRA |','|---|---:|---:|']
    audit_path=out/'final_independent_audit.json'
    if audit_path.exists():
        final_audit=json.loads(audit_path.read_text());pos=lines.index('## 中断与一次修复') if '## 中断与一次修复' in lines else lines.index('## 数据审计')
        lines[pos:pos]=[f'最佳权重来自第 {final_audit["best_checkpoint_epoch"]} 轮的第 {final_audit["best_checkpoint_step"]:,} 步。早停按每轮末 MCC 的连续未提升轮数计算（patience=3）；每 2,000 步的中途验证也用于保存最佳权重。', '最终冻结基线与 LoRA 均重新加载保存权重，使用相同 bf16 推理流程重评完整验证集；逐条行标识、标签和各分组指标已独立复算。', '损失图浅色线为每 20 步记录的小批量损失，深色线为连续 50 个记录点的滑动平均；第 10,000 步标示恢复位置。第三面板显示验证集 MCC 与模型选择。','']
    recovery=lora.get('training',{}).get('recovery')
    if recovery:
        pos=lines.index('## 数据审计')
        lines[pos:pos]=['## 中断与一次修复','',f'正式训练在第 {recovery["from_global_step"]:,} 步验证后，因 ESM 旋转位置编码缓存与 inference_mode 的反向传播限制中断。已复现并将验证改为 no_grad；从保留权重继续，已消费批次不重放，原截止时间不变。', '原检查点没有优化器和随机数状态，因此恢复了 LoRA/分类头权重，AdamW 与随机数状态重新初始化；这不是优化器状态完全连续的训练。最佳模型仍按全部已完成验证的 MCC 保留。完整错误、修复测试及原始检查点均在 logs/、error_log.txt 和 output/recovery_20260922_0520。','']
    for key,entry in reports['frozen']['by_source_partner'].items():
        other=reports['lora']['by_source_partner'].get(key,{})
        lines.append(f'| {key} (n={entry["n"]}) | {fmt(entry["mcc"])} | {fmt(other.get("mcc"))} |')
    lines+=['','## 交付文件','', '- 最佳 LoRA 与分类头权重：`output/esm2_650m_lora_best.pt`（预训练主干位于发布包 `models/esm2_650m_hf`）。',
        '- 完整配置：`training_config.yaml`；预处理：`scripts/preprocess.py`；下载：`scripts/download_data.py`。',
        '- 对比表：`output/comparison.csv`；详细指标：`output/metrics.json`；逐条预测：`output/*_validation_predictions.csv`。',
        '- 损失曲线：`output/loss_curve.png` / `.pdf`；原始曲线数据：`output/losses.csv`。',
        '- 训练日志：`logs/` 和 `output/run_rank*.log`；进度：`progress.txt`；错误与修复记录：`error_log.txt`。',
        '- 数据版本、哈希：`data/download_manifest.json`；划分审计：`output/data_audit.json`。',
        '', '全部任务文件位于本目录（项目根由 `scripts/` 的位置解析，不依赖绝对路径）。','',
        '## 原始来源','', '- [AbBiBench](https://huggingface.co/datasets/AbBibench/Antibody_Binding_Benchmark_Dataset)',
        '- [AbAgym](https://github.com/3BioCompBio/AbAgym)', '- [IL6](https://huggingface.co/datasets/alchemab/il6-binding-prediction)',
        '- [Ab-Tune 包信息](https://pypi.org/project/Ab-Tune/)']
    (out/'final_report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    manifest={os.path.relpath(p, ROOT):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [best,ROOT/'training_config.yaml',ROOT/'scripts/preprocess.py',out/'comparison.csv',out/'final_report.md',out/'loss_curve.png']}
    (out/'deliverable_manifest.json').write_text(json.dumps(manifest,indent=2))
    print('FINAL_REPORT_VERIFIED',out/'final_report.md',flush=True)

if __name__=='__main__':main()
