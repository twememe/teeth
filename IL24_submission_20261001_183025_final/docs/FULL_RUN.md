# 完整上游计算和训练

本页命令属于长流程，不在顶层默认演示中执行。先按`envs/README.md`准备相应环境；首次下载需要联网。所有重算使用新工作副本，避免覆盖随包参考结果。`R`只是在当前shell读取的发布根目录，不是硬编码服务器路径。

## 阶段一：从序列/结构到证据、区段与热点

使用已安装conda。先按原规格建立核心与BepiPred环境（历史锁定依赖保留，安装未在本轮重做）：

```bash
R="$PWD"
conda env create -f "$R/modules/phase1ab/envs/il24-core.yml"
conda create -y -n il24-bp3 -c conda-forge python=3.10.20 pip
conda run -n il24-bp3 python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
conda run -n il24-bp3 python -m pip install -r "$R/envs/phase1_bepipred.txt"
conda create -y -n il24-dt3-phase1b -c conda-forge python=3.14.7 gcc_linux-64=16.1.0 gxx_linux-64=16.1.0
conda run -n il24-dt3-phase1b python -m pip install torch==2.10.0 --index-url https://download.pytorch.org/whl/cu128
conda run -n il24-dt3-phase1b python -m pip install -r "$R/envs/phase1_discotope.txt"
conda run -n il24-dt3-phase1b python -m pip install --no-deps "$R/models/third_party/DiscoTope-3.0"
python tools/get_models.py --group phase1 --destination "$R/external_models/phase1/hub/checkpoints"
export TORCH_HOME="$R/external_models/phase1"
mkdir -p runs
cp -a modules/phase1ab runs/phase1_full
cd "$R/runs/phase1_full"
conda run -n il24-core python scripts/sequence_qc.py --reference inputs/reference_181.fasta --pdb inputs/AF-Q925S4-F1-model_v6.pdb --uniprot inputs/uniprot_Q925S4.fasta --results-dir results/01_sequence_qc --report reports/step01_sequence_qc.md
conda run -n il24-core python scripts/structure_qc.py --reference inputs/reference_181.fasta --pdb inputs/AF-Q925S4-F1-model_v6.pdb --output-dir results/02_structure_qc --report reports/step02_structure_qc.md
conda run -n il24-core python scripts/surface_analysis.py --reference inputs/reference_181.fasta --pdb inputs/AF-Q925S4-F1-model_v6.pdb --output-dir results/04_surface --report reports/step04_surface.md --failure-log logs/surface_full.log
conda run -n il24-core python scripts/conservation.py --reference inputs/reference_181.fasta --uniprot-json results/06_conservation/raw/uniprot_il24_mammalia.json --output-dir results/06_conservation --report reports/step06_conservation.md --retrieval-date 2026-08-21
BP_PY=$(conda run -n il24-bp3 python -c 'import sys;print(sys.executable)')
conda run -n il24-core python scripts/run_bepipred.py --bp3-python "$BP_PY" --repo "$R/models/third_party/BepiPred-3.0" --input "$PWD/inputs/reference_181.fasta" --raw-dir "$PWD/results/03_bepipred/raw" --scores-csv results/03_bepipred/bepipred_residue_scores.csv --profile results/03_bepipred/bepipred_profile.png --report reports/step03_bepipred.md --log logs/bepipred_inference.log
DT="$R/models/third_party/DiscoTope-3.0"
conda run -n il24-dt3-phase1b /usr/bin/time -v python "$DT/discotope3/main.py" --pdb_or_zip_file "$PWD/inputs/AF-Q925S4-F1-model_v6.pdb" --struc_type alphafold --out_dir "$PWD/results/05b_discotope_restored/raw" --models_dir "$DT/models" --verbose 2 > logs/discotope_full.log 2>&1
nvidia-smi --query-gpu=timestamp,name,memory.used,memory.free,utilization.gpu --format=csv > logs/discotope_full_gpu.csv
conda run -n il24-core python scripts/phase1b_discotope_outputs.py --raw-csv results/05b_discotope_restored/raw/AF-Q925S4-F1-model_v6/AF-Q925S4-F1-model_v6_A_discotope3.csv --reference inputs/reference_181.fasta --pdb inputs/AF-Q925S4-F1-model_v6.pdb --repo "$DT" --esm-checkpoint "$TORCH_HOME/hub/checkpoints/esm_if1_gvp4_t16_142M_UR50.pt" --gpu-log logs/discotope_full_gpu.csv --inference-log logs/discotope_full.log --output-dir results/05b_discotope_restored --report reports/step05b_discotope_restored.md
conda run -n il24-core python scripts/integration_phase1b.py
conda run -n il24-core python scripts/candidates_phase1b.py
```

旧Phase1A冻结校验属于原完整入口；如果主动改变上游输入/模型字节，冻结校验会拒绝冒充旧结果。官方源的获取记录在`inputs/manifest.json`、`methods_registry.csv`和`results/06_conservation/raw/uniprot_query_log.json`；上面保守性命令使用已固定的真实数据库响应，避免重查询日期改变输入。

实际区段选择入口（在上述工作副本中）：

```bash
cd "$R/runs/phase1_full/highlevel_phase1a_construct_v2"
conda run -n il24-core python scripts/input_nterm.py
conda run -n il24-core python scripts/hydrophobic_boundary.py
conda run -n il24-core python scripts/structural_boundary.py
conda run -n il24-core python scripts/final_decision.py
```

`input_nterm.py`明确读取随包既有DeepSig真实GFF和UniProt注释，不声称重新推理DeepSig；其他三步执行实际疏水/结构边界/决策计算。DeepSig原环境Python3.8.20、deepsig0.9、TF2.2.0、Keras2.4.3、numpy1.18.5、biopython1.78、protobuf3.20.3。原模型目录已随包放在`models/third_party/DeepSig`，SHA见`models/deepsig_assets.json`。官方来源https://github.com/BolognaBiocomp/deepsig；历史模型revision未记录，按随包原字节固定。实际上游重推理命令如下，运行后再执行上述input_nterm.py读取新GFF：

```bash
conda create -y -n il24-deepsig38 -c conda-forge python=3.8.20 pip
conda run -n il24-deepsig38 python -m pip install deepsig-biocomp==0.9 tensorflow==2.2.0 keras==2.4.3 numpy==1.18.5 biopython==1.78 protobuf==3.20.3
export DEEPSIG_ROOT="$R/models/third_party/DeepSig"
conda run -n il24-deepsig38 deepsig -f "$PWD/inputs/reference_181.fasta" -o "$PWD/results/01_n_terminal/deepsig_raw/deepsig_il24.gff3" -k euk
```

该参数与[官方使用说明](https://raw.githubusercontent.com/BolognaBiocomp/deepsig/master/README.md)一致。SignalP6没有被使用，未提供受限权重。本版不把缺失授权的软件补成假模型。

## 阶段二：正式预测→接触→PRODIGY→排序

激活Phase2环境，从发布根目录：

```bash
R="$PWD"
python tools/get_models.py --group boltz --destination "$R/external_models/boltz"
python tools/get_models.py --group antibody --destination "$R/external_models/antibody"
export BOLTZ_CACHE="$R/external_models/boltz"
mkdir -p runs/phase2_full
cp -a modules/phase2/scripts modules/phase2/data modules/phase2/vendor modules/phase2/hmmer_bin runs/phase2_full/
mkdir -p runs/phase2_full/results/12_complex_prediction
cp -a modules/phase2/results/12_complex_prediction/boltz_inputs runs/phase2_full/results/12_complex_prediction/
cp modules/phase2/results/12_complex_prediction/phase1_prior_deduplicated.csv runs/phase2_full/results/12_complex_prediction/
cd "$R/runs/phase2_full"
export PATH="$PWD/hmmer_bin/bin:$PATH"
python scripts/run_imgt_numbering.py
mkdir -p results/11_antibody_structure
python -c 'from ImmuneBuilder import ABodyBuilder2;from pathlib import Path; import sys;sys.path.insert(0,"scripts");from step11_antibody_structure import parse_fasta,identify_chain;_,h,_,l=identify_chain(parse_fasta(Path("data/VH_VL_clean.fasta")));m=ABodyBuilder2(weights_dir=Path("../../external_models/antibody"));a=m.predict({"H":h,"L":l});a.save_all("results/11_antibody_structure")'
# 当前正式输入已含固定MSA。若有意从远程服务重新获取，下面是真实入口：
# python scripts/generate_colabfold_msa.py
bash scripts/run_boltz_formal.sh
bash scripts/run_boltz_formal_msa.sh
python scripts/analyze_boltz_complexes.py --set formal
python scripts/analyze_boltz_complexes.py --set formal_msa
python scripts/run_prodigy_formal.py --set formal
python scripts/run_prodigy_formal.py --set formal_msa
python scripts/rank_formal_msa_models.py
python scripts/run_phase2_external_analysis.py --orthologs "$R/modules/phase1ab/results/06_conservation/orthologs.fasta"
cd "$R"
python run.py --module export --input runs/phase2_full/results/14_quality/model_ranking_formal_msa.csv --source-root runs/phase2_full --output runs/phase2_full_export
```

导出自定义ranking时必须用`--source-root`指定同一工作副本；序列、指标和PDB从该副本共同读取。未指定且ranking在默认模块之外时入口拒绝混用。糖基化分支保留为探索性敏感性。

Full181独立完整入口（工作副本，不使用已有run_state跳过新预测）：

```bash
cd "$R"
mkdir -p runs/full181_full/results/full181_supplement
cp -a modules/phase2/full181_supplement/src modules/phase2/full181_supplement/data runs/full181_full/
cp -a modules/phase2/full181_supplement/results/full181_supplement/inputs modules/phase2/full181_supplement/results/full181_supplement/msa modules/phase2/full181_supplement/results/full181_supplement/runtime_sources runs/full181_full/results/full181_supplement/
cd runs/full181_full
python src/prepare_full181_supplement.py
python src/run_full181_supplement.py --condition empty --gpu 0 --boltz "$(command -v boltz)" --checkpoint "$BOLTZ_CACHE/boltz1_conf.ckpt" --cache "$BOLTZ_CACHE"
python src/run_full181_supplement.py --condition formal_msa --gpu 1 --boltz "$(command -v boltz)" --checkpoint "$BOLTZ_CACHE/boltz1_conf.ckpt" --cache "$BOLTZ_CACHE"
python src/analyze_full181_supplement.py
python src/run_full181_supplement.py --prodigy "$(command -v prodigy)"
```

其比较/代表选择规则在`full181_ensemble.py`，历史完整比较数据与baseline输入仍保留在模块中；不将Full181不同条件强塞进原主清单。

## ESM-2：固定数据→提取→训练→加载评估→报告

先激活ESM环境，从发布根目录准备冻结官方基座（仅转换，不训练）：

```bash
python tools/prepare_base.py
```

配置默认使用`external_models/esm2_650m_hf`；离线官方文件或已转换目录用法见BASE_MODEL.md。

本轮交付的固定划分直接可用。若重新取得并预处理原始公开数据，只在副本中执行（下载锁定到原revision）：

```bash
R="$PWD"
mkdir -p runs
cp -a modules/inovate runs/inovate_preprocess
python runs/inovate_preprocess/scripts/download_data.py
python runs/inovate_preprocess/scripts/preprocess.py
```

原标签及分割不可为提高结果而改变。许可与缺口见DATA_SOURCES。

激活ESM环境。从发布根目录进行从头完整训练（显式执行才训练，默认run.sh不进入此流程）。以下按原超参数进行全量训练，不能保证原中断恢复过程的逐位轨迹：

```bash
R="$PWD"
mkdir -p runs/esm2_training
cp modules/inovate/output/data_audit.json runs/esm2_training/
python -c 'from datetime import datetime,timedelta,timezone;import json;from pathlib import Path;t=datetime.now(timezone.utc);Path("runs/esm2_training/run_state.json").write_text(json.dumps({"experiment_started_at":t.isoformat(),"hard_deadline_at":(t+timedelta(hours=13)).isoformat()}))'
torchrun --standalone --nproc_per_node=2 modules/inovate/scripts/train.py --mode extract --config configs/esm2_train.json --max-seconds 10800
python modules/inovate/scripts/train.py --mode baseline --config configs/esm2_train.json --max-seconds 3600
torchrun --standalone --nproc_per_node=2 modules/inovate/scripts/train.py --mode lora --config configs/esm2_train.json --max-seconds 43200
python modules/inovate/scripts/train.py --mode evaluate --config configs/esm2_train.json --checkpoint runs/esm2_training/esm2_650m_frozen_best.pt --split validation --max-seconds 3600
python modules/inovate/scripts/train.py --mode evaluate --config configs/esm2_train.json --checkpoint runs/esm2_training/esm2_650m_lora_best.pt --split validation --max-seconds 3600
python modules/inovate/scripts/report.py --output-dir runs/esm2_training --run-state runs/esm2_training/run_state.json
```

训练batch16/GPU×2，LoRA学习率1e-4、seed42、bf16，最多15轮、每2000步验证、epoch patience3，最大验证MCC选点，阈值0.5。原基线有独立optimizer日程（学习率0.001、batch4096、最多30轮、patience3）。报告默认原500次分组bootstrap只是已有报告流程；本发布复核不运行它。历史原记录的重报告可加`--reuse-uncertainty`读取保存区间，不新增统计。

完整评估最终保存权重，不重新训练：

```bash
python modules/inovate/scripts/train.py --mode evaluate --config configs/esm2_train.json --checkpoint modules/inovate/output/esm2_650m_lora_best.pt --split validation --max-seconds 3600
# 测试集仅在用户明确开展独立测试时使用：同入口改 --split test；本轮未执行
```

此配置输出到新`runs/esm2_training`。少量推理使用README的`run.py --module esm2`，两者目的不同。
