# 模型资产

ESM官方基座已外置，按`base_manifest.json`固定来源下载并由`../tools/prepare_base.py`转换/校验。`base_reference.json`来自本轮重建前独立加载的原正式基座；`base_config`保留原config、vocab和tokenizer。项目训练的LoRA与分类头仍完整保存在`../modules/inovate/output/esm2_650m_lora_best.pt`，epoch3/step14000、阈值0.5，checkpoint字节不变；冻结基线head及恢复checkpoint也保留。默认准备位置`../external_models/esm2_650m_hf`，该目录及下载缓存不放入ZIP。实际加载比较张量/配置/tokenizer指纹，不把官方源SHA或旧序列化SHA当成转换模型身份。命令见`../docs/BASE_MODEL.md`。

第三方完整上游大模型不必重复放入包：Boltz `boltz1_conf.ckpt`、`ccd.pkl`以及ABodyBuilder2四个模型的大小、SHA和来源在`external_assets.json`。下载工具从该表读取真实官方URL并验证；已有文件只做SHA核对，不重复下载：

```bash
python tools/get_models.py --group boltz --destination ./external_models/boltz
python tools/get_models.py --group antibody --destination ./external_models/antibody
```

BepiPred与DiscoTope的小模型/仓库代码在`third_party`；其大ESM/IF1官方URL和SHA来自原运行记录，见`phase1_assets.json`，可用同一工具`--group phase1`下载到`external_models/phase1/hub/checkpoints`。完整环境和调用命令在`docs/FULL_RUN.md`。原始下载日期和缺失元数据如实记录；本轮没有重新下载这些大权重。
